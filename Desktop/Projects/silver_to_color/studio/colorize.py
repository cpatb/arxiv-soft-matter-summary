"""Local restoration: decimate + deflicker a B&W master, then colorize it on this laptop.

Two outputs, frame-aligned one to one, which is what makes every before/after effect
downstream (split screen, wipe, Shorts stack) trivially exact:

    bw_master.mp4      decimated to the film's true rate, deflickered, working resolution
    color_master.mp4   the same frames, colorized

WHY OPENVINO
------------
This machine has no CUDA GPU. DDColor-tiny (512x512) measured on one real frame:

    ONNX Runtime CPU, fp16 weights      ~12.0 s/frame  (fp16 is emulated on CPU)
    ONNX Runtime CPU, fp32 weights       ~8.5 s/frame
    OpenVINO CPU, fp32                   ~4.6 s/frame
    OpenVINO Arc iGPU, f16 inference     ~0.44 s/frame

so the integrated GPU is ~20x faster than any CPU path and a 12-minute silent film
(~11,500 unique frames) colorizes in about 1.5 hours. ONNX Runtime stays as a fallback.

The model is the edgetools export of DDColor-tiny (Apache-2.0, upstream piddnad/DDColor)
with its float16 weights cast to float32 by `convert_fp16_model()`, because OpenVINO's
GPU plugin picks its own inference precision and the fp16 graph's Cast nodes only cost.

WHY OUR OWN DECIMATION
----------------------
The Market Street 2K file is 24 fps with every second frame doubled, but it is an h264
re-encode, so the duplicates are not bit-identical: consecutive-frame mean |diff| runs
~1.25 grey levels for a duplicate against 5-11 for real motion (measured, 48 frames at
t=100 s). ffmpeg's mpdecimate at its default thresholds kept 47 of 48 frames. Here a
frame is a duplicate when its difference from the previous frame is below
`DUP_RATIO` x the local median difference, which adapts to grain and motion level.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import cv2

from . import _TOOLS  # noqa: F401  (puts ../tools on sys.path)
import temporal_chroma as tc
import pipeline as pl

MODEL_DIR = Path(os.environ.get("SILVER_MODELS", r"E:\silver_to_color\models"))
DEFAULT_MODEL = MODEL_DIR / "ddcolor-tiny-fp32.onnx"
FP16_SOURCE = MODEL_DIR / "ddcolor-tiny-fp16.onnx"
FP16_URL = "https://huggingface.co/edgetools/ddcolor/resolve/main/ddcolor-tiny-fp16.onnx"

INPUT_SIZE = 512            # fixed by the export
DUP_RATIO = 0.35            # duplicate if diff < DUP_RATIO * local median diff
DUP_WINDOW = 15             # frames in the rolling median for that test
PROBE_SIZE = (256, 192)     # resolution of the duplicate-detection pass


# ----------------------------------------------------------------------------- model

def convert_fp16_model(src: Path = FP16_SOURCE, dst: Path = DEFAULT_MODEL) -> Path:
    """Cast an fp16-weight ONNX graph to fp32 (initializers, constants, casts, value_info).

    Checked against the fp16 original through ONNX Runtime on a random input: max |diff|
    0.21 on an output range of +-301, mean 0.009, i.e. fp16 rounding and nothing else.
    """
    import onnx
    from onnx import numpy_helper, TensorProto

    if not src.exists():
        raise FileNotFoundError(f"{src} missing; download it from {FP16_URL}")
    model = onnx.load(str(src))
    f16, f32 = TensorProto.FLOAT16, TensorProto.FLOAT

    def fix(tensor):
        if tensor.data_type == f16:
            arr = numpy_helper.to_array(tensor).astype(np.float32)
            tensor.CopyFrom(numpy_helper.from_array(arr, tensor.name))

    for init in model.graph.initializer:
        fix(init)
    for node in model.graph.node:
        for attr in node.attribute:
            if node.op_type == "Cast" and attr.name == "to" and attr.i == f16:
                attr.i = f32
            if attr.type == onnx.AttributeProto.TENSOR:
                fix(attr.t)
    for vi in list(model.graph.value_info) + list(model.graph.input) + list(model.graph.output):
        if vi.type.tensor_type.elem_type == f16:
            vi.type.tensor_type.elem_type = f32
    onnx.save(model, str(dst))
    return dst


class Colorizer:
    """DDColor-tiny chroma predictor. Input: grey uint8 frames. Output: (512,512,2) ab."""

    def __init__(self, model_path: Path = DEFAULT_MODEL, device: str = "GPU",
                 saturation: float = 1.0, parallel: int = 4):
        model_path = Path(model_path)
        if not model_path.exists():
            convert_fp16_model(dst=model_path)
        self.saturation = float(saturation)
        self.backend = None
        try:
            import openvino as ov
            core = ov.Core()
            if device not in core.available_devices:
                device = "CPU"
            hint = "f16" if device == "GPU" else "f32"
            self._compiled = core.compile_model(
                core.read_model(str(model_path)), device,
                {"INFERENCE_PRECISION_HINT": hint, "PERFORMANCE_HINT": "THROUGHPUT"})
            self._parallel = max(1, int(parallel))
            self.backend = f"openvino:{device}:{hint}"
        except Exception as exc:  # pragma: no cover - depends on the machine
            import onnxruntime as ort
            opts = ort.SessionOptions()
            opts.log_severity_level = 3
            self._session = ort.InferenceSession(str(model_path), opts,
                                                 providers=["CPUExecutionProvider"])
            self.backend = f"onnxruntime:cpu (openvino unavailable: {exc})"

    @staticmethod
    def model_input(gray_u8: np.ndarray) -> np.ndarray:
        """Grey frame -> (1,3,512,512) float32 RGB rebuilt from L with a=b=0."""
        small = cv2.resize(gray_u8, (INPUT_SIZE, INPUT_SIZE), interpolation=cv2.INTER_AREA)
        bgr = cv2.cvtColor(small, cv2.COLOR_GRAY2BGR).astype(np.float32) / 255.0
        lab = cv2.cvtColor(bgr, cv2.COLOR_BGR2Lab)
        lab[..., 1:] = 0.0
        rgb = cv2.cvtColor(lab, cv2.COLOR_LAB2RGB)
        return np.ascontiguousarray(rgb.transpose(2, 0, 1)[None], dtype=np.float32)

    def submit(self, frames_gray: list[np.ndarray]):
        """Start inference without blocking (OpenVINO); call the returned function to collect.

        The point is overlap: the iGPU infers while the CPU computes optical flow.
        """
        inputs = [self.model_input(f) for f in frames_gray]
        outs: list[np.ndarray | None] = [None] * len(inputs)
        if self.backend.startswith("openvino"):
            import openvino as ov
            queue = ov.AsyncInferQueue(self._compiled, self._parallel)

            def done(request, index):
                outs[index] = request.get_output_tensor(0).data[0].copy()
            queue.set_callback(done)
            for i, x in enumerate(inputs):
                queue.start_async({0: x}, i)

            def collect():
                queue.wait_all()
                return self._finish(outs)
            return collect
        for i, x in enumerate(inputs):
            outs[i] = self._session.run(None, {"input": x})[0][0]
        return lambda: self._finish(outs)

    def predict_ab(self, frames_gray: list[np.ndarray]) -> list[np.ndarray]:
        return self.submit(frames_gray)()

    def _finish(self, outs) -> list[np.ndarray]:
        if any(o is None for o in outs):
            raise RuntimeError("inference returned fewer outputs than inputs")
        result = []
        for o in outs:
            ab = np.ascontiguousarray(o.transpose(1, 2, 0), dtype=np.float32)
            if self.saturation != 1.0:
                ab *= self.saturation
            result.append(ab)
        return result


def recombine(gray_full: np.ndarray, ab_native: np.ndarray) -> np.ndarray:
    """Full-resolution L from the B&W frame + bilinear-upsampled ab -> BGR uint8."""
    h, w = gray_full.shape
    bgr = cv2.cvtColor(gray_full, cv2.COLOR_GRAY2BGR).astype(np.float32) / 255.0
    lab = cv2.cvtColor(bgr, cv2.COLOR_BGR2Lab)
    lab[..., 1:] = tc.upsample_chroma(ab_native, h, w)
    out = cv2.cvtColor(lab, cv2.COLOR_LAB2BGR)
    return np.clip(out * 255.0 + 0.5, 0, 255).astype(np.uint8)


# ----------------------------------------------------------------------------- io

def read_gray(path, width: int, height: int, start_s: float | None = None,
              frames: int | None = None, duration_s: float | None = None):
    """Yield grey uint8 frames (height, width) from any video ffmpeg can read."""
    cmd = ["ffmpeg", "-v", "error"]
    if start_s is not None:
        cmd += ["-ss", f"{start_s:.6f}"]
    cmd += ["-i", str(path)]
    if duration_s is not None:
        cmd += ["-t", f"{duration_s:.6f}"]
    if frames is not None:
        cmd += ["-frames:v", str(frames)]
    cmd += ["-vf", f"scale={width}:{height}:flags=lanczos,format=gray",
            "-f", "rawvideo", "-pix_fmt", "gray", "-"]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    size = width * height
    try:
        while True:
            buf = proc.stdout.read(size)
            if len(buf) < size:
                break
            yield np.frombuffer(buf, np.uint8).reshape(height, width)
    finally:
        proc.stdout.close()
        proc.wait()


# YUV matrix for written video. swscale's RGB->YUV default is BT.601 and libx264 leaves
# the stream untagged, while YouTube/HD players assume BT.709: reds and skin shift
# (review 2026-10-01: pure red 200 -> Y 67, the 601 value). "bt601" exists only so a
# run started under the old default can be resumed consistently and tagged honestly.
MATRIX = "bt709"
TAGS = {"bt709": ("bt709", "bt709", "bt709"), "bt601": ("smpte170m", "bt709", "bt709")}
H273 = {"bt709": 1, "bt601": 6}


def color_tags(matrix: str) -> list[str]:
    cs, prim, trc = TAGS[matrix]
    return ["-colorspace", cs, "-color_primaries", prim, "-color_trc", trc]


def writer(path, width: int, height: int, fps: float, pix_fmt_in: str,
           vf: str | None = None, crf: int = 12, matrix: str = MATRIX):
    cmd = ["ffmpeg", "-v", "error", "-y", "-f", "rawvideo", "-pix_fmt", pix_fmt_in,
           "-s", f"{width}x{height}", "-r", f"{fps:.6f}", "-i", "-"]
    conv = f"scale=out_color_matrix={matrix}:out_range=tv,format=yuv420p"
    cmd += ["-vf", f"{vf},{conv}" if vf else conv]
    cmd += ["-c:v", "libx264", "-crf", str(crf), "-preset", "medium",
            "-pix_fmt", "yuv420p", "-g", str(int(round(fps)))] + color_tags(matrix) + [
            "-movflags", "+faststart", str(path)]
    return subprocess.Popen(cmd, stdin=subprocess.PIPE)


def duplicate_mask(diffs: np.ndarray, ratio: float = DUP_RATIO,
                   window: int = DUP_WINDOW) -> np.ndarray:
    """diffs[i] = mean |frame[i+1] - frame[i]|. Returns mask over frames 1..N-1 (dup = True).

    A duplicate is a frame much closer to its predecessor than the local norm. The local
    norm is a centred rolling median, which is dominated by real motion as long as fewer
    than half the frames in the window are duplicates (1/3 for the 16->24 recipe).
    """
    diffs = np.asarray(diffs, dtype=np.float64)
    n = len(diffs)
    half = window // 2
    mask = np.zeros(n, dtype=bool)
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        local = np.median(diffs[lo:hi])
        mask[i] = diffs[i] < ratio * local
    return mask


PERIOD = 3                  # 16 -> 24 fps pulldown: one frame in three is a repeat
PHASE_WINDOW = 48           # +- diffs used to decide the local pulldown phase
PHASE_CONTRAST = 0.5        # a phase counts only if its median diff < this x the others'


def pulldown_mask(diffs: np.ndarray, period: int = PERIOD, window: int = PHASE_WINDOW,
                  contrast: float = PHASE_CONTRAST) -> np.ndarray | None:
    """Phase-locked duplicate detection for a periodic pulldown. None if no pulldown.

    diffs[i] = mean |frame[i+1] - frame[i]|; returns a mask over frames 1..N-1.

    The scanner repeated frames on a fixed cycle, so the right question is not "is this
    frame similar to the last?" (which fails on static shots, where every diff is small)
    but "which phase of the 3-cycle is the repeat?". Measured on Market Street: diffs at
    one phase are 99.1 % below 3 grey levels, the other two phases ~5 %. The phase is
    estimated in a sliding window, which follows the phase SHIFT the film has at a splice
    near t = 150 s, and is carried into windows with no motion to judge by (the black
    leader) from the nearest window that has.
    """
    diffs = np.asarray(diffs, dtype=np.float64)
    n = len(diffs)
    phase = np.arange(n) % period
    best = np.full(n, -1)
    for i in range(n):
        lo, hi = max(0, i - window), min(n, i + window + 1)
        seg, ph = diffs[lo:hi], phase[lo:hi]
        meds = np.array([np.median(seg[ph == p]) if np.any(ph == p) else np.inf
                         for p in range(period)])
        p = int(np.argmin(meds))
        rest = np.delete(meds, p)
        if np.isfinite(rest).all() and meds[p] < contrast * np.median(rest):
            best[i] = p
    known = np.nonzero(best >= 0)[0]
    if len(known) < 0.25 * n:           # no consistent pulldown anywhere: not this recipe
        return None
    # Carry each confident phase into neighbouring unconfident stretches.
    idx = np.searchsorted(known, np.arange(n))
    left = known[np.clip(idx - 1, 0, len(known) - 1)]
    right = known[np.clip(idx, 0, len(known) - 1)]
    nearest = np.where(np.abs(np.arange(n) - left) <= np.abs(right - np.arange(n)), left, right)
    filled = np.where(best >= 0, best, best[nearest])
    return phase == filled


@dataclass
class PrepareReport:
    source: str
    start_s: float | None
    duration_s: float
    source_fps: float
    frames_in: int
    duplicates: int
    frames_kept: int
    fps_out: float
    width: int
    height: int
    method: str
    dup_ratio: float
    deflicker: str
    seconds: float


def frame_diffs(source, start_s=None, duration_s=None) -> np.ndarray:
    """Mean |diff| between consecutive frames at PROBE_SIZE, streamed (no film in RAM)."""
    diffs, prev = [], None
    for frame in read_gray(source, *PROBE_SIZE, start_s=start_s, duration_s=duration_s):
        cur = frame.astype(np.float32)
        if prev is not None:
            diffs.append(float(np.abs(cur - prev).mean()))
        prev = cur
    return np.asarray(diffs)


def prepare(source, out, height: int = 1080, start_s: float | None = None,
            duration_s: float | None = None, deflicker: bool = True,
            decimate: bool = True, diffs_cache=None) -> PrepareReport:
    """Source -> B&W master at the film's true rate. Writes `<out>.json` beside it.

    The diff trace and the keep mask are saved to `<out>.diffs.npz` for audit; pass it
    back as `diffs_cache` to skip pass 1 (only valid for the same source/start/duration).
    """
    t0 = time.time()
    info = pl.probe_video(str(source))
    src_fps = info.avg_frame_rate
    total = info.duration if duration_s is None else duration_s
    width = int(round(height * info.width / info.height / 2) * 2)

    # Pass 1: duplicate detection at low resolution.
    if diffs_cache and Path(diffs_cache).exists():
        diffs = np.load(diffs_cache)["diffs"]
    else:
        diffs = frame_diffs(source, start_s, duration_s)
    frames_in = len(diffs) + 1
    if frames_in < 2:
        raise RuntimeError(f"read {frames_in} frames from {source}")
    keep = np.ones(frames_in, dtype=bool)
    method = "none"
    if decimate:
        mask = pulldown_mask(diffs)
        method = "pulldown-phase-locked"
        if mask is None:
            mask, method = duplicate_mask(diffs), "adaptive-threshold"
        keep[1:] = ~mask
    kept = int(keep.sum())
    dropped_moving = int(((~keep[1:]) & (diffs > 3.0)).sum())
    print(f"[prepare] {method}: dropped {frames_in - kept} of {frames_in}; "
          f"{dropped_moving} dropped frames differ from their predecessor by > 3 grey levels",
          flush=True)
    np.savez(str(out) + ".diffs.npz", diffs=diffs, keep=keep)

    # True rate: what the kept frames must play at to fill the same running time. Snap to
    # a nominal rate when within 2 %, otherwise keep the measured value and say so.
    measured = kept / (frames_in / src_fps)
    nominal = min((12.0, 16.0, 18.0, 20.0, 24.0), key=lambda r: abs(r - measured))
    fps_out = nominal if abs(nominal - measured) / nominal < 0.02 else measured
    print(f"[prepare] frames in {frames_in}, duplicates {frames_in - kept}, kept {kept} "
          f"-> {measured:.3f} fps measured, writing at {fps_out:g} fps", flush=True)

    # Pass 2: write kept frames at working resolution.
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    vf = "deflicker=size=7:mode=am" if deflicker else None
    proc = writer(out, width, height, fps_out, "gray", vf=vf)
    written = 0
    for i, frame in enumerate(read_gray(source, width, height, start_s=start_s,
                                        duration_s=duration_s)):
        if i < frames_in and keep[i]:
            proc.stdin.write(frame.tobytes())
            written += 1
    proc.stdin.close()
    if proc.wait() != 0:
        raise RuntimeError("ffmpeg failed writing the B&W master")
    if written != kept:
        raise RuntimeError(f"pass 2 wrote {written} frames, pass 1 kept {kept}")

    report = PrepareReport(str(source), start_s, float(total), float(src_fps), frames_in,
                           frames_in - kept, kept, float(fps_out), width, height, method,
                           DUP_RATIO,
                           "deflicker=size=7:mode=am" if deflicker else "none",
                           round(time.time() - t0, 1))
    Path(str(out) + ".json").write_text(json.dumps(asdict(report), indent=2), "utf-8")
    return report


# ----------------------------------------------------------------------------- temporal

FLOW_SIZE = 256   # flow is estimated at 256x256 and scaled up to the 512x512 chroma grid


def pair_flows(small_gray: list[np.ndarray], half: int, threads: int = 8) -> dict:
    """Optical flow for every ordered pair (i, j) with 0 < |i-j| <= half.

    Each pair is computed once per chunk instead of once per centre frame (a 2x saving
    over calling tc.flow_warped_median frame by frame), at FLOW_SIZE, on a thread pool:
    OpenCV releases the GIL, so this runs on the CPU while the iGPU does inference.
    Vectors are rescaled to the 512 grid, so the forward-backward tolerance stays in
    512-grid pixels as in tc.DEFAULT_FB_THRESHOLD_PX.
    """
    from concurrent.futures import ThreadPoolExecutor
    import threading

    scale = INPUT_SIZE / FLOW_SIZE
    tiny = [cv2.resize(g, (FLOW_SIZE, FLOW_SIZE), interpolation=cv2.INTER_AREA)
            for g in small_gray]
    local = threading.local()

    def one(pair):
        if not hasattr(local, "dis"):
            local.dis = tc.make_dis_flow("medium")
        i, j = pair
        f = local.dis.calc(tiny[i], tiny[j], None)
        f = cv2.resize(f, (INPUT_SIZE, INPUT_SIZE), interpolation=cv2.INTER_LINEAR) * scale
        return pair, f.astype(np.float32)

    n = len(small_gray)
    pairs = [(i, j) for i in range(n) for j in range(max(0, i - half), min(n, i + half + 1))
             if i != j]
    with ThreadPoolExecutor(threads) as pool:
        return dict(pool.map(one, pairs))


def temporal_median(ab: list[np.ndarray], flows: dict, i: int, half: int,
                    fb_px: float = tc.DEFAULT_FB_THRESHOLD_PX) -> np.ndarray:
    """Same estimator as tc.flow_warped_median, reading precomputed flows.

    Neighbours are warped into frame i by flow i->j, rejected where the round trip
    i->j->i misses by more than fb_px, and the per-pixel median is taken over survivors
    plus frame i itself (always valid, so no pixel can be empty).
    """
    n = len(ab)
    samples, valid = [ab[i]], [np.ones(ab[i].shape[:2], dtype=bool)]
    for j in range(max(0, i - half), min(n, i + half + 1)):
        if j == i:
            continue
        fw, bw = flows[(i, j)], flows[(j, i)]
        samples.append(tc.warp_by_flow(ab[j], fw))
        valid.append(tc.fb_consistency_mask(fw, bw, fb_px))
    return masked_median(np.stack(samples), np.stack(valid))


def masked_median(stack: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Per-pixel median of stack (N,H,W,C) over entries where mask (N,H,W) is True.

    Equal to np.nanmedian on the NaN-masked stack (tested), ~30x faster: np.sort puts
    NaN last, so with k valid samples the median is the mean of sorted[(k-1)//2] and
    sorted[k//2].
    """
    data = np.where(mask[..., None], stack, np.nan).astype(np.float32)
    data.sort(axis=0)
    k = np.broadcast_to(mask.sum(axis=0)[..., None], data.shape[1:])
    lo = np.take_along_axis(data, ((k - 1) // 2)[None], axis=0)[0]
    hi = np.take_along_axis(data, (k // 2)[None], axis=0)[0]
    return ((lo + hi) * 0.5).astype(np.float32)


# ----------------------------------------------------------------------------- colorize

def count_packets(path) -> int:
    """Frame count from the container's packets: no decoding, so seconds instead of the
    ~20 min `-count_frames` spent decoding a 1.7 GB grainy master. Exact for the
    intra/inter-coded streams written here (one packet per frame, no B-frame reorder
    changes the count); test_chunked_colorize_equals_single_chunk cross-checks it."""
    out = subprocess.run(["ffprobe", "-v", "error", "-select_streams", "v:0", "-count_packets",
                          "-show_entries", "stream=nb_read_packets", "-of", "csv=p=0",
                          str(path)], capture_output=True, text=True, check=True).stdout
    return int(out.strip().split(",")[0])


def _chunk_name(workdir: Path, index: int) -> Path:
    return workdir / f"chunk_{index:04d}.mp4"


def colorize(bw_master, out, colorizer: Colorizer | None = None, window: int = 5,
             chunk_frames: int = 160, workdir=None, crf: int = 14) -> dict:
    """B&W master -> colour master, chunked and resumable.

    Each chunk is computed with a halo of window//2 frames on each side so the temporal
    median is identical whether a film is done in one session or ten. A finished chunk is
    recorded in `manifest.json` only after its frame count has been verified.
    """
    bw_master, out = Path(bw_master), Path(out)
    workdir = Path(workdir) if workdir else out.with_suffix("").with_name(out.stem + "_chunks")
    workdir.mkdir(parents=True, exist_ok=True)
    info = pl.probe_video(str(bw_master))
    fps, w, h = info.avg_frame_rate, info.width, info.height
    n = count_packets(bw_master)
    half = window // 2
    manifest_path = workdir / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    if manifest.get("source_frames") not in (None, n):
        raise RuntimeError(f"manifest is for {manifest['source_frames']} frames, master has {n}")
    colorizer = colorizer or Colorizer()
    if manifest.get("done"):
        # Resuming: every setting that changes pixels must match the finished chunks, or
        # the master silently splices two different restorations. Manifests from before
        # these keys existed were written as BT.601.
        manifest.setdefault("matrix", "bt601")
        want = {"window": window, "chunk_frames": chunk_frames,
                "saturation": colorizer.saturation, "backend": colorizer.backend}
        for key, val in want.items():
            if key in manifest and manifest[key] != val:
                raise RuntimeError(f"resume mismatch: {key} was {manifest[key]!r}, now {val!r}")
    matrix = manifest.get("matrix", MATRIX)
    manifest.update(source=str(bw_master), source_frames=n, fps=fps, window=window,
                    chunk_frames=chunk_frames, saturation=colorizer.saturation,
                    backend=colorizer.backend, matrix=matrix)
    manifest.setdefault("done", {})
    n_chunks = (n + chunk_frames - 1) // chunk_frames
    t_start = time.time()
    done_now = 0

    for c in range(n_chunks):
        s, e = c * chunk_frames, min(n, (c + 1) * chunk_frames)
        if str(c) in manifest["done"] and _chunk_name(workdir, c).exists():
            if manifest["done"][str(c)].get("first") != s:
                raise RuntimeError(f"chunk {c} was made with a different chunk layout")
            continue
        lo, hi = max(0, s - half), min(n, e + half)
        # Seek a quarter frame early so rounding can never skip the first wanted frame.
        start = max(0.0, (lo - 0.25) / fps) if lo > 0 else None
        frames = list(read_gray(bw_master, w, h, start_s=start, frames=hi - lo))
        if len(frames) != hi - lo:
            raise RuntimeError(f"chunk {c}: read {len(frames)} frames, wanted {hi - lo}")
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(1) as gpu:
            pending = gpu.submit(colorizer.predict_ab, frames)   # iGPU starts here ...
            small = [cv2.resize(f, (INPUT_SIZE, INPUT_SIZE), interpolation=cv2.INTER_AREA)
                     for f in frames]
            flows = pair_flows(small, half)                     # ... CPU does flow meanwhile
            ab = pending.result()
        proc = writer(_chunk_name(workdir, c), w, h, fps, "bgr24", crf=crf, matrix=matrix)

        def finish(k):
            i = k - lo
            return recombine(frames[i], temporal_median(ab, flows, i, half))

        with ThreadPoolExecutor(6) as pool:
            for frame in pool.map(finish, range(s, e)):   # map preserves order
                proc.stdin.write(frame.tobytes())
        proc.stdin.close()
        if proc.wait() != 0:
            raise RuntimeError(f"ffmpeg failed on chunk {c}")
        got = count_packets(_chunk_name(workdir, c))
        if got != e - s:
            raise RuntimeError(f"chunk {c}: wrote {got} frames, expected {e - s}")
        manifest["done"][str(c)] = {"frames": got, "first": s}
        manifest_path.write_text(json.dumps(manifest, indent=2))
        done_now += e - s
        rate = (time.time() - t_start) / max(1, done_now)
        left = sum(min(n, (k + 1) * chunk_frames) - k * chunk_frames
                   for k in range(n_chunks) if str(k) not in manifest["done"])
        print(f"[colorize] chunk {c + 1}/{n_chunks} ok ({got} frames), "
              f"{rate:.2f} s/frame, ~{left * rate / 60:.0f} min left", flush=True)

    chunks = [_chunk_name(workdir, c) for c in range(n_chunks)]
    listing = workdir / "concat.txt"
    listing.write_text("".join(f"file '{p.resolve().as_posix()}'\n" for p in chunks), "utf-8")
    subprocess.run(["ffmpeg", "-v", "error", "-y", "-f", "concat", "-safe", "0",
                    "-i", str(listing), "-c", "copy",
                    # Stream copy ignores -colorspace; write the VUI in the bitstream.
                    # H.273 codes: matrix 1 = BT.709, 6 = SMPTE 170M (601); prim/trc 1 = 709.
                    "-bsf:v", f"h264_metadata=matrix_coefficients={H273[matrix]}:"
                              "colour_primaries=1:transfer_characteristics=1",
                    "-movflags", "+faststart", str(out)],
                   check=True)
    total = count_packets(out)
    if total != n:
        raise RuntimeError(f"colour master has {total} frames, B&W master has {n}")
    manifest["output"] = str(out)
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"[colorize] {out} : {total} frames, aligned with {bw_master}", flush=True)
    return manifest


def still(image_path, out_path, colorizer: Colorizer | None = None) -> Path:
    """Colorize one image (for thumbnails and quick looks)."""
    img = cv2.imread(str(image_path), cv2.IMREAD_GRAYSCALE)
    if img is None:
        raise FileNotFoundError(image_path)
    colorizer = colorizer or Colorizer()
    cv2.imwrite(str(out_path), recombine(img, colorizer.predict_ab([img])[0]))
    return Path(out_path)


def main(argv=None) -> int:
    import argparse
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("prepare", help="source -> decimated, deflickered B&W master")
    a.add_argument("source"); a.add_argument("out")
    a.add_argument("--height", type=int, default=1080)
    a.add_argument("--start", type=float); a.add_argument("--duration", type=float)
    a.add_argument("--no-deflicker", action="store_true")
    a.add_argument("--no-decimate", action="store_true")
    a.add_argument("--diffs", help="reuse a <out>.diffs.npz from an earlier run")
    b = sub.add_parser("color", help="B&W master -> colour master (resumable)")
    b.add_argument("bw"); b.add_argument("out")
    b.add_argument("--window", type=int, default=5)
    b.add_argument("--chunk", type=int, default=160)
    b.add_argument("--saturation", type=float, default=1.0)
    b.add_argument("--device", default="GPU")
    c = sub.add_parser("still", help="colorize one image")
    c.add_argument("image"); c.add_argument("out")
    args = p.parse_args(argv)
    if args.cmd == "prepare":
        r = prepare(args.source, args.out, args.height, args.start, args.duration,
                    not args.no_deflicker, not args.no_decimate, args.diffs)
        print(json.dumps(asdict(r), indent=2))
    elif args.cmd == "color":
        colorize(args.bw, args.out, Colorizer(device=args.device, saturation=args.saturation),
                 window=args.window, chunk_frames=args.chunk)
    else:
        print(still(args.image, args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
