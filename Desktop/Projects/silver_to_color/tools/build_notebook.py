"""Generate market_street_pipeline.ipynb.

The notebook is generated rather than hand-written so the cell sources stay readable,
reviewable in diffs, and syntax-checkable. Run:

    python build_notebook.py --out ../market_street_pipeline.ipynb
"""

from __future__ import annotations

import argparse
import ast
import json
from pathlib import Path

CELLS: list[tuple[str, str]] = []


def md(text: str) -> None:
    CELLS.append(("markdown", text.strip("\n")))


def code(text: str) -> None:
    CELLS.append(("code", text.strip("\n")))


# ===========================================================================
md(r"""
# Silver to Color, restoration pipeline

Chunked and manifest-resumable, sized for a free Kaggle T4/P100 (16 GB, ~30 GPU-hr/week,
12 h sessions, ~20 GB persistent disk).

**This is a reference implementation.** Diff it against your working notebook and take what
you want; it is not meant to replace yours wholesale.

## What changed, and why

Every item links to the reasoning in `research/pipeline-optimization.md`.

| # | Change | Why |
|---|---|---|
| 0 | **Decimate duplicate frames first** | The 4K master pads 16 fps to 24 fps by doubling every second frame. **One frame in three is an exact duplicate.** Measured: 33.3% of GPU work saved. (§5.6) |
| 1 | **Bilinear chroma upsample** | DDColor's `F.interpolate` has no `mode`, so it defaults to `nearest`. From 512² to 4K that is an 8×6 px block per chroma sample, and sub-pixel drift makes output pixels *snap* between samples, chroma jitter created by the resize, then blamed on the model. (§2.1) |
| 2 | **Motion-compensated chroma median** | The naive temporal median has no notion of motion, so colour lags behind a razor-sharp luminance edge. Flow-warping first: **3.4× lower chroma error at moving edges, measured.** (§2.4) |
| 3 | **Deflicker with `mode=am`, not `median`** | Measured across three flicker profiles: `am`/`gm` cut flicker 78-88%, `median` managed −0.5% to 6.5%. The plausible outlier-rejection argument does not survive measurement. (§2.3) |
| 4 | **Stop upscaling** | A 4K scan of 1906 stock has no sub-grain detail to recover. Deliver 2160p for the **bitrate** (YouTube gives 2160p 35-45 Mbps vs 1080p 8 Mbps), not the resolution. (§0.2) |
| 5 | **Segment from the URL** | The master is **85.4 GB**; `/kaggle/working` is ~20 GB. It cannot be downloaded. (§3) |
| 6 | **fp16 via explicit capability check** | `torch.cuda.is_bf16_supported()` returns **True on a T4** and bf16 is *emulated in software*, silently slower than fp32. (§5.1) |
| 7 | **Per-shot processing** | Deflicker and the temporal median both misbehave across a hard cut. (§2.5) |

## Still open, decide these deliberately

- **`-dn` grain level.** Two sources disagree on which direction the knob runs (§5.5). Render one shot at 1.0 / 0.5 / 0.2 / 0.0 and pick by eye. This is the film's whole look.
- **RIFE interpolation.** A taste call, not a technical one (§6). Silent film was hand-cranked at 12-40 fps; there is no single "true" rate.
- **DDColor `input_size` above 512.** Plumbed through, untested. DeOldify's `render_factor` note suggests it is a *stability* lever, not just quality.
""")

code(r'''
# ---------------------------------------------------------------------------
# CONFIG -- every knob in one place.
# ---------------------------------------------------------------------------
from pathlib import Path

CONFIG = {
    # --- source -------------------------------------------------------------
    # The 4K master is 4096x3072 / 85.4 GB and CANNOT be staged on free tier.
    # The 2K version is already cropped (no sprocket holes), already contrast
    # adjusted, 2.5 GB, and the uploader recommends it. Start here; graduate to
    # the master by URL segmentation once the pipeline is proven.
    "archive_id": "MarketStreet19064KScan20181016",
    "source_file": "MarketStreet_4K_to_2K_cropped_higher_contrast.mov",
    # "source_file": "Market_Street_1906_4K_scan_20181016.mov",   # 85.4 GB master

    # --- timing -------------------------------------------------------------
    # The archive description states the scan was made AT 16 fps, with every
    # second frame doubled for 24 fps playback. So: decimate, then present at 16.
    "native_fps": 16.0,
    "decimate_duplicates": True,

    # --- chunking -----------------------------------------------------------
    "chunk_seconds": 30,
    "snap_to_shots": True,
    "shot_threshold": 0.3,

    # --- deflicker ----------------------------------------------------------
    # mode="am" on measured evidence, not on argument. See research section 2.3.
    "deflicker": True,
    "deflicker_size": 7,
    "deflicker_mode": "am",

    # --- colorization -------------------------------------------------------
    "ddcolor_model": "damo/cv_ddcolor_image-colorization",
    "ddcolor_input_size": 512,      # raising this is an untested stability lever
    "ddcolor_batch": 4,             # ramp up and watch nvidia-smi; est. 8 fp32 / 16 fp16
    "temporal_window": 5,           # odd. Longer = steadier, more lag at flow failures.
    "fb_threshold_px": 1.0,         # forward-backward occlusion tolerance

    # --- super-resolution ---------------------------------------------------
    # Cleanup at delivery resolution, NOT upscaling. outscale=1 because the source
    # is already beyond what the emulsion recorded.
    "run_esrgan": True,
    "esrgan_outscale": 1,
    "esrgan_denoise": 0.5,          # A/B this at 1.0/0.5/0.2/0.0 -- see section 5.5
    "esrgan_tile": 512,             # 0 to disable; try 0 first and watch nvidia-smi
    "esrgan_tile_pad": 16,

    # --- delivery -----------------------------------------------------------
    "deliver_height": 2160,         # 2160p for the bitrate ladder, not for detail
    "crf": 16,
    "preset": "slow",

    # --- paths --------------------------------------------------------------
    # /kaggle/working persists via Save Version but is only ~20 GB.
    # /kaggle/temp is ~60 GB and does NOT persist -- fine for intermediates.
    "work_dir": "/kaggle/working/stc",
    "temp_dir": "/kaggle/temp/stc",
    "tools_dir": ".",               # where ia_meta.py / pipeline.py / temporal_chroma.py live
}

WORK = Path(CONFIG["work_dir"])
TEMP = Path(CONFIG["temp_dir"]) if Path(CONFIG["temp_dir"]).parent.exists() else WORK / "temp"
for d in (WORK, TEMP):
    d.mkdir(parents=True, exist_ok=True)

MANIFEST_PATH = WORK / "manifest.json"   # lives on the PERSISTENT disk, deliberately
print(f"work: {WORK}\ntemp: {TEMP}")
''')

md(r"""
## Environment

The capability check matters more than it looks. `torch.cuda.is_bf16_supported()` returns
**True on a T4**, it defaults to `including_emulation=True` and merely checks that a bf16
tensor can be allocated. Turing has no bf16 hardware; PyTorch emulates it in software, so
bf16 autocast runs **silently slower than fp32**. Read the compute capability directly.
""")

code(r'''
import subprocess, sys, shutil, json, os

def sh(cmd, **kw):
    return subprocess.run(cmd, shell=isinstance(cmd, str), capture_output=True, text=True, **kw)

print("python :", sys.version.split()[0])
print("ffmpeg :", "yes" if shutil.which("ffmpeg") else "MISSING")

import torch
print("torch  :", torch.__version__, "| cuda:", torch.cuda.is_available())

AMP_DTYPE = None
if torch.cuda.is_available():
    name = torch.cuda.get_device_name(0)
    major, minor = torch.cuda.get_device_capability()
    total_gib = torch.cuda.get_device_properties(0).total_memory / 1024**3
    print(f"gpu    : {name} (sm_{major}{minor}), {total_gib:.2f} GiB")

    # Do NOT use torch.cuda.is_bf16_supported() -- it reports True on Turing.
    if major >= 8:
        AMP_DTYPE = torch.bfloat16
        note = "Ampere+: bf16 has hardware support"
    else:
        AMP_DTYPE = torch.float16
        note = ("Turing/Pascal: fp16 only. bf16 would be EMULATED in software "
                "and slower than fp32")
    print(f"amp    : {AMP_DTYPE} ({note})")

    if major == 6:
        print("         NOTE: P100 has no tensor cores -- expect ~1.2-1.5x from fp16, not 3x.")
    elif major == 7:
        print("         T4: 65 TFLOPS fp16 vs 8.1 fp32. Usable VRAM ~15,360 MiB (ECC reserves 6.25%).")

    torch.backends.cudnn.benchmark = True   # safe: DDColor's input shape is always 512x512
else:
    print("gpu    : NONE -- CPU-only stages will still run")

import shutil as _sh
for label, path in (("work", CONFIG["work_dir"]), ("temp", str(TEMP))):
    try:
        u = _sh.disk_usage(path)
        print(f"disk {label}: {u.free/1e9:.1f} GB free of {u.total/1e9:.1f} GB")
    except OSError:
        pass
''')

code(r'''
# Dependencies. Skip anything already present so a resumed session starts fast.
def ensure(pkg, pip_name=None):
    try:
        __import__(pkg)
        return True
    except ImportError:
        r = sh([sys.executable, "-m", "pip", "install", "-q", pip_name or pkg])
        if r.returncode != 0:
            print(f"  ! {pkg} failed:\n{r.stderr[-600:]}")
            return False
        return True

for pkg, pip_name in [("cv2", "opencv-python-headless"),
                      ("numpy", "numpy"),
                      ("modelscope", "modelscope"),
                      ("realesrgan", "realesrgan"),
                      ("basicsr", "basicsr")]:
    print(f"{pkg:12s}", "ok" if ensure(pkg, pip_name) else "FAILED")

# The project's own tested modules. Attach tools/ as a Kaggle dataset, or upload them.
sys.path.insert(0, CONFIG["tools_dir"])
import ia_meta, pipeline as pl, temporal_chroma as tc
print("\nlocal modules loaded: ia_meta, pipeline, temporal_chroma")
''')

md(r"""
## Stage 0, source, verification, and a flicker measurement

Three things happen here and none should be skipped.

**Verify the frame rate rather than assuming it.** The 4K master reports 404.5 s while its own
derivative reports 716.8 s, because the master pads 16 fps to 24 fps. A wrong rate is baked
into every frame.

**Measure the flicker before spending a pass on it.** `deflicker=bypass=1` writes per-frame
metadata without touching a single pixel (verified bit-identical). If the source is already
steady, skip the pass.

**Record provenance now**, while it is cheap. "Public domain" is not one of YouTube's
enumerated Content ID dispute reasons, you file under "misidentified", which means naming the
item, its rights metadata and the film's year.
""")

code(r'''
item = ia_meta.fetch_item(CONFIG["archive_id"])
meta = item["metadata"]
print(meta.get("title"))
print("collection:", ia_meta.as_list(meta.get("collection")))

# archive.org carries NO rights metadata for this item. Absence is not evidence either
# way -- PD status here rests on the 1906 publication date alone.
for f in ("licenseurl", "rights", "possible-copyright-status"):
    print(f"  {f:28s} {meta.get(f)}")

videos = ia_meta.video_files(item)
print("\navailable video files:")
for v in videos[:8]:
    res = f"{v['width']}x{v['height']}" if v["width"] else "res not stated"
    print(f"  {v['name'][:52]:52s} {res:>13}  {ia_meta.human_size(v['size']):>9}  {v['source']}")

SOURCE_URL = ia_meta.download_url(CONFIG["archive_id"], CONFIG["source_file"])
print("\nusing:", CONFIG["source_file"])
''')

code(r'''
# Probe WITHOUT downloading -- ffprobe reads the URL with range requests.
info = pl.probe_video(SOURCE_URL)
print(f"{info.width}x{info.height}  {info.codec}")
print(f"avg_frame_rate : {info.avg_frame_rate}")
print(f"duration       : {info.duration:.2f}s" if info.duration else "duration: unknown")
print(f"nb_frames      : {info.nb_frames}")

if info.duration is None:
    raise SystemExit("no duration -- cannot plan chunks. Try a different source file.")

# Cross-check the container's own numbers against each other.
implied = info.implied_frame_count
if implied and info.nb_frames:
    drift = abs(implied - info.nb_frames) / info.nb_frames
    print(f"frame-count cross-check: nb_frames={info.nb_frames}, "
          f"duration*fps={implied:.0f} ({drift:.1%} drift)")
    if drift > 0.02:
        print("  ! These disagree. Padded/duplicated frames are the usual cause.")
''')

code(r'''
# Flicker measurement -- modifies nothing. On a long source this reads the whole file,
# so sample a couple of minutes rather than the lot.
SAMPLE = TEMP / "sample.mkv"
if not SAMPLE.exists():
    pl.extract_chunk(SOURCE_URL, SAMPLE, start=0, duration=min(120, info.duration))

flicker = pl.flicker_summary(pl.measure_flicker(str(SAMPLE)))
print(f"frames measured    : {flicker['frames']}")
print(f"mean |rel change|  : {flicker['mean_abs']:.4f}")
print(f"p95  |rel change|  : {flicker['p95_abs']:.4f}")
print(f"worth deflickering : {flicker['worth_deflickering']}")

if flicker["worth_deflickering"] is False:
    print("\n-> Source is already steady. Consider CONFIG['deflicker'] = False and")
    print("   spend the time on defect repair (RRTN) instead.")
''')

code(r'''
# Provenance record, written before anything is rendered.
import datetime as _dt

best = next((v for v in videos if v["name"] == CONFIG["source_file"]), None)
provenance = {
    "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
    "archive_identifier": CONFIG["archive_id"],
    "archive_details_url": ia_meta.details_url(CONFIG["archive_id"]),
    "archive_title": meta.get("title"),
    "archive_uploader": meta.get("uploader"),
    "uploader_reuse_statement": (
        "Anyone may reproduce or reuse this scan. Please attribute it to its source: "
        "Prelinger Archives."
    ),
    "film_year": 1906,
    "pd_basis": "Published 1906, before the 1931 cutoff (95-year term, calendar 2026).",
    "archive_rights_metadata": {k: meta.get(k) for k in
                                ("licenseurl", "rights", "possible-copyright-status")},
    "rights_metadata_note": "archive.org carries NO rights metadata for this item.",
    "source_file": best,
    "source_file_url": SOURCE_URL,
    "probe": {"width": info.width, "height": info.height,
              "avg_frame_rate": info.avg_frame_rate, "duration": info.duration},
    "flicker": flicker,
    "dispute_note": ("'Public domain' is NOT an enumerated YouTube dispute reason. File under "
                     "'misidentified or an error was made' and cite this item URL plus the "
                     "film year; PD footage is ineligible as a Content ID reference "
                     "(support.google.com/youtube/answer/2605065)."),
}
(WORK / "provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
print("wrote", WORK / "provenance.json")
''')

md(r"""
## Stage 1, plan the chunks

Chunk boundaries are snapped to shot boundaries where one falls nearby, because both the
deflicker and the temporal median assume their window lies inside a single continuous shot.
Across a hard cut, deflicker reads a legitimate brightness change as flicker and pumps against
it, and the median blends colour from two unrelated scenes.

For a single-take 1906 street film this may find no cuts at all, which is a useful answer
rather than a failure.

The manifest lives on the **persistent** disk and is written atomically, so a session killed
mid-write cannot lose the record of what is already done.
""")

code(r'''
shots = []
if CONFIG["snap_to_shots"]:
    print("detecting shot boundaries (reads the file once)...")
    shots = pl.detect_shots(str(SAMPLE), threshold=CONFIG["shot_threshold"])
    print(f"  {len(shots)} boundaries in the sampled region: {[round(t,1) for t in shots[:12]]}")

manifest = pl.Manifest(MANIFEST_PATH)
if not manifest.all_chunks():
    chunks = pl.plan_chunks(info.duration, CONFIG["chunk_seconds"], shot_boundaries=shots)
    manifest.set_chunks(chunks)
    manifest.data["source"] = SOURCE_URL
    manifest.data["config"] = {k: v for k, v in CONFIG.items() if not k.endswith("_dir")}
    manifest.save()
    print(f"planned {len(chunks)} chunks")
else:
    print(f"resuming existing manifest: {MANIFEST_PATH}")

p = manifest.progress()
print(f"progress: {p['done']}/{p['total']} done ({p['percent']:.1f}%), "
      f"{p['failed']} failed, {p['pending']} pending")
''')

md(r"""
## Stage 2, CPU preparation, per chunk

Order matters. **Decimate first**: the master repeats one frame in three, and there is no
reason to deflicker, colorize or upscale a duplicate. Measured saving on a clip built to the
archive's own 16→24 recipe: **48 frames in, 32 out, 33.3% of GPU work removed**, recovering
exactly the unique frame count.

Then deflicker, because DDColor's *only* input is luminance, flicker is a direct perturbation
of the model's input, frame by frame, on a network with no temporal state.

One subtlety worth not rediscovering: the re-timing after decimation happens with `setpts`
**inside** the filter graph plus `-fps_mode passthrough`. Using an output `-r` instead makes
ffmpeg re-time a stream whose timestamps still describe the old rate, and it drops real frames.
""")

code(r'''
import time

def prepare_chunk(chunk):
    """Extract -> decimate -> deflicker. All CPU; safe to run without a GPU."""
    tag = f"{chunk.index:04d}"

    raw = manifest.stage_done(chunk.index, "extract")
    if not raw:
        raw = str(pl.extract_chunk(SOURCE_URL, TEMP / f"raw{tag}.mkv",
                                   chunk.start, chunk.duration))
        manifest.mark_stage(chunk.index, "extract", raw)

    current = raw
    if CONFIG["decimate_duplicates"]:
        dec = manifest.stage_done(chunk.index, "decimate")
        if not dec:
            dec = str(pl.decimate_chunk(current, TEMP / f"dec{tag}.mkv",
                                        frame_rate=CONFIG["native_fps"]))
            manifest.mark_stage(chunk.index, "decimate", dec)
        current = dec

    if CONFIG["deflicker"]:
        df = manifest.stage_done(chunk.index, "deflicker")
        if not df:
            df = str(pl.deflicker_chunk(current, TEMP / f"df{tag}.mkv",
                                        size=CONFIG["deflicker_size"],
                                        mode=CONFIG["deflicker_mode"]))
            manifest.mark_stage(chunk.index, "deflicker", df)
        current = df

    return current


pending = manifest.pending()
print(f"{len(pending)} chunks to prepare\n")
started = time.time()
for n, chunk in enumerate(pending, 1):
    try:
        out = prepare_chunk(chunk)
        if n == 1 and CONFIG["decimate_duplicates"]:
            before = pl.count_frames(manifest.stage_done(chunk.index, "extract"))
            after = pl.count_frames(manifest.stage_done(chunk.index, "decimate"))
            print(f"  decimation on chunk 0: {before} -> {after} frames "
                  f"({100*(before-after)/before:.1f}% of GPU work saved)\n")
        print(f"  [{n}/{len(pending)}] chunk {chunk.index} prepared")
    except pl.PipelineError as exc:
        manifest.mark_failed(chunk.index, str(exc))
        print(f"  [{n}/{len(pending)}] chunk {chunk.index} FAILED: {str(exc)[:200]}")
print(f"\nelapsed {time.time()-started:.0f}s")
''')

md(r"""
## Stage 3, colorization

Two patches to DDColor, both small and both consequential.

**Return the raw 512×512 `ab` tensor** instead of an upsampled frame. All temporal work then
happens on `(2, 512, 512)` arrays rather than `(3, 4096, 3072)`, 32× less arithmetic to filter
what is, after upsampling, a replica of the same 512² samples. Optical flow only needs computing
at 512×512 for the same reason, which is what makes motion compensation nearly free.

**Upsample bilinearly at the end.** The shipped pipeline calls `F.interpolate` with no `mode`,
so PyTorch defaults to `nearest`.

Note what is structurally protected: `orig_l` is taken at **full resolution** and never enters
the network. No amount of chroma smoothing can soften an edge, so the usual "temporal filtering
costs sharpness" objection does not apply here.

DDColor has no fp16 support of its own. `torch.autocast` is the right wrapper rather than
`.half()`, because it keeps the LayerNorms and softmaxes in fp32. **Untested, validate on ~20
frames against an fp32 reference before committing a whole film.**
""")

code(r'''
import numpy as np, cv2, torch

class DDColorChroma:
    """DDColor wrapper that exposes raw 512x512 ab and upsamples bilinearly.

    Mirrors ColorizationPipeline.process() with three changes:
      1. returns the raw (2, S, S) ab tensor so temporal filtering happens at 512x512
      2. batches frames instead of the shipped .unsqueeze(0) batch-of-one
      3. bilinear (not nearest) when the chroma finally goes to full resolution
    """

    def __init__(self, model_path, input_size=512, device=None, amp_dtype=None):
        from modelscope.outputs import OutputKeys          # noqa: F401
        from modelscope.pipelines import pipeline as ms_pipeline
        from modelscope.utils.constant import Tasks

        self.input_size = int(input_size)
        self.device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        self.amp_dtype = amp_dtype
        self._ms = ms_pipeline(Tasks.image_colorization, model=model_path)
        self.model = getattr(self._ms, "model", None)
        if self.model is not None:
            self.model.eval()

    @staticmethod
    def to_lab_l(bgr_float):
        """Full-resolution L. This never enters the network."""
        return cv2.cvtColor(bgr_float, cv2.COLOR_BGR2Lab)[:, :, :1]

    def _grey_rgb_batch(self, frames_bgr_float):
        """Model input: L only, with a and b explicitly zeroed."""
        batch = []
        for img in frames_bgr_float:
            small = cv2.resize(img, (self.input_size, self.input_size),
                               interpolation=cv2.INTER_AREA)
            l = cv2.cvtColor(small, cv2.COLOR_BGR2Lab)[:, :, :1]
            grey_lab = np.concatenate([l, np.zeros_like(l), np.zeros_like(l)], axis=-1)
            batch.append(cv2.cvtColor(grey_lab, cv2.COLOR_LAB2RGB).transpose(2, 0, 1))
        return torch.from_numpy(np.stack(batch)).float().to(self.device)

    @torch.no_grad()
    def predict_ab(self, frames_bgr_uint8):
        """-> (N, S, S, 2) float32 chroma at NATIVE 512x512. No upsampling yet."""
        floats = [f.astype(np.float32) / 255.0 for f in frames_bgr_uint8]
        tensor = self._grey_rgb_batch(floats)
        if self.amp_dtype is not None and self.device == "cuda":
            with torch.autocast("cuda", dtype=self.amp_dtype):
                out = self.model(tensor)
        else:
            out = self.model(tensor)
        return out.float().cpu().numpy().transpose(0, 2, 3, 1)

    @staticmethod
    def recombine(orig_l, ab_native):
        """Bilinear chroma upsample, then rejoin the untouched full-resolution L."""
        h, w = orig_l.shape[:2]
        ab = tc.upsample_chroma(ab_native, h, w)          # INTER_LINEAR, not nearest
        lab = np.concatenate([orig_l, ab], axis=-1)
        bgr = cv2.cvtColor(lab, cv2.COLOR_LAB2BGR)
        return np.clip(bgr * 255.0, 0, 255).astype(np.uint8)


print("DDColorChroma defined. Instantiate when you are ready to burn GPU time:")
print('  colorizer = DDColorChroma(CONFIG["ddcolor_model"],')
print('                            input_size=CONFIG["ddcolor_input_size"],')
print("                            amp_dtype=AMP_DTYPE)")
''')

code(r'''
def colorize_chunk(in_path, out_path, colorizer, window=None, batch=None):
    """Colorize one chunk with motion-compensated temporal chroma smoothing.

    Frames stream in and out through ffmpeg pipes -- no PNG folders. That matters on
    Kaggle beyond speed: notebook output is capped at 500 files, and a chunk of 4K PNGs
    blows through it immediately.
    """
    window = window or CONFIG["temporal_window"]
    batch = batch or CONFIG["ddcolor_batch"]

    probe = pl.probe_video(str(in_path))
    w, h = probe.width, probe.height
    fps = probe.avg_frame_rate or CONFIG["native_fps"]

    reader = subprocess.Popen(
        ["ffmpeg", "-v", "error", "-i", str(in_path),
         "-f", "rawvideo", "-pix_fmt", "bgr24", "-"],
        stdout=subprocess.PIPE)
    writer = subprocess.Popen(
        ["ffmpeg", "-v", "error", "-y", "-f", "rawvideo", "-pix_fmt", "bgr24",
         "-s", f"{w}x{h}", "-r", f"{fps:.6f}", "-i", "-",
         "-an", "-c:v", "ffv1", "-level", "3", str(out_path)],
        stdin=subprocess.PIPE)

    frame_bytes = w * h * 3
    flow_engine = tc.make_dis_flow("medium")
    half = window // 2

    buf_bgr, buf_ab, buf_gray, buf_l = [], [], [], []
    emitted = 0

    def flush(upto):
        """Emit frames whose full temporal window is available."""
        nonlocal emitted
        while len(buf_ab) > upto:
            centre = min(half, len(buf_ab) - 1)
            if len(buf_ab) < window and upto > 0:
                break
            lo = max(0, centre - half)
            hi = min(len(buf_ab), centre + half + 1)
            ab = tc.flow_warped_median(
                np.stack(buf_ab[lo:hi]), np.stack(buf_gray[lo:hi]),
                centre - lo, flow_engine=flow_engine,
                fb_threshold_px=CONFIG["fb_threshold_px"])
            writer.stdin.write(DDColorChroma.recombine(buf_l[centre], ab).tobytes())
            emitted += 1
            for b in (buf_bgr, buf_ab, buf_gray, buf_l):
                b.pop(0)

    try:
        while True:
            frames = []
            for _ in range(batch):
                raw = reader.stdout.read(frame_bytes)
                if not raw or len(raw) < frame_bytes:
                    break
                frames.append(np.frombuffer(raw, np.uint8).reshape(h, w, 3))
            if not frames:
                break

            ab_batch = colorizer.predict_ab(frames)
            for frame, ab in zip(frames, ab_batch):
                f32 = frame.astype(np.float32) / 255.0
                buf_bgr.append(frame)
                buf_ab.append(ab.astype(np.float32))
                buf_l.append(DDColorChroma.to_lab_l(f32))
                # Flow is estimated at the chroma's native size, on luminance.
                buf_gray.append(cv2.resize(
                    cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY),
                    (ab.shape[1], ab.shape[0]), interpolation=cv2.INTER_AREA))
            flush(window - 1)
        flush(0)          # drain the tail
    finally:
        if writer.stdin:
            writer.stdin.close()
        writer.wait()
        reader.stdout.close()
        reader.wait()

    return emitted


print("colorize_chunk defined.")
print("Run it over manifest.pending() the same way Stage 2 does, marking a 'colorize' stage.")
''')

md(r"""
## Stage 4, Real-ESRGAN as a cleanup pass

Not as an upscaler. `outscale=1` because a 4K scan of 1906 stock has no sub-grain detail to
recover, the emulsion never recorded any, so what a 4× network reads as texture to sharpen is
silver-halide grain and dust. The specialists agree by example: the leading old-film restoration
paper sets `scale: 1`.

Two gotchas:

- **`RealESRGANer.__init__` defaults `half=False`.** The CLI defaults to fp16, but a notebook
  calling the library gets fp32 unless it passes `half=True`, silently halving throughput.
- **`--outscale 2` saves nothing.** The network always runs 4×, then resizes on the CPU.

`denoise_strength` is the film-look control. Two sources disagree about which direction it runs,
so label your A/B renders by the flag value rather than by "more grain".
""")

code(r'''
def build_upsampler():
    from basicsr.archs.srvgg_arch import SRVGGNetCompact
    from realesrgan import RealESRGANer
    from realesrgan.archs.srvgg_arch import SRVGGNetCompact as _SRVGG   # noqa: F401

    model = SRVGGNetCompact(num_in_ch=3, num_out_ch=3, num_feat=64,
                            num_conv=32, upscale=4, act_type="prelu")

    base = ("https://github.com/xinntao/Real-ESRGAN/releases/download/v0.2.5.0/"
            "realesr-general-x4v3.pth")
    wdn = ("https://github.com/xinntao/Real-ESRGAN/releases/download/v0.2.5.0/"
           "realesr-general-wdn-x4v3.pth")

    dn = float(CONFIG["esrgan_denoise"])
    if dn != 1.0:
        # Deep network interpolation: a weight-space blend of the two checkpoints,
        # evaluated once at load. Free at runtime; both files must be present.
        model_path, dni = [base, wdn], [dn, 1 - dn]
    else:
        model_path, dni = base, None

    return RealESRGANer(
        scale=4, model_path=model_path, dni_weight=dni, model=model,
        tile=int(CONFIG["esrgan_tile"]), tile_pad=int(CONFIG["esrgan_tile_pad"]),
        pre_pad=0,
        half=(AMP_DTYPE == torch.float16),   # the library default is False -- pass it
    )


def esrgan_chunk(in_path, out_path, upsampler):
    probe = pl.probe_video(str(in_path))
    w, h = probe.width, probe.height
    fps = probe.avg_frame_rate or CONFIG["native_fps"]
    out_scale = float(CONFIG["esrgan_outscale"])
    ow, oh = int(round(w * out_scale)), int(round(h * out_scale))

    reader = subprocess.Popen(
        ["ffmpeg", "-v", "error", "-i", str(in_path),
         "-f", "rawvideo", "-pix_fmt", "bgr24", "-"], stdout=subprocess.PIPE)
    writer = subprocess.Popen(
        ["ffmpeg", "-v", "error", "-y", "-f", "rawvideo", "-pix_fmt", "bgr24",
         "-s", f"{ow}x{oh}", "-r", f"{fps:.6f}", "-i", "-",
         "-an", "-c:v", "ffv1", "-level", "3", str(out_path)], stdin=subprocess.PIPE)

    frame_bytes, n = w * h * 3, 0
    try:
        while True:
            raw = reader.stdout.read(frame_bytes)
            if not raw or len(raw) < frame_bytes:
                break
            frame = np.frombuffer(raw, np.uint8).reshape(h, w, 3)
            out, _ = upsampler.enhance(frame, outscale=out_scale)
            writer.stdin.write(np.ascontiguousarray(out).tobytes())
            n += 1
    finally:
        if writer.stdin:
            writer.stdin.close()
        writer.wait()
        reader.stdout.close()
        reader.wait()
    return n


print("build_upsampler / esrgan_chunk defined.")
print(f"denoise_strength = {CONFIG['esrgan_denoise']}  <- A/B at 1.0 / 0.5 / 0.2 / 0.0")
''')

md(r"""
## Stage 5, assemble and compare

`assemble()` refuses to run if any chunk is missing rather than quietly producing a short film, 
the exact failure the manifest exists to prevent.

The before/after render is not decoration. It makes the transformation visible to a YouTube
reviewer within seconds, which is the reused-content policy's "viewers can tell that there's a
meaningful difference" test, and it is the short-form hook.
""")

code(r'''
FINAL_STAGE = "esrgan" if CONFIG["run_esrgan"] else "colorize"

outputs, missing = [], []
for chunk in manifest.all_chunks():
    path = manifest.stage_done(chunk.index, FINAL_STAGE)
    (outputs if path else missing).append(path or chunk.index)

print(f"{len(outputs)}/{len(manifest.all_chunks())} chunks have a '{FINAL_STAGE}' output")
if missing:
    print(f"still to do: {missing[:20]}{' ...' if len(missing) > 20 else ''}")
else:
    final = WORK / "market_street_1906_restored.mp4"
    pl.assemble(outputs, final, frame_rate=CONFIG["native_fps"],
                crf=CONFIG["crf"], preset=CONFIG["preset"])
    out_info = pl.probe_video(str(final))
    print(f"\nwrote {final}")
    print(f"  {out_info.width}x{out_info.height}  {out_info.duration:.1f}s  "
          f"{final.stat().st_size/1e6:.0f} MB")

    # Duration must match the source, allowing for decimation changing the frame count
    # but NOT the running time.
    if out_info.duration and info.duration:
        drift = abs(out_info.duration - info.duration) / info.duration
        print(f"  duration drift vs source: {drift:.2%}")
        if drift > 0.05:
            print("  ! Over 5% drift. Check the decimate re-timing and the assemble frame rate.")
''')

code(r'''
# Before/after, from the same time range in both files.
COMPARE_START, COMPARE_SECONDS = 0, 20

before_clip = TEMP / "cmp_before.mkv"
after_clip = TEMP / "cmp_after.mkv"
final = WORK / "market_street_1906_restored.mp4"

if final.exists():
    pl.extract_chunk(SOURCE_URL, before_clip, COMPARE_START, COMPARE_SECONDS)
    pl.extract_chunk(str(final), after_clip, COMPARE_START, COMPARE_SECONDS)

    for mode in ("sidebyside", "wipe"):
        out = WORK / f"compare_{mode}.mp4"
        pl.make_comparison(before_clip, after_clip, out, mode=mode)
        print("wrote", out)
else:
    print("assemble the film first")
''')

md(r"""
## Stage 6, quality control

Two measurements that catch the failures worth catching, both cheap.

**Residual flicker**, did the deflicker actually work, and did colorization introduce new
instability?

**Chroma stability**, mean frame-to-frame change in the a/b channels. This is the number the
temporal median exists to reduce, so it is the one that tells you whether the window length is
right. Judge it on *moving* footage; a static clip makes any setting look fine.
""")

code(r'''
def chroma_stability(path, max_frames=300, sample_size=256):
    """Mean frame-to-frame |delta| in Lab a/b. Lower is steadier.

    Measured on a downscaled copy: chroma is a smooth signal and the statistic is about
    temporal behaviour, not spatial detail, so full resolution buys nothing here.
    """
    proc = subprocess.Popen(
        ["ffmpeg", "-v", "error", "-i", str(path), "-vf",
         f"scale={sample_size}:{sample_size}", "-f", "rawvideo", "-pix_fmt", "bgr24", "-"],
        stdout=subprocess.PIPE)
    n_bytes = sample_size * sample_size * 3
    prev, deltas = None, []
    try:
        for _ in range(max_frames):
            raw = proc.stdout.read(n_bytes)
            if not raw or len(raw) < n_bytes:
                break
            bgr = np.frombuffer(raw, np.uint8).reshape(sample_size, sample_size, 3)
            ab = cv2.cvtColor(bgr.astype(np.float32) / 255.0, cv2.COLOR_BGR2Lab)[::, 1:]
            if prev is not None:
                deltas.append(float(np.abs(ab - prev).mean()))
            prev = ab
    finally:
        proc.stdout.close()
        proc.wait()
    return float(np.mean(deltas)) if deltas else float("nan")


final = WORK / "market_street_1906_restored.mp4"
if final.exists():
    src_sample = SAMPLE
    print(f"{'':22s} {'flicker':>10s} {'chroma drift':>14s}")
    for label, path in (("source", src_sample), ("restored", final)):
        fl = pl.flicker_summary(pl.measure_flicker(str(path)))["mean_abs"]
        ch = chroma_stability(path)
        print(f"{label:22s} {fl:10.4f} {ch:14.4f}")
    print("\nFlicker should fall. Chroma drift on the restored file is the number to tune")
    print("CONFIG['temporal_window'] against -- raise it if the drift is high, but check")
    print("moving shots for colour lagging behind edges before you do.")
else:
    print("assemble the film first")
''')

md(r"""
## What this notebook deliberately does not do

- **Defect repair.** Nothing here removes a scratch, blotch or splice mark, which is what a
  viewer of a 1906 print actually notices. RRTN (`mountln/RRTN-old-film-restoration`, MIT,
  torch 2.1.2, weights on GitHub Releases) is the right tool and the highest-value addition, 
  but budget a whole session for the `mim install mmcv` and JIT `deform_attn` build, then cache
  the compiled artifacts as a Kaggle dataset. Run it **before** colorization, at ≤1080p, and
  transfer the repaired luma back. Note `--max_recursion 4`: budget up to 4× the single-pass
  cost.
- **Exemplar colorization.** CMNET2 is the cure that this notebook's temporal median only
  treats the symptom of. Measured 2.63 fps at 960×730 on 16 GB, 5.80 with `--max_side 512`.
- **Stabilisation.** `vidstab` is worth running if the print has gate weave, not for its own
  sake, but because weave makes optical flow harder, and Stage 3 depends on flow.
- **Frame interpolation.** A taste decision (§6). If you do it, do it **last**, RIFE will
  happily interpolate a scratch into a smooth moving object.
""")
# ===========================================================================


def build(out_path: Path) -> None:
    cells = []
    for kind, source in CELLS:
        lines = source.splitlines(keepends=True)
        if kind == "code":
            cells.append({
                "cell_type": "code", "execution_count": None,
                "metadata": {}, "outputs": [], "source": lines,
            })
        else:
            cells.append({"cell_type": "markdown", "metadata": {}, "source": lines})

    notebook = {
        "cells": cells,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python",
                           "name": "python3"},
            "language_info": {"name": "python", "version": "3.10"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(notebook, indent=1), encoding="utf-8")


def validate(out_path: Path) -> int:
    """Parse the notebook back and syntax-check every code cell."""
    notebook = json.loads(out_path.read_text(encoding="utf-8"))
    code_cells = [c for c in notebook["cells"] if c["cell_type"] == "code"]
    md_cells = [c for c in notebook["cells"] if c["cell_type"] == "markdown"]

    failures = 0
    for index, cell in enumerate(code_cells):
        source = "".join(cell["source"])
        try:
            ast.parse(source)
        except SyntaxError as exc:
            failures += 1
            print(f"  SYNTAX ERROR in code cell {index}: line {exc.lineno}: {exc.msg}")
            for line in (source.splitlines()[
                    max(0, (exc.lineno or 1) - 3):(exc.lineno or 1) + 1]):
                print(f"      {line}")

    print(f"cells: {len(notebook['cells'])} "
          f"({len(code_cells)} code, {len(md_cells)} markdown)")
    print(f"syntax: {len(code_cells) - failures}/{len(code_cells)} code cells parse")
    return failures


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="../market_street_pipeline.ipynb")
    args = parser.parse_args()

    out_path = Path(args.out)
    build(out_path)
    failures = validate(out_path)
    print(f"wrote {out_path} ({out_path.stat().st_size/1024:.0f} KB)")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
