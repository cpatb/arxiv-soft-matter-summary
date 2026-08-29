"""Chunked, manifest-resumable film restoration pipeline: the non-GPU half.

Everything here runs without a GPU and is covered by `test_pipeline.py`, so the notebook
can orchestrate tested code instead of embedding untested cells. The GPU stages
(DDColor, Real-ESRGAN, RRTN) live in the notebook because their dependencies only exist
on Kaggle/Colab.

Design constraints this exists to satisfy:

  * The 4K master of the first film is 85.4 GB. Kaggle's /kaggle/working is ~20 GB and
    /kaggle/temp is ~60 GB and does not persist. The source cannot be downloaded, so
    segmentation reads byte ranges straight from the URL.
  * Kaggle sessions cap at 12 h and recycle after 60 min idle, and the weekly GPU quota
    is ~30 h. Work must resume at chunk granularity across sessions.
  * Notebook output is capped at 500 files, so each chunk is ONE container file, never a
    frame folder.

Standard library plus a working `ffmpeg`/`ffprobe` on PATH.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

# Chunk length in seconds. 30 s at ~17 fps is ~510 frames, which is a few minutes of GPU
# work -- short enough that losing one to a session timeout is cheap, long enough that
# per-chunk overhead stays negligible.
DEFAULT_CHUNK_SECONDS = 30

# Lossless intermediate. ffv1 compresses better; utvideo encodes ~10x faster. On Kaggle
# the scarce quota-limited resource is GPU time, not disk, so speed usually wins -- but
# ffv1 is the default here because the 20 GB working dir is the tighter constraint.
LOSSLESS_CODEC = "ffv1"

# A trailing chunk shorter than this fraction of the nominal length is merged into its
# predecessor instead of standing alone.
MIN_CHUNK_FRACTION = 0.25

STATUS_PENDING = "pending"
STATUS_DONE = "done"
STATUS_FAILED = "failed"


class PipelineError(RuntimeError):
    """A pipeline stage could not complete."""


# ---------------------------------------------------------------------------
# ffmpeg / ffprobe helpers
# ---------------------------------------------------------------------------

def _run(cmd: list[str], timeout: float = 3600.0,
         cwd: str | Path | None = None) -> subprocess.CompletedProcess:
    """Run a subprocess, raising PipelineError with stderr on failure."""
    try:
        proc = subprocess.run(
            cmd, capture_output=True, timeout=timeout,
            cwd=str(cwd) if cwd is not None else None,
        )
    except FileNotFoundError as exc:
        raise PipelineError(f"{cmd[0]} not found on PATH") from exc
    except subprocess.TimeoutExpired as exc:
        raise PipelineError(f"{cmd[0]} timed out after {timeout}s") from exc
    if proc.returncode != 0:
        tail = proc.stderr.decode("utf-8", "replace")[-2000:]
        raise PipelineError(f"{cmd[0]} failed (exit {proc.returncode}):\n{tail}")
    return proc


def _absolute_source(source: str) -> str:
    """Absolutise a local path so it survives a cwd change; pass URLs through."""
    if "://" in source:
        return source
    path = Path(source)
    return str(path.resolve()) if path.exists() else source


def have_tool(name: str) -> bool:
    return shutil.which(name) is not None


def _parse_fraction(text: str) -> float | None:
    """ffprobe reports frame rates as 'num/den'."""
    if not text or text in ("0/0", "N/A"):
        return None
    if "/" in text:
        num, _, den = text.partition("/")
        try:
            numerator, denominator = float(num), float(den)
        except ValueError:
            return None
        return numerator / denominator if denominator else None
    try:
        return float(text)
    except ValueError:
        return None


@dataclass
class VideoInfo:
    """What ffprobe knows about a video stream."""
    width: int | None = None
    height: int | None = None
    r_frame_rate: float | None = None
    avg_frame_rate: float | None = None
    nb_frames: int | None = None
    duration: float | None = None
    codec: str | None = None

    @property
    def implied_frame_count(self) -> float | None:
        """Frames implied by duration x rate, for cross-checking nb_frames."""
        if self.duration and self.avg_frame_rate:
            return self.duration * self.avg_frame_rate
        return None


def probe_video(source: str, timeout: float = 120.0) -> VideoInfo:
    """Read stream properties. Works on a local path or an http(s) URL.

    Probing the URL avoids downloading the file, which matters when the file is 85 GB.
    """
    if not source:
        raise ValueError("source is empty")
    # `format=duration` is requested as well as the stream entry because several
    # containers -- Matroska notably -- carry no stream-level duration at all, and a
    # silent None here propagates into chunk planning as "nothing to do".
    proc = _run(
        [
            "ffprobe", "-v", "error", "-select_streams", "v:0",
            "-show_entries",
            "stream=width,height,r_frame_rate,avg_frame_rate,nb_frames,duration,codec_name"
            ":format=duration",
            "-of", "json", source,
        ],
        timeout=timeout,
    )
    try:
        payload = json.loads(proc.stdout.decode("utf-8"))
    except json.JSONDecodeError as exc:
        raise PipelineError(f"ffprobe returned unparseable JSON for {source}") from exc

    streams = payload.get("streams") or []
    if not streams:
        raise PipelineError(f"no video stream found in {source}")
    stream = streams[0]

    def as_int(key):
        try:
            return int(stream[key])
        except (KeyError, TypeError, ValueError):
            return None

    def as_float(key):
        try:
            return float(stream[key])
        except (KeyError, TypeError, ValueError):
            return None

    rate = (_parse_fraction(stream.get("avg_frame_rate", ""))
            or _parse_fraction(stream.get("r_frame_rate", "")))
    frames = as_int("nb_frames")

    # Duration, most to least trustworthy: stream -> container -> frames / rate.
    duration = as_float("duration")
    if duration is None:
        try:
            duration = float(payload["format"]["duration"])
        except (KeyError, TypeError, ValueError):
            duration = None
    if duration is None and frames and rate:
        duration = frames / rate

    return VideoInfo(
        width=as_int("width"),
        height=as_int("height"),
        r_frame_rate=_parse_fraction(stream.get("r_frame_rate", "")),
        avg_frame_rate=rate,
        nb_frames=frames,
        duration=duration,
        codec=stream.get("codec_name"),
    )


def measure_flicker(source: str, size: int = 5, timeout: float = 3600.0) -> list[float]:
    """Per-frame relative luminance change, WITHOUT modifying any pixels.

    Uses ffmpeg's deflicker filter in `bypass=1` mode, which was verified to leave the
    output bit-identical to the source. Run this before deciding to spend a deflicker
    pass: if the values are already small, the pass is not worth the time.

    Returns the `lavfi.deflicker.relative_change` series.
    """
    source = _absolute_source(source)
    with tempfile.TemporaryDirectory() as tmp:
        report = Path(tmp) / "flicker.txt"
        # The report filename is passed BARE and ffmpeg is run from the temp directory.
        # An absolute Windows path cannot be used here: ffmpeg's filter parser treats
        # ":" as its option separator, so "file=C:/..." fails to parse. Escaping the
        # colon is possible but fragile across shells; changing directory is not.
        _run(
            [
                "ffmpeg", "-v", "error", "-i", source,
                "-vf", f"deflicker=size={size}:bypass=1,"
                       f"metadata=mode=print:file=flicker.txt",
                "-an", "-f", "null", "-",
            ],
            timeout=timeout, cwd=tmp,
        )
        if not report.exists():
            return []
        values: list[float] = []
        for line in report.read_text(encoding="utf-8", errors="replace").splitlines():
            if "relative_change=" in line:
                try:
                    values.append(float(line.split("relative_change=")[1].strip()))
                except (IndexError, ValueError):
                    continue
        return values


def flicker_summary(values: list[float]) -> dict:
    """Reduce a relative_change series to a decision-ready summary."""
    if not values:
        return {"frames": 0, "mean_abs": None, "p95_abs": None, "worth_deflickering": None}
    magnitudes = sorted(abs(v) for v in values)
    mean_abs = sum(magnitudes) / len(magnitudes)
    p95 = magnitudes[min(len(magnitudes) - 1, int(0.95 * len(magnitudes)))]
    return {
        "frames": len(values),
        "mean_abs": mean_abs,
        "p95_abs": p95,
        # 1% mean frame-to-frame luminance change is around where flicker becomes
        # visible on a large display. A threshold, not a law -- look at the footage.
        "worth_deflickering": mean_abs > 0.01,
    }


def detect_shots(source: str, threshold: float = 0.3, timeout: float = 3600.0) -> list[float]:
    """Shot-boundary timestamps in seconds.

    Both the deflicker pass and the temporal chroma median assume their window lies
    inside one continuous shot. Across a hard cut, deflicker reads a legitimate
    brightness change as flicker and pumps against it, and the median blends colour from
    two unrelated scenes. Make shot boundaries chunk boundaries where you can.

    An empty list is a meaningful answer for a single-take film, not a failure.
    """
    if not 0.0 < threshold < 1.0:
        raise ValueError("threshold must be in (0, 1)")
    source = _absolute_source(source)
    with tempfile.TemporaryDirectory() as tmp:
        report = Path(tmp) / "shots.txt"
        # Bare filename + cwd, for the same colon-parsing reason as measure_flicker.
        _run(
            [
                "ffmpeg", "-v", "error", "-i", source,
                "-vf", f"select='gt(scene,{threshold})',"
                       f"metadata=mode=print:file=shots.txt",
                "-an", "-f", "null", "-",
            ],
            timeout=timeout, cwd=tmp,
        )
        if not report.exists():
            return []
        times: list[float] = []
        for line in report.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith("frame:") and "pts_time:" in line:
                try:
                    times.append(float(line.split("pts_time:")[1].split()[0]))
                except (IndexError, ValueError):
                    continue
        return times


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------

@dataclass
class ChunkState:
    index: int
    start: float
    duration: float
    status: str = STATUS_PENDING
    stages: dict = field(default_factory=dict)
    error: str | None = None
    updated: float = 0.0


class Manifest:
    """Resumable per-chunk state, persisted as JSON.

    Written atomically (temp file + os.replace) so a session killed mid-write cannot
    leave a truncated manifest -- which would otherwise lose the record of every
    completed chunk, the one thing this file exists to protect.
    """

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.data: dict = {"version": 1, "source": None, "chunks": {}, "config": {}}
        if self.path.exists():
            self.load()

    def load(self) -> None:
        try:
            self.data = json.loads(self.path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise PipelineError(
                f"manifest at {self.path} is corrupt: {exc}. "
                "Move it aside to restart, or repair it by hand -- do not delete it "
                "without checking which chunks it says are done."
            ) from exc
        self.data.setdefault("chunks", {})
        self.data.setdefault("config", {})

    def save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        handle, tmp_name = tempfile.mkstemp(
            dir=str(self.path.parent), prefix=".manifest-", suffix=".tmp"
        )
        try:
            with os.fdopen(handle, "w", encoding="utf-8") as stream:
                json.dump(self.data, stream, indent=2)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(tmp_name, self.path)
        except BaseException:
            # Never leave a stray temp file behind on failure.
            try:
                os.unlink(tmp_name)
            except OSError:
                pass
            raise

    # -- chunk access ------------------------------------------------------

    def set_chunks(self, chunks: list[ChunkState]) -> None:
        self.data["chunks"] = {str(c.index): asdict(c) for c in chunks}
        self.save()

    def get(self, index: int) -> ChunkState | None:
        raw = self.data["chunks"].get(str(index))
        return ChunkState(**raw) if raw else None

    def all_chunks(self) -> list[ChunkState]:
        return [
            ChunkState(**raw)
            for _, raw in sorted(self.data["chunks"].items(), key=lambda kv: int(kv[0]))
        ]

    def mark_stage(self, index: int, stage: str, output: str) -> None:
        """Record that one stage finished for one chunk."""
        key = str(index)
        if key not in self.data["chunks"]:
            raise KeyError(f"chunk {index} is not in the manifest")
        entry = self.data["chunks"][key]
        entry.setdefault("stages", {})[stage] = output
        entry["updated"] = time.time()
        entry["error"] = None
        self.save()

    def mark_done(self, index: int) -> None:
        key = str(index)
        if key not in self.data["chunks"]:
            raise KeyError(f"chunk {index} is not in the manifest")
        self.data["chunks"][key]["status"] = STATUS_DONE
        self.data["chunks"][key]["updated"] = time.time()
        self.save()

    def mark_failed(self, index: int, error: str) -> None:
        key = str(index)
        if key not in self.data["chunks"]:
            raise KeyError(f"chunk {index} is not in the manifest")
        self.data["chunks"][key]["status"] = STATUS_FAILED
        self.data["chunks"][key]["error"] = str(error)[:2000]
        self.data["chunks"][key]["updated"] = time.time()
        self.save()

    def pending(self) -> list[ChunkState]:
        """Chunks still to do. Failed chunks are retried; done ones are skipped."""
        return [c for c in self.all_chunks() if c.status != STATUS_DONE]

    def stage_done(self, index: int, stage: str) -> str | None:
        """The output path for a completed stage, if it exists on disk.

        A manifest entry whose file has vanished (a non-persistent /kaggle/temp, say)
        must NOT count as done -- otherwise assembly silently produces a short film.
        """
        chunk = self.get(index)
        if not chunk:
            return None
        output = chunk.stages.get(stage)
        if output and Path(output).exists():
            return output
        return None

    def progress(self) -> dict:
        chunks = self.all_chunks()
        done = sum(1 for c in chunks if c.status == STATUS_DONE)
        failed = sum(1 for c in chunks if c.status == STATUS_FAILED)
        return {
            "total": len(chunks),
            "done": done,
            "failed": failed,
            "pending": len(chunks) - done,
            "percent": (100.0 * done / len(chunks)) if chunks else 0.0,
        }


# ---------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------

def plan_chunks(
    duration: float, chunk_seconds: float = DEFAULT_CHUNK_SECONDS,
    shot_boundaries: list[float] | None = None,
) -> list[ChunkState]:
    """Split a duration into chunks, snapping to shot boundaries where possible.

    Snapping matters because deflicker and the temporal median both misbehave across a
    hard cut. A boundary is used when it falls within +/-50% of the nominal chunk length,
    so chunks stay roughly uniform rather than degenerating to one-per-shot.
    """
    if duration <= 0:
        raise ValueError(f"duration must be positive, got {duration}")
    if chunk_seconds <= 0:
        raise ValueError(f"chunk_seconds must be positive, got {chunk_seconds}")

    boundaries = sorted(t for t in (shot_boundaries or []) if 0 < t < duration)
    chunks: list[ChunkState] = []
    start = 0.0
    index = 0

    while start < duration - 1e-6:
        target = start + chunk_seconds
        if target >= duration - 1e-6:
            end = duration
        else:
            window_low = start + chunk_seconds * 0.5
            window_high = start + chunk_seconds * 1.5
            candidates = [t for t in boundaries if window_low <= t <= window_high]
            end = min(candidates, key=lambda t: abs(t - target)) if candidates else target
        chunks.append(ChunkState(index=index, start=start, duration=end - start))
        start = end
        index += 1

    # Absorb a runt tail into the previous chunk. Container durations rarely divide
    # evenly, so the last chunk is often a fraction of a second -- and seeking to a
    # point that close to EOF yields a truncated, undecodable file rather than an
    # empty one, which a size check does not catch.
    if len(chunks) > 1 and chunks[-1].duration < chunk_seconds * MIN_CHUNK_FRACTION:
        runt = chunks.pop()
        chunks[-1].duration += runt.duration

    return chunks


def extract_chunk(
    source: str, output: str | Path, start: float, duration: float,
    codec: str = LOSSLESS_CODEC, extra_filters: str | None = None,
    timeout: float = 7200.0,
) -> Path:
    """Extract one chunk to a lossless file.

    `-ss` before `-i` seeks first, so for an http(s) source ffmpeg issues a byte-range
    request rather than reading the whole file. That is what makes an 85 GB master
    workable on a platform with 20 GB of disk.
    """
    if duration <= 0:
        raise ValueError(f"duration must be positive, got {duration}")
    if start < 0:
        raise ValueError(f"start must be non-negative, got {start}")

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    cmd = ["ffmpeg", "-v", "error", "-y", "-ss", f"{start:.3f}", "-i", source,
           "-t", f"{duration:.3f}"]
    if extra_filters:
        cmd += ["-vf", extra_filters]
    cmd += ["-an", "-c:v", codec]
    if codec == "ffv1":
        cmd += ["-level", "3"]
    cmd += [str(output)]

    _run(cmd, timeout=timeout)
    if not output.exists() or output.stat().st_size == 0:
        raise PipelineError(f"chunk extraction produced no output at {output}")

    # A non-empty file is not the same as a usable one: a seek very close to EOF can
    # write a header with no decodable frames. Probing here fails loudly at extraction
    # rather than three stages later during assembly.
    try:
        info = probe_video(str(output), timeout=120.0)
    except PipelineError as exc:
        raise PipelineError(
            f"chunk at {output} was written but does not decode "
            f"(start={start:.3f}s duration={duration:.3f}s): {exc}"
        ) from exc
    if not info.width or not info.height:
        raise PipelineError(f"chunk at {output} has no usable video stream")
    return output


def deflicker_chunk(
    source: str | Path, output: str | Path, size: int = 7, mode: str = "am",
    codec: str = LOSSLESS_CODEC, timeout: float = 7200.0,
) -> Path:
    """Global luminance deflicker, before colorization.

    `mode` defaults to "am" on measured evidence, not on argument. Across three synthetic
    flicker profiles (periodic, log-normal, isolated spikes) and eight seeds, `am` and
    `gm` both cut the flicker metric by 78-88% while `median` managed -0.5% to 6.5% --
    including on the isolated-spike profile where outlier rejection should have been its
    strongest case. The plausible "median rejects blown frames" argument does not survive
    measurement. See research/pipeline-optimization.md section 2.3.

    This corrects GLOBAL exposure pumping only: the filter computes one scalar per frame
    and applies a single multiplicative gain, so spatially-varying density flicker passes
    straight through. It touches only the luma plane, so it is safe before colorization.
    """
    valid_modes = {"am", "gm", "hm", "qm", "cm", "pm", "median"}
    if mode not in valid_modes:
        raise ValueError(f"mode must be one of {sorted(valid_modes)}, got {mode!r}")
    if not 2 <= size <= 129:
        raise ValueError(f"size must be in [2, 129], got {size}")

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["ffmpeg", "-v", "error", "-y", "-i", str(source),
           "-vf", f"deflicker=size={size}:mode={mode}", "-an", "-c:v", codec]
    if codec == "ffv1":
        cmd += ["-level", "3"]
    cmd += [str(output)]
    _run(cmd, timeout=timeout)
    if not output.exists() or output.stat().st_size == 0:
        raise PipelineError(f"deflicker produced no output at {output}")
    return output


def decimate_chunk(
    source: str | Path, output: str | Path, frame_rate: float | None = None,
    codec: str = LOSSLESS_CODEC, timeout: float = 7200.0,
) -> Path:
    """Drop duplicated frames before the GPU stages.

    Scans of hand-cranked film are routinely padded to a modern playback rate by
    repeating frames. The Market Street 4K master says so in its own archive.org
    description: it was *"scanned to 16fps ... every second frame is doubled for a
    playback framerate of 24fps."* Two source frames become three, so **one frame in
    three is a duplicate**, and colorizing or upscaling those is pure waste.

    Measured on a clip built to that exact 16->24 recipe: 48 frames in, 32 out, matching
    the unique source count exactly -- a 33.3% reduction in GPU work.

    Run this BEFORE colorization and super-resolution, and re-time afterwards to the
    film's true rate. `setpts` restamps the survivors so the output is not left with
    gaps where the dropped frames were.
    """
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    # hi/lo/frac are mpdecimate's defaults, stated explicitly so the behaviour is
    # visible and tunable rather than implicit. A doubled frame is bit-identical, so
    # even a strict threshold catches it; these values also tolerate re-encoded sources
    # where duplicates differ by a little compression noise.
    decimate = "mpdecimate=hi=64*12:lo=64*5:frac=0.1"
    if frame_rate:
        # Space the survivors at the TARGET rate inside the filter graph, then tell the
        # muxer to pass them through untouched. Doing this with an output `-r` instead
        # makes ffmpeg re-time a stream whose timestamps still describe the old rate,
        # and it drops real frames to fit -- 48 doubled frames came out as 23 rather
        # than the correct 32 before this was fixed.
        filters = f"{decimate},setpts=N/({frame_rate:.6f}*TB)"
    else:
        filters = f"{decimate},setpts=N/FRAME_RATE/TB"

    cmd = ["ffmpeg", "-v", "error", "-y", "-i", str(source), "-vf", filters,
           "-fps_mode", "passthrough", "-an", "-c:v", codec]
    if codec == "ffv1":
        cmd += ["-level", "3"]
    cmd += [str(output)]
    _run(cmd, timeout=timeout)

    if not output.exists() or output.stat().st_size == 0:
        raise PipelineError(f"decimate produced no output at {output}")
    return output


def count_frames(source: str | Path, timeout: float = 3600.0) -> int:
    """Exact decoded frame count. Slower than reading nb_frames, but always right."""
    proc = _run(
        ["ffprobe", "-v", "error", "-select_streams", "v:0", "-count_frames",
         "-show_entries", "stream=nb_read_frames", "-of", "csv=p=0", str(source)],
        timeout=timeout,
    )
    text = proc.stdout.decode("utf-8").strip().split(",")[0]
    try:
        return int(text)
    except ValueError as exc:
        raise PipelineError(f"could not read frame count from {source!r}: {text!r}") from exc


def assemble(
    chunk_paths: list[str | Path], output: str | Path,
    frame_rate: float | None = None, crf: int = 16, preset: str = "slow",
    pix_fmt: str = "yuv420p", timeout: float = 14400.0,
) -> Path:
    """Concatenate processed chunks into the delivery file.

    Refuses to run on a missing chunk rather than silently producing a short film -- the
    failure mode this whole manifest exists to prevent.
    """
    if not chunk_paths:
        raise ValueError("no chunks to assemble")
    missing = [str(p) for p in chunk_paths if not Path(p).exists()]
    if missing:
        raise PipelineError(
            f"{len(missing)} chunk(s) missing, refusing to assemble a short film: "
            + ", ".join(missing[:5]) + (" ..." if len(missing) > 5 else "")
        )

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        listing = Path(tmp) / "concat.txt"
        listing.write_text(
            "".join(
                "file '{}'\n".format(Path(p).resolve().as_posix().replace("'", r"'\''"))
                for p in chunk_paths
            ),
            encoding="utf-8",
        )
        cmd = ["ffmpeg", "-v", "error", "-y", "-f", "concat", "-safe", "0",
               "-i", str(listing)]
        if frame_rate:
            cmd += ["-r", f"{frame_rate:.6f}"]
        cmd += ["-c:v", "libx264", "-crf", str(crf), "-preset", preset,
                "-pix_fmt", pix_fmt, "-movflags", "+faststart", str(output)]
        _run(cmd, timeout=timeout)

    if not output.exists() or output.stat().st_size == 0:
        raise PipelineError(f"assembly produced no output at {output}")
    return output


def make_comparison(
    before: str | Path, after: str | Path, output: str | Path,
    mode: str = "sidebyside", timeout: float = 7200.0,
) -> Path:
    """Render a before/after comparison.

    Worth building into the pipeline rather than doing by hand: it makes the
    transformation visible to a YouTube reviewer within the first seconds, which is the
    reused-content policy's "viewers can tell that there's a meaningful difference" test.
    It is also the short-form hook.

    `mode` is "sidebyside" (hstack) or "wipe" (left half before, right half after).
    """
    if mode not in ("sidebyside", "wipe"):
        raise ValueError(f"mode must be 'sidebyside' or 'wipe', got {mode!r}")
    for path in (before, after):
        if not Path(path).exists():
            raise PipelineError(f"missing input for comparison: {path}")

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    if mode == "sidebyside":
        filtergraph = "[0:v][1:v]hstack=inputs=2[v]"
    else:
        # Take the left half of `before` and the right half of `after`, then stack.
        filtergraph = (
            "[0:v]crop=iw/2:ih:0:0[l];"
            "[1:v]crop=iw/2:ih:iw/2:0[r];"
            "[l][r]hstack=inputs=2[v]"
        )

    _run(
        ["ffmpeg", "-v", "error", "-y", "-i", str(before), "-i", str(after),
         "-filter_complex", filtergraph, "-map", "[v]",
         "-c:v", "libx264", "-crf", "18", "-preset", "medium",
         "-pix_fmt", "yuv420p", str(output)],
        timeout=timeout,
    )
    if not output.exists() or output.stat().st_size == 0:
        raise PipelineError(f"comparison produced no output at {output}")
    return output


def disk_report(path: str | Path = ".") -> dict:
    """Free space, so a job can refuse to start rather than dying at 90%."""
    usage = shutil.disk_usage(str(path))
    return {
        "total_gb": usage.total / 1e9,
        "used_gb": usage.used / 1e9,
        "free_gb": usage.free / 1e9,
    }
