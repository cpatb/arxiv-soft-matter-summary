"""End-to-end verification of the non-GPU pipeline against a real synthetic clip.

This builds an actual video with ffmpeg, segments it, deflickers it, assembles it, and
checks the manifest survives a simulated session kill. Nothing here is mocked -- the
point is to exercise the real ffmpeg invocations, because that is where the bugs are.

Run:  python test_pipeline.py
"""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

import pipeline as pl

WIDTH, HEIGHT = 192, 144
FPS = 16
SECONDS = 6
N_FRAMES = FPS * SECONDS


def build_clip(path: Path, flicker: bool = True, cut_at: int | None = None) -> None:
    """A moving box, optional brightness flicker, optional hard cut."""
    rng = np.random.default_rng(4)
    frames = []
    for i in range(N_FRAMES):
        if cut_at is not None and i >= cut_at:
            base = np.full((HEIGHT, WIDTH), 40, np.uint8)      # different scene
            base[20:60, 100:150] = 200
        else:
            base = np.full((HEIGHT, WIDTH), 120, np.uint8)
            left = 10 + (i * 2) % (WIDTH - 50)
            base[50:100, left:left + 40] = 210
        gain = float(np.exp(rng.normal(0, 0.10))) if flicker else 1.0
        frames.append(np.clip(base.astype(np.float32) * gain, 0, 255).astype(np.uint8))

    raw = b"".join(np.repeat(f[:, :, None], 3, axis=2).tobytes() for f in frames)
    subprocess.run(
        ["ffmpeg", "-y", "-v", "error", "-f", "rawvideo", "-pix_fmt", "rgb24",
         "-s", f"{WIDTH}x{HEIGHT}", "-r", str(FPS), "-i", "-",
         "-c:v", "ffv1", "-level", "3", str(path)],
        input=raw, check=True, capture_output=True,
    )


# ---------------------------------------------------------------------------

def test_plan_chunks_covers_duration_exactly():
    chunks = pl.plan_chunks(100.0, 30.0)
    assert len(chunks) == 4, len(chunks)
    total = sum(c.duration for c in chunks)
    assert abs(total - 100.0) < 1e-6, total
    assert abs(chunks[0].start) < 1e-9
    # No gaps and no overlaps.
    for previous, following in zip(chunks, chunks[1:]):
        assert abs((previous.start + previous.duration) - following.start) < 1e-9
    assert abs(chunks[-1].start + chunks[-1].duration - 100.0) < 1e-6
    print(f"  4 chunks, total {total:.1f}s, contiguous")


def test_plan_chunks_snaps_to_shots():
    # A boundary at 34 s is inside the +/-50% window around a 30 s target.
    chunks = pl.plan_chunks(100.0, 30.0, shot_boundaries=[34.0, 71.0])
    starts = [round(c.start, 3) for c in chunks]
    assert 34.0 in starts, starts
    assert 71.0 in starts, starts
    assert abs(sum(c.duration for c in chunks) - 100.0) < 1e-6
    # A boundary far outside the window must be ignored, not forced.
    ignored = pl.plan_chunks(100.0, 30.0, shot_boundaries=[2.0])
    assert round(ignored[0].duration, 3) == 30.0, ignored[0].duration
    print(f"  snapped to shots at {starts}")


def test_plan_chunks_absorbs_runt_tail():
    """A sub-second trailing chunk must be merged, not left to decode as garbage."""
    chunks = pl.plan_chunks(6.0625, 2.0)      # would otherwise leave a 0.0625 s tail
    assert len(chunks) == 3, [round(c.duration, 4) for c in chunks]
    assert abs(sum(c.duration for c in chunks) - 6.0625) < 1e-6
    assert chunks[-1].duration > 2.0, chunks[-1].duration

    # A tail above the threshold is kept as its own chunk.
    kept = pl.plan_chunks(100.0, 30.0)        # 10 s tail, above 30 * 0.25
    assert len(kept) == 4, len(kept)

    # Indices must stay contiguous from zero after a merge.
    assert [c.index for c in chunks] == [0, 1, 2]
    print(f"  runt absorbed: durations {[round(c.duration, 4) for c in chunks]}")


def test_plan_chunks_rejects_bad_input():
    for duration, size in [(0, 30), (-5, 30), (100, 0), (100, -1)]:
        try:
            pl.plan_chunks(duration, size)
        except ValueError:
            pass
        else:
            raise AssertionError(f"accepted duration={duration} chunk={size}")
    print("  bad durations rejected")


def test_manifest_roundtrip_and_resume(tmp: Path):
    manifest_path = tmp / "manifest.json"
    manifest = pl.Manifest(manifest_path)
    manifest.set_chunks(pl.plan_chunks(90.0, 30.0))

    assert manifest.progress()["total"] == 3
    assert manifest.progress()["done"] == 0

    stage_output = tmp / "chunk0.mkv"
    stage_output.write_bytes(b"x")
    manifest.mark_stage(0, "deflicker", str(stage_output))
    manifest.mark_done(0)
    manifest.mark_failed(1, "boom")

    # Reopen from disk, as a fresh session would.
    reopened = pl.Manifest(manifest_path)
    assert reopened.progress()["done"] == 1
    assert reopened.progress()["failed"] == 1
    assert [c.index for c in reopened.pending()] == [1, 2]
    assert reopened.stage_done(0, "deflicker") == str(stage_output)

    # A stage whose file has vanished must NOT count as done.
    stage_output.unlink()
    assert reopened.stage_done(0, "deflicker") is None, (
        "a manifest entry pointing at a deleted file must not read as complete"
    )
    print("  manifest resumes; vanished stage output correctly not trusted")


def test_manifest_write_is_atomic(tmp: Path):
    manifest_path = tmp / "atomic.json"
    manifest = pl.Manifest(manifest_path)
    manifest.set_chunks(pl.plan_chunks(60.0, 30.0))
    manifest.mark_done(0)

    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert payload["chunks"]["0"]["status"] == pl.STATUS_DONE
    # No temp files left behind.
    strays = [p.name for p in manifest_path.parent.glob(".manifest-*")]
    assert not strays, f"stray temp files: {strays}"
    print("  atomic write leaves no temp files")


def test_manifest_rejects_corruption(tmp: Path):
    bad = tmp / "corrupt.json"
    bad.write_text("{not json", encoding="utf-8")
    try:
        pl.Manifest(bad)
    except pl.PipelineError as exc:
        assert "corrupt" in str(exc)
    else:
        raise AssertionError("corrupt manifest was accepted")
    print("  corrupt manifest refused loudly")


def test_manifest_rejects_unknown_chunk(tmp: Path):
    manifest = pl.Manifest(tmp / "m2.json")
    manifest.set_chunks(pl.plan_chunks(30.0, 30.0))
    for call in (lambda: manifest.mark_done(99),
                 lambda: manifest.mark_failed(99, "x"),
                 lambda: manifest.mark_stage(99, "s", "o")):
        try:
            call()
        except KeyError:
            pass
        else:
            raise AssertionError("unknown chunk index accepted")
    print("  unknown chunk index refused")


def test_probe_reads_real_video(tmp: Path):
    clip = tmp / "probe.mkv"
    build_clip(clip)
    info = pl.probe_video(str(clip))
    assert info.width == WIDTH and info.height == HEIGHT, (info.width, info.height)
    assert info.avg_frame_rate is not None
    assert abs(info.avg_frame_rate - FPS) < 0.5, info.avg_frame_rate
    assert info.codec == "ffv1", info.codec
    print(f"  probed {info.width}x{info.height} @ {info.avg_frame_rate:.2f} fps, {info.codec}")


def test_flicker_probe_detects_and_summarises(tmp: Path):
    noisy, clean = tmp / "noisy.mkv", tmp / "clean.mkv"
    build_clip(noisy, flicker=True)
    build_clip(clean, flicker=False)

    noisy_summary = pl.flicker_summary(pl.measure_flicker(str(noisy)))
    clean_summary = pl.flicker_summary(pl.measure_flicker(str(clean)))

    assert noisy_summary["frames"] > 0, "probe returned no frames"
    assert noisy_summary["mean_abs"] > clean_summary["mean_abs"], (
        f"flickering clip {noisy_summary['mean_abs']:.4f} should exceed "
        f"clean {clean_summary['mean_abs']:.4f}"
    )
    assert noisy_summary["worth_deflickering"] is True
    print(f"  flickering mean|rel| {noisy_summary['mean_abs']:.4f} "
          f"vs clean {clean_summary['mean_abs']:.4f}")


def test_deflicker_actually_reduces_flicker(tmp: Path):
    noisy, fixed = tmp / "n.mkv", tmp / "f.mkv"
    build_clip(noisy, flicker=True)
    pl.deflicker_chunk(noisy, fixed, size=7, mode="am")

    before = pl.flicker_summary(pl.measure_flicker(str(noisy)))["mean_abs"]
    after = pl.flicker_summary(pl.measure_flicker(str(fixed)))["mean_abs"]
    reduction = 100.0 * (1 - after / before)
    print(f"  mean|rel| {before:.4f} -> {after:.4f} ({reduction:.1f}% reduction, mode=am)")
    assert after < before, "deflicker made it worse"
    assert reduction > 40.0, f"only {reduction:.1f}% reduction; expected mode=am to do better"


def test_deflicker_rejects_bad_parameters(tmp: Path):
    clip = tmp / "p.mkv"
    build_clip(clip)
    for kwargs in [{"mode": "nonsense"}, {"size": 1}, {"size": 500}]:
        try:
            pl.deflicker_chunk(clip, tmp / "out.mkv", **kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError(f"accepted {kwargs}")
    print("  bad deflicker parameters rejected")


def test_full_chunked_roundtrip(tmp: Path):
    """Segment -> deflicker -> assemble, with a simulated session kill in the middle."""
    source = tmp / "source.mkv"
    build_clip(source, flicker=True)

    info = pl.probe_video(str(source))
    assert info.duration is not None
    manifest = pl.Manifest(tmp / "run.json")
    chunks = pl.plan_chunks(info.duration, chunk_seconds=2.0)
    manifest.set_chunks(chunks)
    assert len(chunks) >= 3, len(chunks)

    work = tmp / "work"

    def process(limit=None):
        """Process pending chunks, optionally stopping early to simulate a kill."""
        processed = 0
        for chunk in manifest.pending():
            if limit is not None and processed >= limit:
                break
            raw = pl.extract_chunk(str(source), work / f"raw{chunk.index:04d}.mkv",
                                   chunk.start, chunk.duration)
            manifest.mark_stage(chunk.index, "extract", str(raw))
            fixed = pl.deflicker_chunk(raw, work / f"df{chunk.index:04d}.mkv", mode="am")
            manifest.mark_stage(chunk.index, "deflicker", str(fixed))
            manifest.mark_done(chunk.index)
            processed += 1
        return processed

    first = process(limit=2)
    assert first == 2
    mid = manifest.progress()
    assert mid["done"] == 2, mid

    # A new session reopens the manifest and must skip completed work.
    resumed = pl.Manifest(tmp / "run.json")
    assert resumed.progress()["done"] == 2
    remaining = len(resumed.pending())
    assert remaining == len(chunks) - 2

    process()
    final = manifest.progress()
    assert final["done"] == len(chunks), final
    assert final["percent"] == 100.0

    outputs = [manifest.stage_done(c.index, "deflicker") for c in manifest.all_chunks()]
    assert all(outputs), "a stage output went missing"

    delivered = pl.assemble(outputs, tmp / "final.mp4", frame_rate=FPS)
    out_info = pl.probe_video(str(delivered))
    assert out_info.duration is not None
    assert abs(out_info.duration - info.duration) < 0.5, (
        f"assembled duration {out_info.duration:.2f}s != source {info.duration:.2f}s"
    )
    assert out_info.width == WIDTH and out_info.height == HEIGHT
    print(f"  {len(chunks)} chunks, killed after 2, resumed, assembled "
          f"{out_info.duration:.2f}s (source {info.duration:.2f}s)")


def test_decimate_removes_doubled_frames(tmp: Path):
    """Reproduce the archive's own 16->24 fps doubling and check it is undone exactly."""
    unique, doubled = tmp / "u.mkv", tmp / "d.mkv"
    rng = np.random.default_rng(2)
    frames = []
    for i in range(32):
        f = np.full((HEIGHT, WIDTH), 100, np.uint8)
        f[40:90, 5 + i * 4:35 + i * 4] = 215
        frames.append(f)

    emitted = []
    for i, f in enumerate(frames):
        emitted.append(f)
        if i % 2 == 1:                     # double every second frame
            emitted.append(f.copy())
    assert len(emitted) == 48, len(emitted)

    for path, seq, fps in ((unique, frames, 16), (doubled, emitted, 24)):
        raw = b"".join(np.repeat(f[:, :, None], 3, axis=2).tobytes() for f in seq)
        subprocess.run(
            ["ffmpeg", "-y", "-v", "error", "-f", "rawvideo", "-pix_fmt", "rgb24",
             "-s", f"{WIDTH}x{HEIGHT}", "-r", str(fps), "-i", "-",
             "-c:v", "ffv1", "-level", "3", str(path)],
            input=raw, check=True, capture_output=True)

    before = pl.count_frames(doubled)
    out = pl.decimate_chunk(doubled, tmp / "dec.mkv", frame_rate=16)
    after = pl.count_frames(out)
    saved = 100.0 * (before - after) / before
    print(f"  {before} frames -> {after} ({saved:.1f}% of GPU work saved); "
          f"unique source was {len(frames)}")

    assert before == 48, before
    assert after == len(frames), f"expected {len(frames)} unique frames, got {after}"
    assert 30.0 < saved < 36.0, f"expected ~33% saving, got {saved:.1f}%"

    # Decimating an already-unique clip must not throw frames away.
    clean = pl.decimate_chunk(unique, tmp / "clean.mkv", frame_rate=16)
    assert pl.count_frames(clean) == len(frames), "decimate dropped non-duplicate frames"
    _ = rng


def test_assemble_refuses_missing_chunk(tmp: Path):
    clip = tmp / "one.mkv"
    build_clip(clip)
    try:
        pl.assemble([clip, tmp / "does_not_exist.mkv"], tmp / "bad.mp4")
    except pl.PipelineError as exc:
        assert "missing" in str(exc).lower()
    else:
        raise AssertionError("assembled a short film instead of failing")
    print("  assembly refuses to silently produce a short film")


def test_comparison_renders(tmp: Path):
    before, after = tmp / "b.mkv", tmp / "a.mkv"
    build_clip(before, flicker=True)
    build_clip(after, flicker=False)

    side = pl.make_comparison(before, after, tmp / "side.mp4", mode="sidebyside")
    info = pl.probe_video(str(side))
    assert info.width == WIDTH * 2, info.width

    wipe = pl.make_comparison(before, after, tmp / "wipe.mp4", mode="wipe")
    wipe_info = pl.probe_video(str(wipe))
    assert wipe_info.width == WIDTH, wipe_info.width
    print(f"  side-by-side {info.width}px, wipe {wipe_info.width}px")


def test_shot_detection_finds_a_cut(tmp: Path):
    with_cut, without = tmp / "cut.mkv", tmp / "nocut.mkv"
    build_clip(with_cut, flicker=False, cut_at=N_FRAMES // 2)
    build_clip(without, flicker=False)

    cuts = pl.detect_shots(str(with_cut), threshold=0.3)
    none = pl.detect_shots(str(without), threshold=0.3)
    print(f"  clip with a cut: {len(cuts)} boundaries at {[round(t,2) for t in cuts]}; "
          f"single-take clip: {len(none)}")
    assert len(cuts) >= 1, "failed to find an obvious hard cut"
    assert len(none) < len(cuts), "found more cuts in the single-take clip"


def main() -> int:
    if not (pl.have_tool("ffmpeg") and pl.have_tool("ffprobe")):
        print("SKIP: ffmpeg/ffprobe not on PATH")
        return 0

    needs_tmp = {
        "manifest roundtrip + resume": test_manifest_roundtrip_and_resume,
        "manifest atomic write": test_manifest_write_is_atomic,
        "manifest rejects corruption": test_manifest_rejects_corruption,
        "manifest rejects unknown chunk": test_manifest_rejects_unknown_chunk,
        "probe reads real video": test_probe_reads_real_video,
        "flicker probe": test_flicker_probe_detects_and_summarises,
        "deflicker reduces flicker": test_deflicker_actually_reduces_flicker,
        "deflicker validates parameters": test_deflicker_rejects_bad_parameters,
        "full chunked roundtrip + resume": test_full_chunked_roundtrip,
        "decimate removes doubled frames": test_decimate_removes_doubled_frames,
        "assembly refuses missing chunk": test_assemble_refuses_missing_chunk,
        "comparison renders": test_comparison_renders,
        "shot detection": test_shot_detection_finds_a_cut,
    }
    plain = {
        "plan_chunks covers duration": test_plan_chunks_covers_duration_exactly,
        "plan_chunks snaps to shots": test_plan_chunks_snaps_to_shots,
        "plan_chunks validates input": test_plan_chunks_rejects_bad_input,
        "plan_chunks absorbs runt tail": test_plan_chunks_absorbs_runt_tail,
    }

    failures = 0
    for name, fn in plain.items():
        print(f"\n[{name}]")
        try:
            fn()
            print("  PASS")
        except AssertionError as exc:
            failures += 1
            print(f"  FAIL: {exc}")

    root = Path(tempfile.mkdtemp(prefix="stc-test-"))
    try:
        for name, fn in needs_tmp.items():
            print(f"\n[{name}]")
            case_dir = root / name.replace(" ", "_").replace("+", "")
            case_dir.mkdir(parents=True, exist_ok=True)
            try:
                fn(case_dir)
                print("  PASS")
            except (AssertionError, pl.PipelineError) as exc:
                failures += 1
                print(f"  FAIL: {exc}")
    finally:
        shutil.rmtree(root, ignore_errors=True)

    total = len(plain) + len(needs_tmp)
    print(f"\n{total - failures}/{total} passed")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
