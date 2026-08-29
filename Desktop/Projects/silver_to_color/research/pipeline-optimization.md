# Pipeline optimization

Research dates **2026-08-25/26**. Model claims come from repository source, READMEs, GitHub
issues, and papers, never from blog posts. Where no figure is published, this document says so
and labels any estimate as an estimate, with the arithmetic shown.

Two things in here were **measured locally rather than taken on trust**, and one of them
overturns a plausible-sounding recommendation. Those are marked **[MEASURED]**.

Tags: **[QUOTED]**: verbatim from a primary source · **[ESTIMATE]**: derived, arithmetic shown ·
**[NONE]**: no published figure exists · **[MEASURED]** measured on this machine ·
**[TASTE]**: a creative decision, not a technical one.

---

## 0. Two findings that change the shape of the pipeline

### 0.1 The source file does not fit on the platform

Queried live from archive.org on 2026-08-26 via `tools/ia_meta.py`:

| File | Resolution | Size | Duration |
|---|---|---|---|
| `Market_Street_1906_4K_scan_20181016.mov` | **4096x3072** | **85.4 GB** | 404.5 s |
| `Market_Street_1906_4K_scan_20181016.hevc.mp4` | not stated | 14.4 GB | not stated |
| `MarketStreet_4K_to_2K_cropped_higher_contrast.mov` | 2048x1536 | 2.5 GB | 716.8 s |
| `*.mp4` / `*.ogv` derivatives | 640x480 / 400x300 | ~60 MB / 44 MB | 716.8 s |

Item total: **102.4 GB**, 61 files, collection `prelinger`, uploader `footage@panix.com`.

Three consequences, none of which the current plan accounts for:

1. **The master is 4096x3072, not 3840x2160.** That is a 4:3 frame at 12.6 MP, 52% more pixels
   than UHD. Every VRAM and throughput estimate scaled from "4K" is 1.5x optimistic.
2. **85.4 GB will not fit in Kaggle's `/kaggle/working` (~20 GB) and only barely in `/kaggle/temp`
   (~60 GB), which does not persist.** You cannot download this file and work on it. The pipeline
   has to segment straight from the URL, or use the 14.4 GB HEVC, or the 2.5 GB 2K version.
   This is a hard constraint, not a tuning issue.
3. **The duration mismatch is a flag, not a conclusion.** The 4K master reports 404.5 s while the
   2K version reports 716.8 s, giving a ratio of 1.772. If the master is at 30 fps, 404.5 s is 12,135
   frames, and those frames over 716.8 s give **~16.9 fps**. That is consistent with the 2K
   version already being speed-corrected while the 4K master sits at scan rate. It is **not**
   consistent with the "24→16 fps" correction in the current plan. **Verify with `ffprobe`
   before rendering anything.** A wrong frame-rate assumption is baked into every frame and is
   expensive to discover late. (See §7.1.)

### 0.2 Super-resolution is the wrong operation on this source

This is the highest-value change available, and it costs nothing.

A 4K scan of 1906 35 mm **already resolves more than the film holds**. The negative stock was
orthochromatic, grain was coarse, lenses were uncoated, and the access copy is several print
generations from a camera negative that in most cases no longer exists. Archival scanning is done
at 4K precisely *because* it oversamples the grain, it captures grain structure, not additional
image detail. There is no latent high-frequency content for an SR network to recover, because the
emulsion never recorded any. What a 4x network reads as "texture to sharpen" is silver-halide
grain and dust.

The corroboration is that **the specialists in this exact domain do not upscale**: Bringing Old
Films Back to Life (Wan et al., CVPR 2022) sets `scale: 1`, and its successor RRTN inherits that.
When researchers build a model for old film, they build a 1:1 restorer.

Map the operations against what is actually wrong with a 1906 print:

| Defect | Right operation | Current pipeline |
|---|---|---|
| Vertical scratches | defect repair / inpainting | **no** |
| Dust, blotches, splice marks | defect repair | **no** |
| Grain, gate noise | temporal denoise | partly, uncontrolled (see §2.2) |
| Gate weave, frame jitter | stabilisation | **no** |
| Density flicker | deflicker | **no** |
| Softness from print generations | modest sharpen | yes, as an SR side effect |
| Missing resolution | super-resolution | **there is none missing** |

**Deliver 2160p anyway, but for a bitrate reason, not a resolution one.** YouTube's recommended
upload bitrates ([support.google.com/youtube/answer/1722171](https://support.google.com/youtube/answer/1722171),
checked 2026-08-26) **[QUOTED]**: 2160p **35-45 Mbps**, 1080p **8 Mbps**. Uploading at 2160p buys
roughly **4-5x the bitrate allocation** for the same picture. Film grain is the hardest thing a
video encoder faces, high-entropy, temporally uncorrelated noise, so grainy 1906 footage at
1080p/8 Mbps turns to mush while the same footage at 2160p/40 Mbps survives. Upload 2160p to buy
bits for the grain. Do not go to 8K: it costs 2-4x more encode, upload and GPU time for content
with no 8K detail, and 80-160 Mbps of bitrate that the picture cannot fill.

---

## 1. Ranked changes

Ordered by (expected gain) ÷ (cost x risk). Sections below give the detail.

| # | Change | Expected gain | Cost | Risk | Detail |
|---|---|---|---|---|---|
| 0 | **Decimate duplicate frames before the GPU stages** | **33.3% of all GPU work, measured** | one ffmpeg pass | **none** | §5.6 |
| 1 | Bilinear instead of nearest chroma upsample | Removes a whole class of chroma jitter | **one word** | **none** | §2.1 |
| 1b | Use the 2K cropped version for the first film | Removes the 85 GB problem entirely | change one URL | none | §5.6 |
| 2 | Stop upscaling; Real-ESRGAN as cleanup at 2160p | ~4x less SR compute; frees ~3-5 GPU-hr/film | one flag | low | §2.2 |
| 3 | Deflicker before colorization, with `mode=am`, **not** `median` | Attacks colour instability at its cause | one ffmpeg pass | low | §2.3 |
| 4 | Motion-compensated temporal chroma median | **3.4x lower chroma error at moving edges [MEASURED]** | ~80 lines, CPU-side | low | §2.4 |
| 5 | Segment from the URL; never download the 85 GB master | Makes the job possible at all on free tier | pipeline restructure | medium | §3 |
| 6 | Process per shot, not per chunk | Stops deflicker and the median fighting hard cuts | shot detection pass | low | §2.5 |
| 7 | Insert RRTN defect repair at ≤1080p | Fixes what viewers actually notice | a whole session to install | **high** | §4.1 |
| 8 | Confirm fp16 and raise batch size | Throughput | small | low | §5 |
| 9 | Replace the median with CMNET2 exemplar propagation | Cure rather than symptom treatment | large | high | §4.2 |
| 10 | RIFE interpolation | | | | **[TASTE]** §6 |

---

## 2. The cheap wins

### 2.1 The nearest-neighbour chroma upsample is a real bug for video

Read DDColor's own pipeline first, because everything else follows from it. From
`ddcolor/pipeline.py` ([raw source](https://raw.githubusercontent.com/piddnad/DDColor/master/ddcolor/pipeline.py),
checked 2026-08-25):

```python
orig_l = cv2.cvtColor(img, cv2.COLOR_BGR2Lab)[:::1]      # FULL RES, never touched
img_resized = cv2.resize(img, (self.input_size, self.input_size))   # 512x512
output_ab = self.model(tensor_gray_rgb).cpu()                # (1, 2, 512, 512)
output_ab_resized = F.interpolate(output_ab, size=(height, width))[0]
output_lab = np.concatenate((orig_l, output_ab_resized), axis=-1)
```

Six facts that drive every later decision:

1. **DDColor outputs chroma only**: `num_output_channels=2`, the a and b of CIE Lab.
2. **Luminance never enters the network.** `orig_l` is taken at full resolution and concatenated
   back at the end. No colorization-side temporal filter can soften a single edge.
3. **Chroma is computed at a fixed 512x512** regardless of source size. Against a 4096x3072
   master that is **2.1% of the source pixel count**.
4. **The chroma upsample is nearest-neighbour.** `F.interpolate` is called with no `mode`, and
   PyTorch's default is `mode="nearest"` (confirmed in `torch/nn/functional.py`).
5. **Frames are processed strictly one at a time**, `.unsqueeze(0)` makes a batch of exactly 1.
6. **The model's only input is luminance**, with a and b zeroed. Output is a deterministic
   function of L.

Fact 4 is a latent bug that matters far more for video than for stills. Nearest-neighbour
upsampling from 512x512 to 4096x3072 replicates each chroma sample into an **8x6 px block** with
a hard edge. As content drifts sub-pixel between frames, the output-pixel-to-source-sample
mapping is a floor operation, so output pixels **snap** between neighbouring chroma samples at
different moments across the frame. That is high-frequency chroma jitter **generated by the
resize itself**, and it will show up in your flicker measurements attributed to DDColor.

Fix:

```python
F.interpolate(output_ab, size=(height, width), mode='bilinear', align_corners=False)
```

Bilinear is also the perceptually correct choice here: chroma is a smooth signal. (Nearest is
mandatory only for label/mask images, where interpolating invents classes that do not exist.)

**Do this before tuning any temporal filter**, or you will be tuning the filter to suppress an
artifact that has a one-word fix. Verified by `tools/test_temporal_chroma.py::test_upsample_is_bilinear_not_nearest`.

### 2.2 Real-ESRGAN: what you are running, and the flag that is quietly wasteful

Confirmed from `inference_realesrgan.py` (checked 2026-08-26):

```python
elif args.model_name == 'realesr-general-x4v3':  # x4 VGG-style model (S size)
    model = SRVGGNetCompact(num_in_ch=3, num_out_ch=3, num_feat=64, num_conv=32,
                            upscale=4, act_type='prelu')
```

**SRVGGNetCompact, 1,213,296 parameters.** Cross-checked: 1,213,296 x 4 bytes fp32 =
4,853,184 B against a released `.pth` of **4,885,111 B**: a 31,927-byte delta that is exactly
state-dict pickle overhead. Plain VGG-style CNN, no attention, no deformable conv, no recurrence,
which is why it runs anywhere.

**Repo status: frozen.** Last push **2024-08-06**, 646 open issues, no maintainer triage. Not
archived. Frozen is fine here: a 1.2 M-parameter CNN with no exotic dependencies does not rot
the way the mmedit-family repos do (§4.3).

**`-s 2` is the waste.** There is no 2x general model. `--outscale 2` runs the **x4 network** and
then downsamples with cv2. From a 4096x3072 input that means materialising a **16384x12288**
tensor and discarding three quarters of it. It works, and it is a defensible way to get
denoise-plus-sharpen, but it is the most expensive possible way to spend a T4 hour on a mild
cleanup.

**VRAM: [NONE] published.** The repo's only guidance is a shrug at line 150: *"If you encounter
CUDA out of memory, try to set --tile with a smaller number."*

**[ESTIMATE] ~5-7 GB fp16 untiled at 3840x2160; ~9-13 GB fp32.** Arithmetic: SRVGGNetCompact runs
all 33 convolutions at *input* resolution and only upsamples in the final PixelShuffle, so
activation memory scales with input pixels x 64 channels. One 64-channel fp16 activation at
8.29 M px is 1.06 GB; ~2 live at once is 2.1 GB; tail conv (48 ch) 0.80 GB; PixelShuffle output
0.80 GB; the `F.interpolate` residual base another 0.80 GB. Sum ≈ 4.5 GB, plus fragmentation and
cuDNN workspace → 5-7 GB. **Scale by 1.52 for the actual 4096x3072 master: ~8-11 GB fp16.**

**So: drop `--tile` at delivery resolution.** Tiling costs seam risk (`--tile_pad` defaults to a
thin 10 px) for no benefit if you are not OOMing. Check `nvidia-smi` rather than tiling
defensively. fp16 is already the default, `half=not args.fp32`.

**The `-dn` knob is the most consequential parameter in your SR stage, and it is almost certainly
sitting at its default.** It is not a scalar denoiser; it is a **weight-space blend**:

```python
if args.model_name == 'realesr-general-x4v3' and args.denoise_strength != 1:
    model_path = [model_path, wdn_model_path]
    dni_weight = [args.denoise_strength, 1 - args.denoise_strength]
```

Deep network interpolation between `realesr-general-x4v3` and `realesr-general-wdn-x4v3`,
evaluated once at load time, free at runtime, both `.pth` files required. `-dn 1` (default) =
pure x4v3 = **weakest denoise, most grain retained**. `-dn 0` = strongest denoise, smoothest,
most plastic.

**[TASTE]** How much grain to keep is *the* creative decision in this project, and right now it is
being made by accident by an undocumented default. The archival-purist position is that grain is
the photographic record and removing it is falsification; the YouTube-audience position is that
heavy grain reads as a bad rip. Both are defensible. What is not defensible is letting a default
choose. **A/B `-dn` at 1.0 / 0.5 / 0.2 on one difficult 200-frame shot and pick with your eyes.**

**[ESTIMATE] throughput ~0.5-0.8 fps at 3840x2160 on a T4.** 32 body convs x 36,864 MAC/px x
8.29 M px ≈ 2.0e13 FLOP/frame; T4 fp16 tensor peak is 65 TFLOPS and small-channel 3x3 convs
realistically hit 15-25% of it → 1.2-2.0 s/frame. **A 12-minute reel at 16 fps is 11,520 frames =
4-6 GPU-hours for the SR stage alone**, roughly a fifth of the weekly budget, before DDColor runs
at all. Worth knowing before adding stages.

### 2.3 Deflicker before colorization, but not the way the argument suggests

The mechanistic case is unusually clean. From fact 6 in §2.1, DDColor's input is *exactly* source
luminance resized to 512x512. Old-film brightness flicker is by definition a frame-to-frame
perturbation of L. So source flicker is not merely correlated with colour instability; it is a
**direct perturbation of the model's only input**, on a network with no temporal state and a known
tendency to flip between colour hypotheses. Every flickering frame presents DDColor with a
different input, so it is entitled to produce a different answer.

DeOldify's own README says the same thing from the other direction: stability comes from giving
the model better input: *"rendering at a higher resolution ... will increase stability of
colorization decisions ... the model has higher fidelity image information to work with and will
have a greater chance of making the 'right' decision consistently."*

**What ffmpeg's `deflicker` actually does**, read from `libavfilter/vf_deflicker.c`:

1. Computes **one scalar per frame**: the whole-frame luma average.
2. Applies a **single global multiplicative gain**: `dst[x] = av_clip_uint8(src[x] * f)` where
   `f = smoothed_window_average / current_frame_average`.
3. **Only touches the luma plane**; chroma planes are copied through unchanged. Safe to run
   before colorization.
4. **It clips.** When `f > 1`, near-white highlights are crushed to 255 and that information is
   gone.
5. **It corrects global exposure pumping only.** No spatial model, so spatially-varying density
   flicker, very common in old film, passes straight through.

**[MEASURED], and this overturns the recommendation I started from.** The research argument was
that `mode=median` beats the default `mode=am` on old film, by analogy with the multimodality
argument in §2.4. That analogy does not transfer, and measurement says so clearly.

Synthetic clips were generated with three flicker profiles, encoded lossless (ffv1), passed
through each mode, and scored by mean frame-to-frame luma delta: which is precisely what the
filter exists to flatten. ffmpeg 8.0.1:

| Profile | `size=7:mode=am` | `size=7:mode=gm` | `size=7:mode=median` |
|---|---|---|---|
| Period-3 square wave, +25% | **85.4%** | **85.2%** | 4.4% |
| Log-normal wobble, σ=0.10 (realistic density flicker) | **86.3%** | **86.4%** | **−0.5%** |
| Isolated blown frames (median's best case) | **86.4%** | **88.2%** | 6.5% |

(reduction in flicker metric; higher is better)

`mode=median` is **useless across all three**, including the isolated-spike profile where outlier
rejection should be its strongest case. Across eight random seeds `am` and `gm` are
indistinguishable (84-88%, `am` ahead 6/8, always by under 1 point). A clipping check on frames
containing a near-255 region showed saturated pixels *decreasing* (70,200 → 14,700 for `am`),
because the correction scales bright frames down.

I did not diagnose why `median` fails inside the filter, and I am not going to assert a mechanism
I did not verify. The empirical result is consistent and reproducible, which is enough to act on.

**Recommendation: use the default `mode=am` (or `gm`, pick either, they tie).** Do not use
`mode=median` here, despite how reasonable the outlier argument sounds.

Free measurement, no pixels modified, `bypass=1` was confirmed bit-identical to the source:

```bash
ffmpeg -i in.mkv -vf deflicker=size=5:bypass=1,metadata=mode=print:file=flicker.txt -f null -
```

writes `lavfi.deflicker.luminance`, `new_luminance` and `relative_change` per frame. Use it to
decide whether your source even has a flicker problem before spending a pass on it.

Then the correction:

```bash
ffmpeg -i in.mkv -vf deflicker=size=7:mode=am -c:v ffv1 -level 3 out.mkv
```

**`vidstab` is a different thing and the names are often paired misleadingly.**
`vidstabdetect`/`vidstabtransform` is *geometric* stabilisation, camera shake and gate weave. It
does nothing for brightness flicker. It is genuinely valuable for old film for a different
reason: **gate weave makes optical flow harder**, so stabilising first improves the flow estimates
that §2.4 depends on. Run it if you have weave; do not expect it to touch flicker.

### 2.4 Motion-compensate the temporal chroma median

**Keep the median; give it geometry.** The statistic is right and the neighbourhood is wrong.

Colorization flicker is **multimodal**, not Gaussian. DDColor does not jitter around a true
colour; it flips between plausible hypotheses, blue shirt on frame 40, brown on frame 41. This
is the central diagnosis of the Deep Video Prior line of work, whose headline contribution
(iteratively reweighted training) exists specifically to *select one mode* rather than average
modes ([arXiv:2010.11838](https://arxiv.org/abs/2010.11838)).

Given multimodality the statistic matters enormously:

- A temporal **mean** of {blue, blue, brown, blue, brown} lands *between* the modes, a
  desaturated grey-brown that is in neither mode and wrong on every frame.
- A temporal **median** *selects* a mode. It is a 50%-breakdown-point estimator: up to half the
  window can be wrong and the output is still exactly right.

**[MEASURED]** on a synthetic multimodal sequence: temporal mean MAE **28.8** vs temporal median
MAE **9.75**. The mean is 3x worse. This also rules out ffmpeg's `tmix`, which is exactly that
unconditional weighted mean.

**The naive version's failure is that it has no notion of motion.** A per-pixel temporal median
asks *"what colour was this pixel location over the last N frames?"* On moving content the right
question is *"what colour was this surface point?"* With a red car crossing a grey road and a
5-frame window, a pixel on the leading edge holds {road, road, car, car, car} → "car", and one
frame earlier {road, road, road, car, car} → "road". The colour boundary **snaps** across the
window instead of tracking the object, and lags at the trailing edge by about half the window.

This is worse here than the generic intuition suggests, for two compounding reasons:

- **The luminance edge stays perfectly sharp** (fact 2). Normally temporal smearing softens the
  whole object, which at least reads as motion blur. Here you get a razor-sharp luminance edge
  with a colour boundary lagging behind it, which the visual system reads as chromatic fringing, considerably more objectionable than the flicker you removed.
- **The chroma is already blocky** (facts 3-4). A temporal median makes those 8x6 px blocks
  *temporally sticky*, so instead of a smooth smear you get visible blocks of colour persisting
  for N frames along the motion path.

**The fix, and why it is nearly free.** Warp neighbours' chroma into the current frame's
coordinate system with optical flow before taking the median. The crucial cost insight: **flow
only needs computing at 512x512**, because that is the chroma's native resolution and it carries
no information above that scale. Warping 4K chroma would be 32x the arithmetic to filter an
upsampled replica of the same 512² samples. At 512x512, `cv2.DISOpticalFlow` is cheap enough to
run **on the CPU concurrently with GPU inference**, i.e. off the critical path entirely.

RAFT is deliberately not used: it needs the VRAM DDColor and Real-ESRGAN are already competing
for, and it is precisely what OOMs All-In-One-Deflicker on a 16 GB card (§4.4).

**Do not skip the occlusion check.** Flow warping without one drags colour across occlusion
boundaries and produces comet-tail artifacts worse than the flicker. The guard is a
forward, backward consistency test: compute flow t→t+1 and t+1→t, warp one by the other, reject
any pixel whose round-trip error exceeds a threshold, and fall back to the unwarped current-frame
chroma there, trusting DDColor rather than inventing a value.

**[MEASURED]**, implemented in `tools/temporal_chroma.py`, verified in
`tools/test_temporal_chroma.py`, 6/6 passing on a fixture with motion, mode flips and an
occlusion:

| Metric | Flow-warped median | Naive median |
|---|---|---|
| Chroma MAE over the moving object | **2.87** | 9.75 |
| Chroma MAE in the edge band | **2.40** | 8.13 |

**3.4x lower chroma error, concentrated exactly where the naive filter fails.** With the centre
frame's own colour flipped to the wrong mode, the filter recovered the true value exactly
(`[-30, 45]` → `[40, -25]`). The test asserts the fixture is actually hard, it fails if the
naive median's edge error is under 5.0, so a passing run cannot be an artifact of an easy
fixture.

Honest sourcing note: I found no published write-up pairing flow-warped chroma-only smoothing
with a modern per-frame colorizer. Both halves are standard and well-published separately (flow
warping is the canonical temporal-consistency construction in video stylization; flow-based
colour propagation is an established line). **This is a recombination of two known primitives,
not a novel method, and I am not going to cite a paper I did not verify.**

**[TASTE]** Window length. 5 frames is a reasonable start; longer windows buy stability and cost
lag wherever flow fails. Judge it on *moving* footage, a static test clip makes every setting
look fine.

### 2.5 Process per shot, not per chunk

Both the deflicker and the temporal median assume the window is one continuous shot. Run them
across a hard cut and the filter fights the cut: deflicker sees a legitimate brightness change as
flicker and pumps against it, and the median blends colours from two unrelated scenes.

Detect shot boundaries first and make them chunk boundaries:

```bash
ffmpeg -i in.mkv -vf "select='gt(scene,0.3)',metadata=mode=print:file=shots.txt" -an -f null -
```

The 0.3 threshold is a starting point, not a constant, check it against the film. For a
single-take 1906 street film this may find nothing at all, which is itself the useful answer.

---

## 3. Making the job fit on free tier

### 3.1 The disk arithmetic

One 4096x3072 8-bit RGB frame is 4096 x 3072 x 3 = **37,748,736 B (36 MB)**. At 16 fps, ten
minutes is 9,600 frames:

| Format | 4096x3072, 10 min | 3840x2160, 10 min |
|---|---|---|
| Raw RGB | **362 GB** | 239 GB |
| PNG **[ESTIMATE** 70-90% of raw on grainy film, DEFLATE cannot predict grain**]** | 254-326 GB | 167-215 GB |
| FFV1 **[ESTIMATE** 35-60%**]** | 127-217 GB | 84-143 GB |

Against `/kaggle/working` at ~20 GB, every one of these is out of the question for a whole film.
**Chunking is not an optimization here, it is the only way the job exists.**

### 3.2 Platform limits

Kaggle's docs pages are client-side-rendered and returned no body text to either WebFetch or
curl, so these come from Kaggle's own product-feedback and discussion threads, first-party, but
forum posts rather than the formal docs. Marked **[community-reported]** accordingly; verify in
the UI.

| Limit | Value | Source |
|---|---|---|
| `/kaggle/working` | ~**20 GB**, persists via Save Version | [community-reported] |
| `/kaggle/temp` | ~**60 GB**, **does not persist** | [community-reported] |
| GPU quota | **30 h/week**, occasionally floating higher | [product-feedback threads] |
| GPUs | single **P100** or dual **T4x2**, 16 GB each, 32 GB system RAM | [product-feedback] |
| Session runtime | **12 h** (raised from 9 h) | [product update] |
| Idle timeout | **60 min** | [community-reported] |
| **Notebook-output file count** | **500 files** | [multiple threads] |
| Dataset size | 200 GB private | [product announcement] |

**The 500-file output limit kills loose PNG frames as a persistence strategy.** 9,600 frames
exceeds it by 19x. Chunk output must be one container file per chunk.

Colab free tier, from the [official FAQ](https://research.google.com/colaboratory/faq.html)
(checked 2026-08-25) **[OFFICIAL]**: still free; **12 h** maximum session; *"The types of GPUs and
TPUs that are available in Colab vary over time"*, **T4 is not guaranteed**; idle timeout
**[NONE] published**; free-tier disk GB **[NONE] published**.

**Do not mount Google Drive for frame folders.** Reading many small files off mounted Drive
collapses from ~300 samples/s to ~1 sample/s
([colabtools #1691](https://github.com/googlecolab/colabtools/issues/1691)), and folders over
10,000 items can make `drive.mount()` fail outright. Use a single archive or container file.

### 3.3 Frame I/O: pipe, do not dump

Real-ESRGAN's own video script already does the right thing, raw frames over stdin/stdout pipes:

```python
ffmpeg.input(video_path).output('pipe:', format='rawvideo', pix_fmt='bgr24').run_async(...)
```

PNG extraction is available behind an opt-in `--extract_frame_first` flag, confirming piping is
the primary path. DDColor by contrast ships **no video pipeline at all**, documented usage is
single-image inference, so the frame harness is yours to build either way.

Ranked for this project:

1. **For resumability: one lossless container per chunk.** A killed session leaves whole playable
   chunk files, unlike a half-written PNG sequence or an in-flight pipe. Use **FFV1 in Matroska**
   if disk binds (better ratio; a Library of Congress *Preferred Format* for video preservation
   since December 2023), or **UT Video / x264 `-qp 0`** if GPU-hours bind. Benchmark data
   [community-reported, Doom9]: UT Video encoded **10.7-13.4x faster** than FFV1 on synthetic
   content at 22-36% larger files; on *natural film* the size gap collapsed to ~14% while FFV1
   was still ~6x slower. **On Kaggle, GPU time is the scarce quota-limited resource, trade disk
   for CPU-side encode speed.**
2. **Within a session: raw pipes**, no intermediate file I/O at all.
3. **Avoid loose PNG dumps** on Kaggle regardless of speed, the 500-file limit is a dead end.

### 3.4 Chunked resumable structure

Segment and rejoin are both documented ffmpeg muxers; the composition is not a named pattern.

```bash
# split into independently-decodable chunks with a manifest
ffmpeg -i in.mkv -f segment -segment_time 30 -reset_timestamps 1 \
       -segment_list chunks.txt -c:v ffv1 -level 3 chunk%04d.mkv

# ... process one chunk per session ...

ffmpeg -f concat -safe 0 -i done.txt -c copy final.mkv
```

`individual_header_trailer` defaults true, making each segment independently usable, which is
exactly what resumability needs. `concat` needs `-safe 0` for absolute paths.

**For the 85 GB master, segment straight from the URL** rather than downloading. ffmpeg reads
HTTP input with range requests, so `-ss`/`-t` against the archive.org URL fetches only the bytes
for that chunk. That, not local storage, is what makes this film processable on free tier.

---

## 4. Bigger changes

### 4.1 Defect repair is the actual hole, use RRTN, not the famous one

**Nothing in DDColor + Real-ESRGAN removes a scratch, a blotch, or a splice mark**, which is what
a viewer of a 1906 print actually sees.

**Two corrections to the plan's premise.** There is no `microsoft/Bringing-Old-Films-Back-to-Life`, that URL is a 404. The repo is [`raywzy/Bringing-Old-Films-Back-to-Life`](https://github.com/raywzy/Bringing-Old-Films-Back-to-Life),
and the venue is **CVPR 2022, not ICCV 2021**. It is also unusable: it pins **Python 3.6.13 /
torch 1.5.1 / CUDA 10.1** (Python 3.6 reached EOL in December 2021), last commit **2023-02-04**,
with **30 of 32 issues open including crashes on the authors' own test data**. Install pain 5/5.

Its successor is the better bet on every axis:
**[RRTN](https://github.com/mountln/RRTN-old-film-restoration)**, *Restoring Degraded Old Films
with Recursive Recurrent Transformer Networks* (2024). Eleven stars, so almost nobody has heard
of it.

| | Bringing Old Films (raywzy) | RRTN (mountln) |
|---|---|---|
| License | none | **MIT** |
| Weights | author's SharePoint | **GitHub Releases** (41 MB x2 + RAFT 20 MB) |
| Python / torch | 3.6 / 1.5.1 / cu101 | **3.10 / 2.1.2 / cu118** |
| Last commit | 2023-02-04 | **2024-09-01** |
| OOM guidance | none | **explicit, in the README** |

**VRAM: [NONE] as a number, [QUOTED] as a procedure**, *"If GPU memory is not enough, try
reducing `temporal_length`, `temporal_stride`, or the resolution of the input video."*

**The hidden cost is the recursion.** Defaults are `--temporal_length 15`, `--temporal_stride 3`,
`--recursion_threshold 43` (PSNR), `--max_recursion 4`. The model re-runs itself up to four times
until input/output PSNR exceeds 43 dB. **Budget up to 4x the single-pass cost.** That is the "R"
in RRTN and the README does not flag it.

**Install pain 4/5, two specific landmines.** `mim install mmcv`, mmcv's own badges read
`pytorch 1.8~2.0` and `cuda 10.1~11.8`, and Kaggle is well past both; with no matching wheel pip
silently falls back to the source tarball and compiles every CUDA op (20-40 minutes when it
works, and it often does not). Then `deform_attn_cuda_kernel.cu` is JIT-compiled at import via
`torch.utils.cpp_extension.load`, needing `ninja` and a matching `CUDA_HOME`. Neither is fatal on
Kaggle, nvcc is in the image, but **budget a whole session for the install and cache the
compiled artifacts to a Kaggle dataset once it works.**

**Run it before colorization**, at ≤1080p, then transfer the repaired luma back to full
resolution. Doing it first matters twice over: DDColor is currently predicting chroma partly from
dust, and your temporal chroma median is partly fighting artifacts a defect-repair stage removes
at source.

### 4.2 Exemplar propagation, if you want the cure rather than the symptom

The planned ColorMNet upgrade is sound in direction. The current best-matched option is
**[CMNET2](https://github.com/dan64/cmnet2)** (2026), ColorMNet plus XMem++ three-tier permanent
memory, built for long videos.

It is the **only 2026-maintained colorizer in this survey with a measured 16 GB performance
table**, README **[QUOTED]**: *"16 GB VRAM recommended for long videos"*, measured on an RTX
5070 Ti 16 GB at **2.63 fps** at 960x730 with 158 references, or **5.80 fps** with `--max_side
512` plus chroma transfer. Its sliding permanent-memory window is explicitly VRAM-bounded, which
is exactly the property a fixed 16 GB budget needs.

ColorMNet itself publishes only *relative* memory (paper §4: *"only 11.2% of DeepRemaster, 10.0%
of DeepExemplar and 5.4% of BiSTNet"* at 960x536 on an A6000), no absolute GB figure.

**Install pain 4/5 for both**, and it is the same hurdle the plan already anticipated:
`spatial_correlation_sampler` from Pytorch-Correlation-extension needs `python setup.py install`
against CUDA. Same class of pain as RRTN's `deform_attn`.

Note the `--max_side` + chroma-transfer trick, because it is the right pattern for RRTN too: **run
the expensive model small, transfer the result back at full resolution.** That is the same
insight as §2.4, chroma does not need to be computed at delivery resolution.

**Avoid the HAVC / vs-havc stack on Kaggle** despite it integrating everything: VapourSynth plus
Hybrid Windows Addons, prebuilt correlation-sampler wheels only for CUDA 12.4/py3.12, and it
assumes a Windows desktop. Mine its *filter designs* (chroma resize + transfer, B&W tune,
ChromaBound) rather than adopting the stack.

### 4.3 What to avoid despite the hype

Every rejection below has a receipt.

| Model | Why not |
|---|---|
| **BasicVSR++** | **[QUOTED]** issue #10, RTX 3080 **16 GB**, `--max-seq-len=1`: *"GPU OOM occurred when I use input image size is 1920 X 1080."* Also x4-only, and the EOL mmedit toolchain |
| **RealBasicVSR** | Maintainer on 8 GB: *"I guess 8G is not enough"*; loads the whole video into host RAM (a 52 GB attempt was reported) |
| **STAR** (ICCV 2025) | **[QUOTED]** *"4x, 72 frames, width 426, height 240, requires around **39GB**"* |
| **SeedVR2** (ICLR 2026) | **[QUOTED]** *"1 H100-80G can handle videos with 100x720x1280"*. Also wants bf16; **Turing has none** |
| **Vivid-VR** (ICLR 2026) | **[QUOTED]** *"For a 121-frame video, it requires approximately **43GB**"* |
| **FaithDiff** (CVPR 2025) | Markets itself as "Classic Film Rejuvenation", but it is per-image → flicker, and its 5 GB path needs **fp8**, which requires sm_89+. **T4 is sm_75** |
| **Real-CUGAN, APISR** | Fit easily, wrong prior, trained on anime |
| **DAT** | Trained on DIV2K/Flickr2K **bicubic** degradation; fails on real film degradation |

The pattern: **the 2025-26 frontier moved to multi-GPU 80 GB diffusion and left the 16 GB tier
roughly where it was in 2022.** That is not a reason for despair, it means the cheap, well-worn
tools are still the right ones, and the real wins are in §2 rather than in a newer model.

### 4.4 Blind temporal consistency, assessed and mostly rejected

| Method | Per-clip training? | Verdict |
|---|---|---|
| Bonneel et al. 2015 | No | C, no license, last push 2016-06-28. Build risk |
| **Lai et al. ECCV 2018** | **No, pretrained, feed-forward** | **The one worth trying** |
| Deep Video Prior (NeurIPS 2020) | **Yes, every clip** | Dead toolchain. Skip |
| All-In-One-Deflicker (CVPR 2023) | **Yes, neural atlas** | OOMs below 4K. Proxy only |

**Lai et al.** is the realistic learned option: *"At test time, our model does not require
computing optical flow and thus achieves real-time speed even for high-resolution videos"*, and
colorization is in its training distribution. The README's scary FlowNet2 dependency applies to
*training*, in `test_pretrained.py` the `Resample2d` import is on line 13 and **never called
anywhere in the file**. Stub that one import and inference is pure PyTorch. That is the single
biggest practical unblocker for this repo, and it is why people give up on it unnecessarily.
The model is also **residual** (`frame_o2 = frame_p2 + output`), which is a good sign for
sharpness, it nudges rather than resynthesises.

**DVP: skip.** Pins `python=3.6.12`, `tensorflow=1.13.1`, `cudatoolkit=10.0.130`, not installable
on a current Kaggle image. And the cost is prohibitive anyway: **[QUOTED]** *"For an 800x480 video
frame, our approach costs 80 ms for each iteration"*, with colorization needing 50 epochs →
~4 s/frame at 800x480 **on an RTX 2080 Ti**. A 9,600-frame film is ~10.7 GPU-hours at a
resolution you do not want. A community PyTorch port exists but is unverified.

**All-In-One-Deflicker** is on-topic twice (it targets *"old movies and cartoons"* and
colorization flicker) and its filter stage claims *"about 3000M GPU memory"*, but issue #5 shows
it OOMing at 1280x1280 on a **14.75 GiB** card, i.e. exactly this GPU, with the author confirming
*"The preprocess optical flow cannot support high-resolution images."* Usable only on a
downscaled proxy, which, given chroma is 512x512 anyway, is less absurd than it sounds. Worth one
pilot, not a plan.

**ffmpeg-native fallbacks**, ranked: `atadenoise` is the best zero-code option because its A/B
thresholds give it a crude motion guard, the only one of these with any. **Chroma-only is
`p=6`, not `p=2`**: the option is `AV_OPT_TYPE_FLAGS`, a bitmask with default 7; `p=2` would
filter only the U plane and produce a colour cast. `hqdn3d` has no motion adaptation at all and
its parameter order is misleading (`hqdn3d=0:0:4:3` sets *luma* temporal, not chroma). `tmix` is
the temporal mean, ruled out by the measurement in §2.4. `minterpolate` is frame-rate
conversion, not a smoother.

---

## 5. Throughput

### 5.1 The bf16 trap, which is worse than it looks

**T4 is Turing (sm_75): fp16 tensor cores, no bf16, no TF32.** The datasheet gives 320 Turing
Tensor Cores, **65 TFLOPS fp16 against 8.1 TFLOPS fp32**: an 8x peak ratio. TF32 is
Ampere-and-later (PyTorch docs). Usable VRAM is **15,360 MiB**, not 16 GB: GDDR6 ECC reserves
6.25%.

The trap: **`torch.cuda.is_bf16_supported()` returns True on a T4.** It defaults to
`including_emulation=True`, and on sm_7x it merely checks that a bf16 tensor can be *allocated*, 
which succeeds. PyTorch PR #129699 states Volta and Turing "do not have HW support for
bfloat16 ... it is emulated in software." So `torch.autocast('cuda', torch.bfloat16)` on a T4
**runs silently, emulated, and slower than fp32.** No error, no warning, just a slow job.

Check the compute capability directly, not the helper:

```python
major, minor = torch.cuda.get_device_capability()
dtype = torch.float16          # correct for both T4 (7.5) and P100 (6.0)
use_bf16 = major >= 8          # Ampere and later only
```

**Kaggle also offers P100 (Pascal, sm_60), and it behaves differently.** No tensor cores at all;
fp16 is only "up to twice FP32 throughput" via packed CUDA-core operations (GP100 whitepaper), so
expect **~1.2-1.5x from fp16 there, not 3x**. HBM2 ECC is free, so a P100 exposes about 1 GiB
more usable memory than a T4. **Do not assume which card you got, print it.**

### 5.2 fp16 safety, per model

**Real-ESRGAN `realesr-general-x4v3`: safe, and verified rather than assumed.** Peak activation
was probed at **160** even on an adversarial checkerboard input, against fp16's max of 65,504,
that is **x409 headroom**, and fp32 vs fp16 output differs by **0.004**. NaN and black-output
reports do exist (issues #632, #330, #479) but every one is RRDBNet on GTX 16-series or old
cards, or ncnn/ROCm. None is a T4, and none is x4v3.

**One gotcha that will silently cost you half your throughput:** the CLI defaults to fp16
(`half=not args.fp32`), but **`RealESRGANer.__init__` defaults `half=False`**. Calling the
library directly, which a notebook does, gives you fp32 unless you pass `half=True`.

**DDColor has no fp16 support anywhere**, verified across `pipeline.py`, `infer.py`, `model.py`,
the demos, and all 68 issues. Wrapping the forward pass in `torch.autocast('cuda',
dtype=torch.float16)` is the sound approach rather than `.half()`, because autocast keeps
ConvNeXt's LayerNorms and the transformer softmaxes in fp32 while sending convolutions and matmuls
to fp16. `.half()` would cast everything and is materially riskier. **This is untested, validate
on ~20 frames against an fp32 reference before committing a whole film to it.**

### 5.3 Batching

**DDColor is 227.98 M parameters**, cross-checked two ways: the checkpoint is 911,914,869 bytes ÷
4 = 227.98 M, matching the paper's stated 227.9 M. It is not a small model, and it is hard-coded
to batch 1 by the `.unsqueeze(0)` in `process()`.

**[ESTIMATE]** ~1.66 GiB at batch 1, ~834 MiB per additional frame → **batch 8 at fp32, 16 at
fp16**. Caveat: issue #57 reports an unexplained 8 GB, so ramp the batch size and watch
`nvidia-smi` rather than jumping straight to 16.

Adding the batch path is the same ~15-line change that gives you the raw 512x512 `ab` tensor
needed for §2.4, so the two changes are one piece of work.

Real-ESRGAN x4v3 at 1,213,296 parameters (verified by a `strict=True` load) is never the
constraint.

### 5.4 Settings for grainy 4K film

- **`--tile 512 --tile_pad 16`**, about 2.42 GiB. Note this partly contradicts §2.2's "drop
  `--tile`": tiling cannot shrink the 2.24 GiB output buffers, and `pre_process` still uploads
  the whole frame anyway (issue #184). At 4096x3072 rather than UHD, tiling is the safer default.
  Try untiled first and watch `nvidia-smi`; fall back to `--tile 512` if it OOMs.
- **`cudnn.benchmark = True`** is appropriate here, because it only helps when input shapes are
  constant, and DDColor's always are, at 512x512.
- **`--outscale 2` saves nothing.** The network always runs 4x, then a CPU Lanczos resize
  discards the result.

### 5.5 ⚠ A conflict worth knowing about: which way `-dn` runs

Two independent investigations disagreed on the **direction** of the `-dn` knob, and neither is
obviously wrong.

| Source | Claim |
|---|---|
| Repo README + DNI code reading | `-dn 1` = pure `x4v3` = **weakest** denoise, keeps the most grain |
| Empirical probe of both checkpoints on synthetic grain | `x4v3` **strips** grain (HF std 0.009 vs clean 0.010); `wdn` **keeps** it (0.017), so `-dn 1` denoises *more* |

The mechanics are agreed: `dni_weight = [denoise_strength, 1 - denoise_strength]` blends `x4v3`
with `realesr-general-wdn-x4v3`, so `-dn 1` is pure x4v3 and `-dn 0` is pure wdn. What is disputed
is which checkpoint actually smooths more, and the naming ("wdn" = with-denoise, versus a flag
called denoise_strength that *decreases* wdn's weight) is genuinely confusing enough that the
README's own wording may be misleading.

**This does not change what to do, which is reassuring.** Both analyses end at the same
instruction: **render one difficult 200-frame shot at `-dn` 1.0, 0.5, 0.2 and 0.0, and choose with
your eyes.** Just do not assume you know which end is which, label the renders by the flag value,
not by "more grain" or "less grain", and let the footage tell you.

### 5.6 Frame decimation, 33% of GPU time, free **[MEASURED]**

This one is specific to the first film and it is the largest single throughput win available.

The archive.org item description says so in the uploader's own words:

> "It was scanned by Adrianne Finelli at Internet Archive on October 11, 2018 using IA's
> Lasergraphics ScanStation scanner set at 5120 x 3840 and **scanned to 16fps**. Due to Final Cut
> Pro limitations, this full-width scan is 4096 x 3072 pixels wide and **every second frame is
> doubled for a playback framerate of 24fps**."

Two source frames become three, so **one frame in three is an exact duplicate.** Colorizing and
super-resolving those is pure waste.

Verified by building a clip to that exact 16→24 recipe and decimating it:

```
unique source frames 32 -> emitted 48 (ratio 1.50)
after mpdecimate        : 32
duplicates removed      : 16  (33.3% of GPU work saved)
matches unique source   : True
```

Implemented as `pipeline.decimate_chunk()`, covered by
`test_pipeline.py::test_decimate_removes_doubled_frames`, which also asserts that decimating an
already-unique clip drops nothing.

**Run it before the GPU stages and re-time afterwards.** One subtlety cost a debugging round and
is worth repeating: do the re-timing with `setpts` *inside* the filter graph plus
`-fps_mode passthrough`. Combining `mpdecimate` with an output `-r` makes ffmpeg re-time a stream
whose timestamps still describe the old rate, and it drops real frames, 48 doubled frames came
out as 23 instead of the correct 32.

Two further consequences of that description:

- **The film's native rate is 16 fps**, so the "24→16 fps speed correction" in the current plan is
  describing the wrong operation. What the master needs is duplicate removal, then presentation at
  16 fps (or whatever rate you choose, see §6).
- **The 4K master is a full-width scan including the sprocket perforations**, which is why a crop
  step exists. The 2048x1536 `MarketStreet_4K_to_2K_cropped_higher_contrast.mov` is **already
  cropped and contrast-adjusted**, and the uploader recommends it: *"For most people, your best
  download link will be to the 'cropped higher-contrast version.'"* At 2.5 GB it also fits on free
  tier trivially. **Use it for the first pass**; graduate to the 4K master by URL segmentation once
  the pipeline is proven.

Also worth recording in the provenance file, the uploader grants reuse explicitly:

> "Anyone may reproduce or reuse this scan. Please attribute it to its source: Prelinger Archives."

---

## 6. Frame interpolation, a taste decision, not a recommendation

**[TASTE]** [Practical-RIFE](https://github.com/hzwer/Practical-RIFE), MIT, pure PyTorch, no
compiled ops, `python <= 3.11`. Models frozen since 2024-09-21. **[NONE]** published VRAM; the
README says *"If your video has high resolution, such as 4K, we recommend set --scale=0.5"*.

The aesthetic question is genuinely unsettled and I am not going to prescribe an answer. Silent
films were hand-cranked at **12-40 fps**, so there is no single "true" rate to restore, the
premise of "restoring correct motion" is shakier than it sounds. Interpolating to 60 fps produces
the soap-opera effect that many viewers dislike on any material and that archivists object to on
principle; interpolating to ~24 fps is where most restorations land and is far less contentious.

**A/B guidance:** take one 200-frame shot with fast lateral motion (a passing streetcar is ideal),
render it at native, at 24 fps and at 60 fps, and watch all three at full size a day apart. Judge
the fast-motion shot, not a static one, interpolation artifacts hide in slow footage. **If you
do interpolate, do it last**, after defect repair: RIFE will happily interpolate a scratch into a
smooth moving object.

---

## 7. Do these before rendering anything

### 7.1 Verify the frame rate

The 404.5 s / 716.8 s discrepancy in §0.1 means the "24→16 fps" assumption is unverified and
probably wrong for this scan.

```bash
ffprobe -v error -select_streams v:0 \
  -show_entries stream=r_frame_rate,avg_frame_rate,nb_frames,width,height,duration \
  -of default=noprint_wrappers=1 \
  "https://archive.org/download/MarketStreet19064KScan20181016/Market_Street_1906_4K_scan_20181016.mp4"
```

Run it against the master and the 2K version and compare. A wrong rate is baked into every frame.

### 7.2 Measure the flicker before treating it

Run the `bypass=1` metadata probe from §2.3. If `relative_change` is already small, skip the
deflicker pass and spend the time on defect repair instead.

### 7.3 Decide the grain question deliberately

A/B `-dn` at 1.0 / 0.5 / 0.2 (§2.2). This is the film's whole look and it is currently a default.

---

## 8. What is still uncertain

1. **Real-ESRGAN VRAM at 4096x3072**, estimated, not measured. The repo publishes nothing. Check
   `nvidia-smi` on the first run rather than trusting the arithmetic in §2.2.
2. **Whether RRTN installs cleanly on a 2026 Kaggle image.** The mmcv wheel gap and the JIT
   `deform_attn` build are both real risks. Budget a throwaway session to find out.
3. **Why ffmpeg's `deflicker mode=median` fails.** The measurement is reproducible; the mechanism
   was not diagnosed.
4. **Whether DDColor tolerates `input_size` above 512 at inference.** `build_ddcolor_model`
   accepts and forwards the parameter, so it is at least plumbed through, and DeOldify's
   `render_factor` observation suggests raising it is a *stability* lever, not just a quality one.
   Untested here, and worth an A/B because it attacks flicker at the cause.
5. **Kaggle's current limits**, the docs pages are JS-rendered and could not be read directly.
   Every figure in §3.2 is community-reported; confirm in the UI before planning around them.
