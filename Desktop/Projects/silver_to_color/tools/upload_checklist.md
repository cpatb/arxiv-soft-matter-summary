# Upload checklist, master copy

`make_metadata.py` writes a pre-filled copy of this per film. This is the version to edit when
the process changes.

---

## Before you render

- [ ] **Publication year confirmed from a second source.** archive.org's `year` is
      uploader-supplied and unvalidated, a Prelinger short was observed dated 2026. Check the
      Catalog of Copyright Entries or AFI.
- [ ] **Published 1931-1963?** Check renewal, per film. Two corrections worth having:
      the band starts at the PD cutoff (1931 in calendar 2026) and **moves forward every
      year**, and renewals live in the **annual** Catalog of Copyright Entries, Parts 12-13.
      The **cumulative** volumes contain no renewal entries at all, which is the usual reason
      a search comes back empty and reassuring.
- [ ] **Do not assume non-renewal.** Films had the **highest renewal rate of any class**,
      roughly 74% and 65% in the studied samples, against about 7% for books. "Probably never
      renewed" is a reasonable prior for a book and a bad one for a film.
- [ ] **Checked `archive.org/details/movie_research`.** Archive.org staff maintain
      OKAY / NOT-OKAY / Unknown lists with renewal registration numbers. Films on the NOT-OKAY
      list are live on the site today under Creative Commons Public Domain Marks, two of them
      uploaded in 2025. Archive.org states plainly that it "does not make guarantees as to the
      copyright status of items" and that it has no copyright experts.
- [ ] **Item is not dark.** A withdrawn item usually means a claim succeeded against it. The
      tools here raise on `is_dark` rather than treating it as a missing item.
- [ ] **Foreign film?** URAA (17 U.S.C. 104A) restored US copyright in many foreign works that
      had fallen into the public domain. A 1923 French feature is not automatically safe.
- [ ] **Source is the rawest available scan**, not someone else's restoration. Restorers assert
      rights in restoration work; whether that holds is contested, and you do not want to be the
      test case.
- [ ] **Any embedded audio muted.** A modern score on a restored silent film is a live copyright
      in both composition and recording, and that claim is *valid* even though the film is PD.
- [ ] **Score is self-composed or provably clear.** Not just for copyright: videos with music
      claims are **ineligible for channel memberships**, and "a significant number of ineligible
      videos" can disqualify the whole channel.
- [ ] **Frame rate verified with `ffprobe`, not assumed.** It is baked into every frame.
- [ ] **Checked for padded/duplicated frames.** If the scan pads to a modern rate, decimation is
      free GPU time, one frame in three on the Market Street master.
- [ ] **Shot boundaries detected.** Deflicker and the temporal median both misbehave across cuts.

## Render

- [ ] Grain level (`-dn`) chosen **deliberately**, not left at a default. This is the film's whole
      look. Label the A/B renders by flag value, since the two available sources disagree on which
      direction the knob runs.
- [ ] Colorization reviewed for the colours that are known-wrong: skin, sky, signage, uniforms.
- [ ] Before/after rendered for the first 20 seconds.
- [ ] Delivered at 2160p, for the bitrate allocation (35-45 Mbps vs 8 Mbps at 1080p), not
      because there are 8 million pixels of detail.

## Metadata

- [ ] Title under 100 characters.
- [ ] Description pasted from `<slug>.description.txt`.
- [ ] Tags pasted, under 500 characters total.
- [ ] Chapters added if over 8 minutes.
- [ ] Mid-roll ads enabled if over 8 minutes, the second-largest RPM lever after narration.
- [ ] End screen and cards set.

## The two that protect the channel

Everything above costs you a video if you skip it. These two cost you the channel.

- [ ] **"Altered or synthetic content" disclosed in Studio.**
      Upscaling is explicitly exempt from YouTube's disclosure rules; **colorization is not
      addressed anywhere**, and it sits between an exempt bullet ("color adjustment") and a
      triggering one ("altering footage of real events or places"). Disclosure is stated to carry
      **no monetization or reach cost**; persistent non-disclosure risks YPP suspension. The
      asymmetry decides it.

- [ ] **Narration present, and it is the spine of the video, not a caption over footage.**
      This is the single highest-value item on the list. It appears in five of six of the reused
      content policy's allowed-side bullets; it is the fix for the generic/repetitive-content
      risk; it satisfies TikTok's 60-second minimum and Meta's originality carve-out; and it is
      worth **3-5x on RPM** ($7-9 narrated vs $1.50-4 unnarrated, because silent footage is
      language-neutral and attracts low-RPM global traffic).

      The most instructive rejection wording found: *"Appeal didn't demonstrate enough value-add
      editing. The editing shown in your appeal didn't substantially change the original work."*
      The appellant had submitted an editing timelapse and project files. **Proving labour does
      not work. On-screen and audible difference does.**

## After publishing

- [ ] `<slug>.provenance.json` filed somewhere findable in a hurry.
- [ ] Checked for Content ID claims within 24 h. If one lands: "public domain" is **not** an
      enumerated dispute reason, file under *"misidentified or an error was made"* and cite the
      archive item URL, the film's year, and the fact that PD footage is ineligible as a Content
      ID reference (`support.google.com/youtube/answer/2605065`). Note that disputing is not
      risk-free: a claimant can escalate to a takedown, and a copyright strike is disqualifying
      for YPP.
- [ ] Short-form cut is **60 s or longer** and carries the narration, not just the reveal.
- [ ] Tracking row added: archive ID, narrated Y/N, length, mid-roll count, RPM, audience
      geography, claims, GPU-hours.

## Periodically

- [ ] **Applied to YPP before 2027-02-01?** Entry thresholds double for new applicants on that
      date (4,000 → 8,000 watch hours) and existing members are grandfathered. With a ~1-month
      review queue, the real deadline is around the end of December 2026.
- [ ] Advanced-features verification enabled, it is a separate YPP requirement, easy to miss,
      and not instant.
