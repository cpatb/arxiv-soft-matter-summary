# Silver to Color: production slate, episodes 2-11

Compiled 2026-09-30 / 10-01. Episode 1 is *A Trip Down Market Street* (1906) and is not repeated here.

**How to read this.** Every file size, resolution, duration and codec below was read in this
session from the archive.org metadata API (`archive.org/metadata/<id>`) or by `ffprobe` against the
remote file (headers only, nothing downloaded). bpps = size / (width x height x duration_s).
Anything not checked that way is tagged **[UNVERIFIED]**. **No reel was watched.** Content
descriptions come from the holding archive's catalogue text, so the first job on each film is to
watch the 640x480 derivative before committing GPU time.

**Rights in one paragraph.** US public domain as of 2026 covers works *published* in 1930 or
earlier. Three things on this slate need more than a date:
- *Federal works* (US Army Signal Corps, Navy) are PD under 17 USC 105 whatever their date.
- *Actualities with no publication evidence* (local parade footage, home movies, outtakes) are a
  separate case. An unpublished work is not governed by the publication rule. For an anonymous
  pre-1978 unpublished work the term is 120 years from creation (1915 + 120 = 2035), and for a known
  author it is life + 70. This is my reading of the statute **[UNVERIFIED, get a lawyer's view]**.
  The Content ID risk on such footage is small; the legal position is still not the same as Market
  Street's.
- *A Public Domain Mark or "no restrictions" statement* on an item is the uploader's or holding
  institution's assertion, not a determination.

No item below is `is_dark`: every metadata call returned a live `metadata` block.
**Library of Congress could not be checked.** loc.gov returned HTTP 403 (a Cloudflare challenge)
to curl and to WebFetch, and the Chrome extension was not connected. NARA's catalog proxy API did
work and was used.

### A tooling finding that changes how to read bpps
`scout_films.py` flags the Prelinger 2K masters `re-encode-not-master(bpps=0.81)`. It flags the
California Revealed HD files `LIKELY-UPSCALE-NOT-SCAN(bpps=0.45)`. Both flags misfire:
- The 2K Prelinger `.mov` files are **HEVC Rext at about 20.5 Mb/s**, a mezzanine codec (probed:
  `codec_name=hevc|profile=Rext`). Their 0.81 comes from codec efficiency, not from a lossy
  re-encode of something better.
- The California Revealed files are institutional **H.264 access copies at about 10 Mb/s**. Low
  bpps there means compression, not synthetic resolution. Whether they are upscaled is
  **[UNVERIFIED]** until someone looks at the grain.
- The thresholds in `workflow.md` (≈10 for a genuine scan) were calibrated on **ProRes 4444**, and
  bpps also scales with frame rate (16 vs 24 fps).

Compare bpps only within a codec. The scout script's flags should key on codec, not just the
number. This is not fixed here.

---

## Ranked slate

| # | Film | Year | Best usable file | Res | Size | Dur | bpps | Story | Scan | Rights |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | PPIE Opening Parade, San Francisco | 1915 | HEVC full-res | 4096x3072 | 6.82 GB | 8.2 min | 1.11 (HEVC) / 10.19 (ProRes) | A | A+ | B (unpublished?) |
| 2 | Start of Epic Flight (PN-9 to Hawaii) | 1925 | 2 x HEVC Rext .mov | 2048x1536 | 0.48 + 0.42 GB | 5.8 min | 0.81 | A+ | B+ | B+ |
| 3 | America Goes Over (AEF in WWI) | 1918 | NARA 111-M-600E reels | 1920x1080 | 3.7-5.0 GB/reel | 55 min | 3.04 | A | B+ | A |
| 4 | Rambles 'Round Rio (Kineto Review) | 1921 | ProRes 4444 master | 4096x3072 | 112.5 GB | 14.0 min | 10.62 | B | A+ | A |
| 5 | U.S. Navy of 1915 (Lyman H. Howe) | 1915 | H.264 10 Mb/s | 1438x1080 | 0.94 GB | 12.2 min | 0.83 | A- | B | A- |
| 6 | Official Opening of the LA Aqueduct | 1913 event | H.264 access HD | 1920x1440 | 1.54 GB | 20.5 min | 0.45 | A+ | B- | B- |
| 7 | The Mutual Weekly No. 101 (newsreel) | 1916 | H.264 access HD | 1920x1440 | 1.24 GB | 16.8 min | 0.45 | B+ | B- | A- |
| 8 | Allied Expeditionary Force in Russia | 1918-19 footage | HEVC Rext .mov | 2048x1536 | 7.36 GB | 47.9 min | 0.81 | A | B+ | B |
| 9 | Railroads in United States History | 1920 | ProRes 4444 master | 4096x3072 | 123.2 GB | 11.2 min | 14.56 | B | A+ | A |
| 10 | Ford Film Collection selects | 1910s-1930s | HEVC Rext .mov | 2048x1536 | 6.27 GB | 40.8 min | 0.81 | B+ | B+ | B+ |

Ranking = story x scan x rights, with rights weighted heaviest, because a strike costs the channel
and a weak scan costs one video.

---

### 1. Panama-Pacific International Exposition Opening Parade (San Francisco, 1915)

- **Year + evidence:** the event is the PPIE opening, 20 February 1915 (Wikipedia / NPS pages
  surfaced by search). The archive item has **no date field**. Its description reads "Dupe of
  original nitrate film of Panama-Pacific International Exposition (PPIE) opening parade with
  additional scenes around downtown San Francisco." The footage cannot predate the event. Whether
  and when it was *published* is unknown: no maker, no copyright notice, no newsreel attribution.
- **Source:** `200333_Panama-Pacific_International_Exposition_Opening_Parade`.
  https://archive.org/details/200333_Panama-Pacific_International_Exposition_Opening_Parade
  (Prelinger, FFDW-sponsored scan, added 2023-05-18).
- **File to use:** `...Opening_Parade.intros.hevc.mp4`: HEVC, **4096x3072 @16 fps, 489.3 s,
  6.82 GB, bpps 1.11**. It fits `/kaggle/working`.
  Master: `...intros.mov`, ProRes 4444 12-bit, 4096x3072, 479.4 s, 61.46 GB, bpps 10.19
  (segment from URL if needed).
  The HEVC runs 9.9 s longer than the master. Check for padded intro frames before rendering.
- **Rights checked:** no `licenseurl`, `rights` or `possible-copyright-status`. Not dark.
  Prelinger collection. The risk is the unpublished-actuality question above. An AI colorization
  of PPIE footage already exists (SFGate on "Reimagined Reels"), which suggests no one is actively
  claiming PPIE footage **[UNVERIFIED]**.
- **Narration angle:** nine years after Episode 1's street burned, San Francisco throws a world's
  fair to prove it rose. The parade route runs through the rebuilt downtown. The episode can
  intercut a 1906 frame of Market Street with the same corner in 1915 and tell the
  rebuilt-in-record-time story, the Panama Canal's opening that the fair celebrated, and the Tower
  of Jewels.
- **Shorts hooks:**
  - "This street was rubble in 1906. Nine years later, 18,000 people marched down it."
  - "San Francisco built a fake city out of plaster just to prove the real one was back."
- **Verdict:** best scan on the slate and the natural Episode 2, since it is Market Street's
  sequel. The risk is moderate and legal rather than practical: publication is unproven.
  **Go**, after a viewing pass confirms what the parade footage actually shows.
- The parade crowd counts and the Tower of Jewels detail need sourcing **[UNVERIFIED]**.

### 2. Start of Epic Flight: the PN-9 No. 1 flight to Hawaii (1925)

- **Year + evidence:** the flight departed 31 August 1925 (Hawaii DOT Aviation history; USNI
  *Naval History* Feb 2019; Wikipedia, all via search). The archive record gives "circa 1925",
  creator "Pathe News", can description "Newsreel footage depicting the embarkation of 'the longest
  nonstop flight ever undertaken'". A topical newsreel released in 1925 is PD (published before
  1931). The release itself is **[UNVERIFIED]**: no issue number was found.
- **Source:** `200677_Start_Of_Epic_Flight`. https://archive.org/details/200677_Start_Of_Epic_Flight
  (Prelinger / FFDW, added 2024-05-23).
- **Files to use:** `200677_Start_of_Epic_Flight_001_master.intros.mov` (HEVC Rext, 2048x1536
  @16 fps, 186.9 s, 0.48 GB, bpps 0.81) and `..._002_master.intros.mov` (162.5 s, 0.42 GB,
  bpps 0.81). The description notes the content is split across two files "due to scanning issue".
- **Rights checked:** `rights: Military`. No license field. Not dark.
  CriticalPast sells the same takeoff footage as "HD stock video" and has a YouTube upload of it
  (search result). Pathé-branded material can draw claims from rights-holder catalogues. Expect a
  possible Content ID claim and keep the provenance file ready **[UNVERIFIED]**.
- **Narration angle:** the US Navy tries the first nonstop flight to Hawaii. PN-9 No. 1 runs out
  of fuel short of Oahu and vanishes. The crew is given up for dead. They cut sails from the wing
  fabric and sail the flying boat for nine days, until a submarine finds them off Kauai. The film
  is only the departure; the episode is the story of what the audience knows comes next.
- **Shorts hooks:**
  - "These five men took off for Hawaii. Nine days later they arrived by sailboat."
  - "The Navy announced them lost at sea. They were busy turning their airplane into a boat."
- **Verdict:** strongest story per minute on the slate. It is short (5.8 min), so it is ideal for
  Shorts and cheap to render, but it needs stills or maps to fill a long-form episode. **Go.**

### 3. America Goes Over (US Army Signal Corps, 1918)

- **Year + evidence:**
  - Wikipedia (via search): a 1918 documentary by the Army Signal Corps under the Committee on
    Public Information.
  - NARA catalogues it as `(AMERICA GOES OVER.)`, naId 35697, series 111-M, Records of the Office
    of the Chief Signal Officer. NARA gives a production date of 1936, which is the catalogue
    accession, not the film.
  - Prelinger's prints carry can notes "1927 print".
  - A federal work is PD under 17 USC 105 regardless of date.
- **Sources:**
  - NARA https://catalog.archives.gov/id/35697, files `111-m-600e-r1..r5.mp4`.
  - Alternative: Prelinger `202309_America_Goes_Over_P1-2` (4.70 GB, 30.6 min) and
    `202204_America_Goes_Over_P3-5` (7.88 GB, 51.3 min), both HEVC Rext 2048x1536, bpps 0.81.
- **File to use:** NARA reels, H.264 1920x1080 @24 fps, about 50 Mb/s.
  - r1 630 s / 3.97 GB
  - r2 612 s / 3.85 GB
  - r3 801 s / 5.04 GB
  - r4 687 s / 4.32 GB
  - r5 592 s / 3.73 GB
  - bpps 3.04 on every reel
  - Total 55.4 min, 20.9 GB, about 4 GB per reel, so render one reel per session.
  - NARA's copy carries a 384 kb/s AAC track whose content is **[UNVERIFIED]** (silence or added
    music). The 24 fps container probably holds duplicated frames, so check cadence.
  - The Prelinger 2K is the higher-resolution option. Run a side-by-side on one shot before
    choosing.
- **Rights checked:** NARA `useRestriction: Undetermined`, `accessRestriction: Unrestricted`.
  Prelinger items carry no license field and are not dark. The film includes captured German
  footage of a U-boat sinking a freighter, per NARA's scope note. Foreign-origin inserts are the
  one place a federal-work argument could be contested **[UNVERIFIED]**.
- **Narration angle:** how two million Americans were drafted, shipped and fed into the Meuse-Argonne
  in eighteen months, told through the government's own propaganda film. The episode can be honest
  about what the CPI chose to show and what it cut. The content warning on Prelinger's item
  (images of the dead) needs an age-gate decision.
- **Shorts hooks:**
  - "This is the US Army's own film of WWI, the version they wanted you to see."
  - "In 1917 America had fewer than 200,000 soldiers. Watch what happened in one year."
- **Verdict:** safest rights position on the slate. Scan quality is good, with a choice of two
  sources. **WWI colorization is a crowded format** (Peter Jackson's *They Shall Not Grow Old* plus
  many AI uploads), so the differentiation has to be the narration. **Go**, but not as an early
  episode.
- The 200,000 figure in the second hook is from memory **[UNVERIFIED]**.

### 4. Rambles 'Round Rio (Kineto Review, 1921)

- **Year + evidence:** *Moving Picture World*, Dec 1921 (Lantern, `movpicwor532movi_0561`, read in
  this session). It lists "Rambles 'Round Rio — sightseeing with the Marines to the beauty spots
  around Brazil's capital" in the Kineto Company of America's fourth series of thirteen one-reel
  subjects, released through National Exchanges. The archive date (1921) agrees. Charles Urban's
  Kineto Review series ran 1921-23 (charlesurban.com via search).
- **Source:** `200626_Rambles_Round_Rio`. https://archive.org/details/200626_Rambles_Round_Rio
  (Prelinger / FFDW). Title card per description: "Visit to Tijuca, Sylvester, Nicheroy, and
  Paquata Island."
- **File to use:** `200626_Rambles_Round_Rio_master.intros.mov`, ProRes 4444 12-bit,
  **4096x3072 @16 fps, 841.4 s, 112.47 GB, bpps 10.62**. **There is no HEVC or 2K derivative**,
  only a 640x480 h.264. It is too big for any Kaggle disk, so it must be segmented straight from
  the URL (pipeline §3 pattern) or transcoded once to a local 2K/HEVC proxy.
- **Rights checked:** no license or rights fields. Not dark. Published 1921, so PD.
  A colorized "Rio de Janeiro 1920s in color … sightseeing tour" is already on YouTube
  (`watch?v=dRvpLq1yPqo`, via search). It may be this very film **[UNVERIFIED]**. Watch it before
  committing.
- **Narration angle:** US Marines on shore leave in 1921 Rio de Janeiro, then Brazil's capital, a
  year before its independence centennial exposition. The cable-car and boat sequences give a
  ride-along structure like Market Street. The narration covers the city's Belle Époque rebuilding
  (the Pereira Passos reforms) and what US sailors were doing in South America at the time.
- **Shorts hooks:**
  - "Rio de Janeiro, 1921, filmed from a cable car by American Marines."
  - "Before Christ the Redeemer existed: Rio's skyline with an empty Corcovado." The statue dates
    from 1931 [UNVERIFIED that Corcovado appears in the reel].
- **Verdict:** best raw scan on the slate and the rights are clean. The friction is the 112 GB
  master and a possible competing colorization. **Go if the competing upload is a different film.**
- Pereira Passos is named from memory **[UNVERIFIED]**.

### 5. U.S. Navy of 1915 (Lyman H. Howe Company)

- **Year + evidence:** NFPF's screening-room page and the archive description identify it as an
  11-minute surviving reel of Howe's three-reel *U.S. Navy of 1915*. Charles "Buckey" Grimm
  identified it partly from an E-2 class submarine taken out of service in 1915. Edge codes on the
  source print date from 1917. The National Film and Sound Archive of Australia lent the nitrate;
  prints are at the Library of Congress (NFPF page and doughboy.org via search).
- **Source:** `u-s-navy-documentary-1915`. https://archive.org/details/u-s-navy-documentary-1915
  (uploader `info@filmpreservation.org`).
- **File to use:** `US-Navy-1915-10mbps.mp4`, H.264 Main, **1438x1080 @24 fps, 730.0 s,
  0.94 GB, bpps 0.83**. This is an access encode at 10 Mb/s with a 301 kb/s AAC track (probably a
  score; content **[UNVERIFIED]**). It is the best copy on archive.org. LoC prints may exist at
  higher quality, but loc.gov was unreachable.
- **Rights checked:** `licenseurl` is the CC Public Domain Mark (NFPF's assertion). Released 1915,
  so PD if published, and Howe's road-show films were exhibited publicly. NFPF's page notes the
  Navy reserved reuse rights to officially approved films, which is not a restriction on us.
  **Strip the audio track** in case the score is a new copyrighted work **[UNVERIFIED]**.
- **Narration angle:** a "lost" film found in Australia. Secretary Josephus Daniels let a showman's
  camera crew aboard to sell a skeptical, neutral America on a bigger navy, months after the
  Lusitania. The episode covers what the sailors' daily life looked like and how many of those
  ships and men were at war two years later.
- **Shorts hooks:**
  - "This film was lost for a century. It turned up in Australia."
  - "1915: the US Navy hires a showman to convince America it needs a war fleet."
- **Verdict:** good story, clean rights, modest scan (access copy, not a master). **Go**, as a
  medium-effort episode.
- The Lusitania timing is my framing **[UNVERIFIED as a fact about the film]**.

### 6. Official Opening of the Los Angeles Aqueduct (event 5 Nov 1913)

- **Year + evidence:**
  - The event, Mulholland's "There it is. Take it.", took place on 5 November 1913 (Wikipedia, SCV
    History, waterandpower.org via search).
  - The archive item says `date 1916`, creator Southern California Edison.
  - **Calisphere's record (from a search snippet; the page itself would not load) gives "date
    created and/or issued 1913-02-13/1945".** A compilation date as late as 1945 cannot be
    excluded, so pre-1931 publication is **not established**.
- **Source:** `csmha_000346`, Huntington Library via California Revealed.
  https://archive.org/details/csmha_000346
- **File to use:** `csmha_000346_access.HD.mp4`, H.264 High, **1920x1440 @24 fps, 1231.5 s,
  1.54 GB, bpps 0.45**. This is an access copy, and the scout tool flags it as an upscale; see the
  tooling note at the top. True scan resolution is **[UNVERIFIED]**.
- **Rights checked:**
  - `rights: "Public domain. No restrictions on use."` and a CC PD Mark, both the institution's
    assertion.
  - Periscope Film sells the same event as stock footage (`17354 OPENING OF THE LOS ANGELES
    AQUEDUCT 1913`). Periscope-sourced footage may be in Content ID **[UNVERIFIED]**.
- **Narration angle:** the most consequential water project in American history and the true
  history behind *Chinatown*. It covers 233 miles from the Owens Valley, the city that bought a
  valley's water out from under it, and the Owens Valley dynamiting a decade later. Mulholland's
  St. Francis Dam collapse in 1928 is a ready-made epilogue.
- **Shorts hooks:**
  - "'There it is. Take it.' Five words that built Los Angeles and drained a valley."
  - "40,000 people came to watch water fall down a hill. Here's why."
- **Verdict:** story A+. Scan and date are both soft. Before production, find the Huntington
  finding-aid entry that dates the edit. **Conditional go.**

### 7. The Mutual Weekly No. 101 (newsreel, 1916)

- **Year + evidence:** the archive date is 1916. Segments per the catalogue: the Somme, Joffre at
  headquarters, Christmas turkey markets in Philadelphia, an NYC police safety drive, and a Los
  Angeles motorcycle hill climb. Subject headings name Jack London (d. 22 Nov 1916) and racing
  driver Lewis Jackson (d. 1916). Mutual Weekly was a released newsreel series, and the Christmas
  and Somme items put this issue at Dec 1916 **[issue date UNVERIFIED]**. Published 1916, so PD.
- **Source:** `copl_000132`, Oakland Public Library History Room via California Revealed.
  https://archive.org/details/copl_000132
- **File to use:** `copl_000132_access.HD.mp4`, H.264 High, **1920x1440 @16 fps, 1005.7 s,
  1.24 GB, bpps 0.45**. Access copy, same caveat as #6.
- **Rights checked:** `rights: "Public domain. No restrictions on use."` and a CC PD Mark. Not dark.
- **Narration angle:** one ordinary newsreel from the last Christmas before America entered the
  war, watched the way a 1916 audience saw it: the Somme next to turkey markets and motorcycle
  daredevils. It plays as a "time capsule" episode format that can recur (one newsreel, one
  week).
- **Shorts hooks:**
  - "What Americans saw at the movies the week Jack London died."
  - "December 1916: the Somme, Christmas turkeys and a motorcycle death-climb, in one newsreel."
- **Verdict:** clean rights, real scan of uncertain resolution, good recurring format. **Go.**
- The Jack London segment's content (funeral? ranch?) is **[UNVERIFIED]**. Watch before writing
  the first hook.

### 8. Allied Expeditionary Force in Russia (footage 1918-1920)

- **Year + evidence:** the archive item has no date. Description: "Documentary footage of military
  troops sent to Russia after World War I … Silent master reproduced from preprint material at
  National Archives and Records Administration". Can: "ALLIED EXPEDITIONARY FORCE IN RUSSIA …
  ORIG.". The AEF North Russia / Siberia interventions ran 1918-1920 (general knowledge
  **[UNVERIFIED here]**). The related NARA item it cites, naId 4523912 *AEF in Siberia*, is dated
  **1964** and marked "Restricted - Possibly". **This reel may be a later Army compilation of
  1918-20 footage.** If it is Army-made, it is still a federal work and PD. If any third-party
  newsreel material is inside, it is not covered.
- **Source:** `202197_Allied_Expeditionary_Force_in_Russia`.
  https://archive.org/details/202197_Allied_Expeditionary_Force_in_Russia
- **File to use:** `..._master.intros.mov`, HEVC Rext, **2048x1536 @16 fps (assumed from sibling
  files, not probed), 2876.3 s, 7.36 GB, bpps 0.81**. It fits Kaggle.
- **Rights checked:** no license or rights fields. Not dark. NARA flag as above.
- **Narration angle:** America's forgotten invasion of Russia. The "Polar Bears" were sent to
  Arkhangelsk and fought the Red Army in the winter after the Armistice. A second force sat in
  Siberia guarding a railway amid Czechs, Japanese and Cossacks. Almost nobody knows it happened,
  which is the hook.
- **Shorts hooks:**
  - "World War I ended. These American soldiers kept fighting, in Russia."
  - "In 1919 US troops fought the Red Army at 40 below. You weren't taught this."
- **Verdict:** excellent story, good scan. Rights hinge on confirming who compiled the reel.
  **Conditional go**, ranked here for that one open question.
- The second hook's temperature is from memory **[UNVERIFIED]**.

### 9. Railroads in United States History (Society for Visual Education, 1920)

- **Year + evidence:** the archive `rights` field transcribes the on-screen notice: "Copyrighted
  1920 by Society for Visual Education, Inc." The date is 1920. Published with notice in 1920,
  so PD.
- **Source:** `200320_Railroads_In_US_History`.
  https://archive.org/details/200320_Railroads_In_US_History (Prelinger / FFDW).
- **File to use:** `..._master.intros.mov`, ProRes 4444 12-bit, **4096x3072 @24 fps, 672.1 s,
  123.16 GB, bpps 14.56**. No HEVC derivative. Segment from the URL or make a one-time proxy.
- **Rights checked:** the copyright notice is quoted above. No license field. Not dark.
- **Narration angle:** an educational film's own 1920 account of American railroads, compared
  with what it left out (land grants, Chinese labour on the Central Pacific, the strikes). The
  "how they told it then vs. what we know now" frame is the episode.
- **Shorts hooks:**
  - "A 1920 school film explains railroads. Here's the part it skipped."
  - "In 1920 the railroads were America's biggest business. This is them at their peak."
- **Verdict:** superb scan and clean rights. The story is only as good as the footage, and an SVE
  classroom film may be maps and intertitles **[UNVERIFIED, watch first]**. The 123 GB master is
  a pipeline cost. **Reserve-to-go.**
- The second hook's claim is unsourced **[UNVERIFIED]**.

### 10. Ford Film Collection selects (Ford Motor Co. footage, late 1910s-1930s)

- **Year + evidence:** archive description: footage selected from NARA's Ford Film Collection by
  the producers of *The Atomic Café* (a 1982 film, from memory **[UNVERIFIED]**). It covers early
  car culture and Henry Ford, including Ford's acquisition of Lincoln from Henry Leland. A NARA
  search result for the same event, "Ford Family at Lincoln Purchase Ceremony", is dated 1921.
  **Shots are mixed 1910s-1930s and are not dated individually.**
- **Source:** `203684_National_Archives_Cars_Diego`.
  https://archive.org/details/203684_National_Archives_Cars_Diego
- **File to use:** `..._master.intros.mov`, HEVC Rext 2048x1536, **2447.8 s, 6.27 GB, bpps 0.81**.
- **Rights checked:** NARA *Prologue* (2014, "Henry Ford: Movie Mogul?") and the NARA card-catalog
  page, via search, say Ford transferred all rights to its 1.5 million feet of film to the National
  Archives (28 Nov 1962), and the films are "copyright- and royalty-free", with an exception for a
  few films made by Ford contractors. The archive item has no license field and is not dark.
  Because of that exception, cut to pre-1931 material anyway.
- **Narration angle:** Henry Ford ran his own film studio, then the largest industrial film
  operation in the world, and used it to sell the idea of the car itself. The episode covers
  Model T off-roading stunts, the Lincoln takeover and Ford as a media mogul.
- **Shorts hooks:**
  - "Henry Ford ran a movie studio. This is what he filmed."
  - "1920s car ads didn't show the road. They showed cars driving through rivers."
- **Verdict:** good rights position and a good scan. The compilation has no inherent narrative
  shape, so it needs a strong script. **Reserve-to-go.**
- The "largest" claim is unsourced **[UNVERIFIED]**.

---

## Reserves (qualify, ranked below the ten)

- **The Wizardry of Wireless** (GE, 1923). `the-wizardry-of-wireless-1923`, H.264 1440x1080,
  1331.7 s, 1.66 GB, bpps 0.80. NFPF, PD Mark. Covers the birth of broadcast radio and GE's WGY.
  Partly animated, so colorization adds less.
- **[Burton Holmes Travelogue]** (1920s). `204055_Burton_Holmes_Travelogue`, HEVC Rext 2048x1536,
  1231.7 s, 3.15 GB. Includes Hawaiian surf-riding. **Date unknown**, and Holmes copyrighted his
  travelogues, so the date must be established first.
- **Egypt (Part 3, Reel 2)** (1920s). `203152_Egypt_P3_R2`, 2K HEVC, 2.70 GB. Covers Cairo and
  Giza. Unidentified series, undated.

## Rejected, and why

- **Ellis Island** (`200455_Ellis_Island`). The scan is superb: HEVC 4096x3072, 365 s, 3.80 GB,
  plus a 34.7 GB ProRes at bpps 7.77. **It is undated**, described only as an "assembly of
  historical scenes", tagged 1920s/1930s, with no maker. A compilation could postdate 1930 even if
  the shots do not. It would rank in the top three if dated; date it first.
- **New York to San Francisco Through the Panama Canal, 1923** (`201347_…_012623`). This is the
  best scan found: **5120x3840** ProRes, 148 GB, plus an 18.5 GB HEVC. It is a **home movie by
  J. D. Zellerbach**, which means **unpublished**, so the "pre-1931" rule does not apply and
  life + 70 could run past 2030 **[Zellerbach's death year and the rights chain UNVERIFIED]**.
  The can notes mention a Dunningcolor lab, so it may not even be black-and-white
  **[UNVERIFIED]**.
- **The Olympic Games in Paris 1924** (`silent-the-olympic-games-in-paris-1924`). 1438x1080,
  14.6 GB, 174 min, bpps 0.90. The story is great (Nurmi, Weissmuller, *Chariots of Fire*). The
  source appears to be a modern restoration of unknown provenance, and IOC claims on Olympic
  footage are a known risk **[UNVERIFIED]**.
- **Marvels of Motion** (1930) and **The Fight for Life** (UFA, ©1930). 4K Prelinger. They are at
  the edge of the cutoff year, both have voiceover/sound tracks that compete with narration, and
  the stories are thin.
- **Lumière shorts at 1440x1080** (`lourdes-procession-*`, `duel-au-pistolet-1896`, etc.). The
  bpps is about 5, but the uploads are anonymous `opensource_movies` copies, very likely ripped from
  a commercial restoration, and each is under a minute. High claim risk, low narration room.
- **PatheFreres "Early Films 1-3"** (35mmstockfootage). 1080p at bpps about 10, but the files are
  fragments "in poor condition" and the item says to license HD from Getty Images, a claim
  magnet.
- **San Francisco Earthquake Aftermath, 1906** (Prelinger, several parts). Story A+, but only
  640x480 copies exist on archive.org. Revisit if LoC has a better transfer.
- **LoC Paper Print films** (Edison/Biograph actualities: NYC subway 1905, Pan-American Expo
  1901, Galveston 1900). **Not evaluated**: loc.gov blocked automated access (HTTP 403) in this
  session. Archive.org copies of these are SD or DeOldify re-uploads. This is the biggest
  unexplored source and worth one manual browser session.
- **Prelinger 4K fiction reels** (*The Covered Wagon* R4/R6, *Code of the Sea*, *Richthofen*,
  *Daddies*, *Her Boy Friend*, etc.). Genuine 4K scans, but single reels of narrative features
  (different format, weaker for narration). *Richthofen* (1927) also needs its publication
  checked.
- **French History of Aviation** (`202428_…`). Can notes cite a USAF number and the credits
  mention narration and music, so it is probably a post-war compilation. Undated.
- **Panama-Pacific / Kinograms / Gould-can items in SD** (640x480). Superseded by #1 or too
  low-resolution.
- **Everything in `colorized-movies` / DeOldify re-uploads and the 4K Disney/Mickey items**.
  Already colorized, or scout-confirmed upscales (bpps 0.01-0.03).

## Next actions

1. Watch the 640x480 derivatives of #1, #2, #4, #6 and #8 (about 90 min total) and write down
   what is actually on screen before scripting.
2. Date the Ellis Island assembly and the LA Aqueduct edit (Huntington finding aid); either could
   move up.
3. One manual browser session on loc.gov for 1900s paper-print actualities.
4. Fix `scout_films.py` so the bpps verdict depends on codec (ProRes vs HEVC vs H.264 access).
