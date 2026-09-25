# Claims — σ\*, rattling and contact geometry (2026-09-24 overnight)

Each claim: status, the one evidence file, and the referee state. Nothing here is from the bench.

| # | Claim | Status | Evidence | Referee |
|---|---|---|---|---|
| 1 | σ\*/κ = (N_c/V)⟨aᵢaⱼ⟩_c·C/(3qΛ) exactly (Love–Weber, Derjaguin threshold); Λ = 0.994–1.010 per snapshot on 52 sheared DEM cells | **standing** | `derivations/sigma_star_geometry_2026-09-24.md` eq. 7; `out/sigma_star_geometry.txt` | 4 passes, converged |
| 2 | At fixed number-mean size, CV 0 → 0.5 lowers σ\* ×0.73–0.79 and the contact density Zφ carries it (×0.60–0.71); force–size correlation C and tail shape 1/q push back 4–10 % each. The earlier "tail flattening" reading was not the mechanism for the solver's σ\* (Derjaguin, sheared) | **standing** (scoped) | same, §4–5, figure G1(a) | converged |
| 3 | Contact pairs are a configuration model in the grain degrees, ⟨aᵢaⱼ⟩_c = (Σkᵢaᵢ/Σkᵢ)², to ≤ 4 % in dry isotropic, dry sheared and lubricated (LF_DEM) packs; on 135 binary packs the SL count is 1.011× the random-pairing value overall, with a group-median excess up to 5.7 % (ratio 4) and SS/LL down to 0.95; a real 2–3 % disassortative tilt at CV 0.5 (finite-M null 0.9999) | **standing** | σ\* §5b; `derivations/rattler_closure_2026-09-24.md` §3b; `out/decomp_lfsig.txt` (post hoc) | converged |
| 4 | ⇒ σ\*/κ = Zφ ā_k²/(8π⟨a³⟩)·C/q: the PSD enters σ\* only through φ, Z, ⟨a³⟩ and the contact-weighted radius | **standing** (identity given 3) | σ\* eq. 15 | converged |
| 5 | Rattling follows surface area: Poisson counts with mean ∝ a² predict the rattler fraction with no fitted parameter, rms 0.016 out of sample on 702 isotropic records (size-blind 0.128); sharp optimum at p = 2; per class on 135 binaries rms 0.017 / 0.018 (size-blind 0.13 / 0.35); discrete K = 3–5 packs spanning ratio 4 (129 cells) rms 0.11–0.12, under-predicting by ~0.09 (size-blind 0.37–0.39) | **standing** | rattler_closure §4, `out/rattler_shape_compare.txt`; per grain `out/zbysize_check.txt` | 3 passes, converged |
| 6 | The frictional backbone solid fraction does not densify with CV (|Δφ_B| ≤ 0.015 with the a² rattler volume); frictionless backbones do (Δφ_B up to +0.03) | **standing, model-conditional** (ψ is modelled, not measured) | rattler_closure §4b, `out/backbone_psi_compare.txt` | converged |
| 7 | Closed chain Z_nr (length-free), φ, PSD → σ\*(CV): with the a² law for Z and one MEAN-degree exponent p = 1.625 measured on isotropic per-grain dumps, predicted/measured σ\* ratio 0.980–1.018 on every reliable-rung row; the exponent that zeroes each reliable CV-0.5 residual is 1.61–1.67; with p = 2 it captures only 35–58 % of the CV-0.5 drop. C/q still taken measured | **standing, provisional** (exponent from 2 isotropic cells, transferred; binaries give p ≈ 1.45 at ratio ≤ 1.5 but 2.4–2.9 at ratio 4–5, so p ≈ 1.6 is a moderate-polydispersity number, not a law) | σ\* §5c, `out/sigma_star_chain.txt`, `out/degree_exponent_probe.txt` | 6 passes, converged |
| 8 | `cv05margin`: D = φ_m(CV .5, μ .3) − φ₀(CV 0, μ 0) = +0.0018, t 0.99 (n = 9) — whether a weighed φ exists where CV 0 flows and CV 0.5 can shear-jam is below resolution | **UNDETERMINED** | `PREREG_CV05MARGIN_2026-09-24.md` RESULT, `out/grade_cv05margin.txt` | registered test |
| 9 | LF_DEM (lubricated) test of the Yade σ\*(CV) ratio (0.768) — P1–P5 | **INCOMPLETE** (13/28 cells; queue driver killed for low memory; relaunch needs author OK) | `PREREG_LFSIG_2026-09-24.md` Am. 1–4, `out/grade_lfsig.txt` | — |
| 10 | No published measurement of σ\* vs polydispersity at fixed repulsion exists (LF-DEM literature: Mari 2014, Singh 2018, Pednekar 2018 all fixed PSD or no repulsion) | standing, literature | citation agents' reports, DOIs in PREREG_LFSIG Am. 2 and `paper/refs.bib` | — |

| 11 | χ's CV-invariance, read through φ = φ_B/(1−ψ): friction stops the backbone absorbing the fines, which rattle instead; the extra rattler loading (+0.05 … +0.10 in ln φ at CV 0.5, sign and friction trend robust, magnitude closure-dependent) stands in for the frictionless backbone's densification. The near-cancellation itself is forced by the measured χ and is not a finding | **standing, model-conditional** | rattler_closure §4c, `out/backbone_psi_compare.txt` | converged |

**Withdrawn tonight:** "closure C (pair weight aᵢaⱼ(aᵢ+aⱼ)) is what rattling follows" (replaced by 5);
"the contact network is disassortative" (inferred across datasets, refuted in the same cell, now the small
tilt in 3); "the frictionless backbone carries 2/3 of the rise" (closure-dependent, 48–66 %).

**What would move this forward most:** (i) relaunch `lfsig` (`lfdem/relaunch_lfsig.sh`); (ii) one sheared
dump with grain ids to measure the mean degree law under shear (closes claim 7).
