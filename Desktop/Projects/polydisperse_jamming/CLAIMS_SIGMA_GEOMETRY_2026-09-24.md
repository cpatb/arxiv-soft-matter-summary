# Claims — σ\*, rattling and contact geometry (2026-09-24 overnight)

Each claim: status, the one evidence file, and the referee state. Nothing here is from the bench.

| # | Claim | Status | Evidence | Referee |
|---|---|---|---|---|
| 1 | σ\*/κ = (N_c/V)⟨aᵢaⱼ⟩_c·C/(3qΛ) exactly (Love–Weber, Derjaguin threshold); Λ = 0.994–1.010 per snapshot on 52 sheared DEM cells | **standing** | `derivations/sigma_star_geometry_2026-09-24.md` eq. 7; `out/sigma_star_geometry.txt` | 4 passes, converged |
| 2 | At fixed number-mean size, CV 0 → 0.5 lowers σ\* ×0.73–0.79 and the contact density Zφ carries it (×0.60–0.71); force–size correlation C and tail shape 1/q push back 4–10 % each. The earlier "tail flattening" reading was not the mechanism for the solver's σ\* (Derjaguin, sheared) | **standing** (scoped) | same, §4–5, figure G1(a) | converged |
| 3 | Contact pairs are a configuration model in the grain degrees, ⟨aᵢaⱼ⟩_c = (Σkᵢaᵢ/Σkᵢ)², to ≤ 4 % in dry isotropic, dry sheared and lubricated (LF_DEM) packs; a real 2–3 % disassortative tilt at CV 0.5 (finite-M null 0.9999) | **standing** | σ\* §5b; `derivations/rattler_closure_2026-09-24.md` §3b; `out/decomp_lfsig.txt` (post hoc) | converged |
| 4 | ⇒ σ\*/κ = Zφ ā_k²/(8π⟨a³⟩)·C/q: the PSD enters σ\* only through φ, Z, ⟨a³⟩ and the contact-weighted radius | **standing** (identity given 3) | σ\* eq. 15 | converged |
| 5 | Rattling follows surface area: Poisson counts with mean ∝ a² predict the rattler fraction with no fitted parameter, rms 0.016 out of sample on 702 isotropic records (size-blind 0.128); sharp optimum at p = 2 | **standing** | rattler_closure §4, `out/rattler_shape_compare.txt`; per grain `out/zbysize_check.txt` | 3 passes, converged |
| 6 | The frictional backbone solid fraction does not densify with CV (|Δφ_B| ≤ 0.015 with the a² rattler volume); frictionless backbones do (Δφ_B up to +0.03) | **standing, model-conditional** (ψ is modelled, not measured) | rattler_closure §4b, `out/backbone_psi_compare.txt` | converged |
| 7 | Closed chain Z_nr, φ, PSD → σ\*(CV) gets Z to 5 % and CV 0.2 to 3–7 % but captures only 29–58 % of the CV-0.5 drop: the mean degree law is shallower than a² | **standing, open** | σ\* §5c, `out/sigma_star_chain.txt` | pass 5 pending |
| 8 | `cv05margin`: D = φ_m(CV .5, μ .3) − φ₀(CV 0, μ 0) = +0.0018, t 0.99 (n = 9) — whether a weighed φ exists where CV 0 flows and CV 0.5 can shear-jam is below resolution | **UNDETERMINED** | `PREREG_CV05MARGIN_2026-09-24.md` RESULT, `out/grade_cv05margin.txt` | registered test |
| 9 | LF_DEM (lubricated) test of the Yade σ\*(CV) ratio (0.768) — P1–P5 | **INCOMPLETE** (13/28 cells; queue driver killed for low memory; relaunch needs author OK) | `PREREG_LFSIG_2026-09-24.md` Am. 1–4, `out/grade_lfsig.txt` | — |
| 10 | No published measurement of σ\* vs polydispersity at fixed repulsion exists (LF-DEM literature: Mari 2014, Singh 2018, Pednekar 2018 all fixed PSD or no repulsion) | standing, literature | citation agents' reports, DOIs in PREREG_LFSIG Am. 2 and `paper/refs.bib` | — |

**Withdrawn tonight:** "closure C (pair weight aᵢaⱼ(aᵢ+aⱼ)) is what rattling follows" (replaced by 5);
"the contact network is disassortative" (inferred across datasets, refuted in the same cell, now the small
tilt in 3); "the frictionless backbone carries 2/3 of the rise" (closure-dependent, 48–66 %).

**What would move this forward most:** (i) relaunch `lfsig` (`lfdem/relaunch_lfsig.sh`); (ii) one sheared
dump with grain ids to measure the mean degree law under shear (closes claim 7).
