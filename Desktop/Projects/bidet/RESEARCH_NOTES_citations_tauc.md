# Research notes — verified citations & tau_c (2026-07-25)

All fields verified via Crossref direct DOI lookup. None fabricated.

## Flagged citations (Stage 3.5) — VERIFIED, remove `% FLAG FOR AUTHOR: verify`

```latex
\bibitem{BhagatWilson2016} R.~K.~Bhagat and D.~I.~Wilson, Flow in the thin film
created by a coherent turbulent water jet impinging on a vertical wall,
Chemical Engineering Science \textbf{152}, 606 (2016).
% DOI 10.1016/j.ces.2016.06.011 ; pages 606--623

\bibitem{Wilson2014} D.~I.~Wilson, P.~Atkinson, H.~K\"ohler, M.~Mauermann,
H.~Stoye, K.~Suddaby, T.~Wang, J.~F.~Davidson, and J.-P.~Majschak, Cleaning of
soft-solid soil layers on vertical and horizontal surfaces by stationary
coherent impinging liquid jets, Chemical Engineering Science \textbf{109}, 183 (2014).
% DOI 10.1016/j.ces.2014.01.034 ; pages 183--196 ; 9 authors
```

## tau_c — Stage 4 (genuine wall-detachment shear FOUND)

**Best match, in-band (1–10 Pa):** Nejadnik, van der Mei, Busscher & Norde (2008),
parallel-plate flow-displacement chamber, wall shear at attachment/detachment balance
for bacteria on silicone rubber:
- S. epidermidis 2.7±1.1 Pa; S. aureus 1.0±0.2 Pa; P. aeruginosa 1.1±0.2 Pa
- Genuine hydrodynamic WALL detachment shear (not bulk yield stress).
- Appl. Environ. Microbiol. 74(3):916–919, DOI 10.1128/AEM.01557-07

**Upper bracket:** Gião & Keevil (2013), removal of L. monocytogenes biofilm from
steel/PTFE, ~24–144 Pa. J. Appl. Microbiol. 114(1):256–265, DOI 10.1111/jam.12032

```latex
\bibitem{Nejadnik2008} M.~R.~Nejadnik, H.~C.~van der Mei, H.~J.~Busscher, and
W.~Norde, Determination of the shear force at the balance between bacterial
attachment and detachment in weak-adherence systems, using a flow displacement
chamber, Applied and Environmental Microbiology \textbf{74}, 916 (2008).

\bibitem{Giao2013} M.~S.~Gi\~{a}o and C.~W.~Keevil, Hydrodynamic shear stress to
remove Listeria monocytogenes biofilms from stainless steel and
polytetrafluoroethylene surfaces, Journal of Applied Microbiology \textbf{114},
256 (2013).
```

**Recommendation:** cite Nejadnik (1–3 Pa) as a genuine wall-detachment measurement
supporting tau_c ~ O(1) Pa, bracketed above by Gião & Keevil (~24–144 Pa). Keep the
exact tau_c an explicit assumption. Do NOT cite the Wilson soft-solid impinging-jet
papers as a tau_c source (their governing stress is a bulk viscoplastic property —
same conflation as the Woolley bulk-rheology cite we are replacing). Levers are ratios,
Q* ∝ tau_c^{1/3}, so weakly sensitive.
