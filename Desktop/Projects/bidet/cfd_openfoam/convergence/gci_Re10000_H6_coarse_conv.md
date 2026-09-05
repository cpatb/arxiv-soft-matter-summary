# Grid convergence, impinging jet at Re_d = 1e4, H/d = 6 (`--family coarse_conv`)

Three uniformly refined meshes, refinement ratios r21 = 2.000, r32 = 2.000.

| grid | refine | cells | n_t | C_P | tau_at_15mm_Pa | tau_peak_Pa |
|---|---|---|---|---|---|---|
| 1 (fine) | 2 | 122880 | 1.8922 | 1.3928 | 6.2415 | 441.1746 |
| 2  | 1 | 30720 | 1.9023 | 1.3518 | 5.8658 | 389.5219 |
| 3  | 0.5 | 7680 | 1.8871 | 1.3008 | 5.8893 | 341.6188 |

## Iterative convergence

A GCI measures discretisation error and assumes iterative error is far smaller. Residual target in `fvSolution` is 1e-05.

| grid | iterations | `SIMPLE converged` | worst final residual | fields |
|---|---|---|---|---|
| refine 2 | 5105 | yes | **1e-05** | Ux 8.7e-07, Uy 1e-05, Uz 3.3e-06, k 5.2e-06, omega 1.3e-08, p 9.7e-09 |
| refine 1 | 2261 | yes | **1e-05** | Ux 7.6e-07, Uy 1e-05, Uz 1.6e-06, k 2.8e-06, omega 3.7e-08, p 8e-09 |
| refine 0.5 | 823 | yes | **9.9e-06** | Ux 1e-06, Uy 9.9e-06, Uz 6e-06, k 2.6e-06, omega 3.4e-08, p 9.7e-09 |

## Discretisation error

| quantity                | phi1_fine | phi2   | phi3   | p     | phi_extrap | e_a21_% | e_ext21_% | GCI_fine_% | verdict                                                                                                                                                                                                                                                   |
|-------------------------|-----------|--------|--------|-------|------------|---------|-----------|------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| closure exponent n_t    | 1.8922    | 1.9023 | 1.8871 | --    | --         | 0.53    | --        | --         | NOT CONVERGENT: non-monotone in the grid (sign of eps32/eps21 is negative); Richardson extrapolation invalid                                                                                                                                              |
| closure prefactor C_P   | 1.3928    | 1.3518 | 1.3008 | 0.315 | 1.5608     | 2.94    | 10.76     | 15.08      | p = 0.31 OUTSIDE [0.5, 2.5]: not in the asymptotic range, GCI is indicative only                                                                                                                                                                          |
| tau_w at r = 15 mm [Pa] | 6.2415    | 5.8658 | 5.8893 | --    | --         | 6.02    | --        | --         | NOT CONVERGENT: non-monotone in the grid (sign of eps32/eps21 is negative); Richardson extrapolation invalid                                                                                                                                              |
| peak tau_w [Pa]         | 441.17    | 389.52 | 341.62 | --    | --         | 11.71   | --        | --         | NOT CONVERGENT: DIVERGING under refinement: |eps21/eps32| = 1.078 >= 1, i.e. the gap between the two FINEST grids is larger than between the two coarsest. Celik et al. exclude this case; the procedure would report |p| = 0.109 for a true p of -0.109. |

Procedure: Celik et al., J. Fluids Eng. 130(7):078001 (2008); Fs = 1.25 for three grids.
