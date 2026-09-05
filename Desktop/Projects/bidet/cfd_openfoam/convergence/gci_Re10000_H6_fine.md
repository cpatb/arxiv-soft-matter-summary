# Grid convergence, impinging jet at Re_d = 1e4, H/d = 6 (`--family fine`)

Three uniformly refined meshes, refinement ratios r21 = 2.000, r32 = 2.000.

| grid | refine | cells | n_t | C_P | tau_at_15mm_Pa | tau_peak_Pa |
|---|---|---|---|---|---|---|
| 1 (fine) | 4 | 491520 | 1.9338 | 1.5809 | 6.1374 | 486.3032 |
| 2  | 2 | 122880 | 1.9119 | 1.5178 | 6.2448 | 454.4698 |
| 3  | 1 | 30720 | 1.9023 | 1.3518 | 5.8658 | 389.5219 |

## Iterative convergence

A GCI measures discretisation error and assumes iterative error is far smaller. Residual target in `fvSolution` is 1e-05.

| grid | iterations | `SIMPLE converged` | worst final residual | fields |
|---|---|---|---|---|
| refine 4 | 11552 | yes | **1e-05** | Ux 2.6e-06, Uy 1e-05, Uz 4.9e-06, k 5.5e-06, omega 9.5e-08, p 1.6e-08 |
| refine 2 | 3000 | **no (cap)** | **3e-05** ⚠ above target | Ux 7.4e-06, Uy 3e-05, Uz 2.1e-05, k 2e-05, omega 7.6e-07, p 2.9e-08 |
| refine 1 | 2261 | yes | **1e-05** | Ux 7.6e-07, Uy 1e-05, Uz 1.6e-06, k 2.8e-06, omega 3.7e-08, p 8e-09 |

## Discretisation error

| quantity                | phi1_fine | phi2   | phi3   | p     | phi_extrap | e_a21_% | e_ext21_% | GCI_fine_% | verdict                                                                                                                                                                                                                                                   |
|-------------------------|-----------|--------|--------|-------|------------|---------|-----------|------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| closure exponent n_t    | 1.9338    | 1.9119 | 1.9023 | --    | --         | 1.13    | --        | --         | NOT CONVERGENT: DIVERGING under refinement: |eps21/eps32| = 2.279 >= 1, i.e. the gap between the two FINEST grids is larger than between the two coarsest. Celik et al. exclude this case; the procedure would report |p| = 1.188 for a true p of -1.188. |
| closure prefactor C_P   | 1.5809    | 1.5178 | 1.3518 | 1.397 | 1.6195     | 3.99    | 2.38      | 3.05       | OK                                                                                                                                                                                                                                                        |
| tau_w at r = 15 mm [Pa] | 6.1374    | 6.2448 | 5.8658 | --    | --         | 1.75    | --        | --         | NOT CONVERGENT: non-monotone in the grid (sign of eps32/eps21 is negative); Richardson extrapolation invalid                                                                                                                                              |
| peak tau_w [Pa]         | 486.3     | 454.47 | 389.52 | 1.029 | 516.9      | 6.55    | 5.92      | 7.87       | OK                                                                                                                                                                                                                                                        |

Procedure: Celik et al., J. Fluids Eng. 130(7):078001 (2008); Fs = 1.25 for three grids.
