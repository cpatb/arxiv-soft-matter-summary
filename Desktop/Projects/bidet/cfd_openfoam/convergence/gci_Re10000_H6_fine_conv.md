# Grid convergence, impinging jet at Re_d = 1e4, H/d = 6 (`--family fine_conv`)

Three uniformly refined meshes, refinement ratios r21 = 2.000, r32 = 2.000.

| grid | refine | cells | n_t | C_P | tau_at_15mm_Pa | tau_peak_Pa |
|---|---|---|---|---|---|---|
| 1 (fine) | 4 | 491520 | 1.9338 | 1.5809 | 6.1374 | 486.3032 |
| 2  | 2 | 122880 | 1.8922 | 1.3928 | 6.2415 | 441.1746 |
| 3  | 1 | 30720 | 1.9023 | 1.3518 | 5.8658 | 389.5219 |

## Iterative convergence

A GCI measures discretisation error and assumes iterative error is far smaller. Residual target in `fvSolution` is 1e-05.

| grid | iterations | `SIMPLE converged` | worst final residual | fields |
|---|---|---|---|---|
| refine 4 | 11552 | yes | **1e-05** | Ux 2.6e-06, Uy 1e-05, Uz 4.9e-06, k 5.5e-06, omega 9.5e-08, p 1.6e-08 |
| refine 2 | 5105 | yes | **1e-05** | Ux 8.7e-07, Uy 1e-05, Uz 3.3e-06, k 5.2e-06, omega 1.3e-08, p 9.7e-09 |
| refine 1 | 2261 | yes | **1e-05** | Ux 7.6e-07, Uy 1e-05, Uz 1.6e-06, k 2.8e-06, omega 3.7e-08, p 8e-09 |

## Discretisation error

| quantity                | phi1_fine | phi2   | phi3   | p     | phi_extrap | e_a21_% | e_ext21_% | GCI_fine_% | verdict                                                                                                      |
|-------------------------|-----------|--------|--------|-------|------------|---------|-----------|------------|--------------------------------------------------------------------------------------------------------------|
| closure exponent n_t    | 1.9338    | 1.8922 | 1.9023 | --    | --         | 2.15    | --        | --         | NOT CONVERGENT: non-monotone in the grid (sign of eps32/eps21 is negative); Richardson extrapolation invalid |
| closure prefactor C_P   | 1.5809    | 1.3928 | 1.3518 | 2.198 | 1.6333     | 11.90   | 3.21      | 4.14       | OK                                                                                                           |
| tau_w at r = 15 mm [Pa] | 6.1374    | 6.2415 | 5.8658 | --    | --         | 1.69    | --        | --         | NOT CONVERGENT: non-monotone in the grid (sign of eps32/eps21 is negative); Richardson extrapolation invalid |
| peak tau_w [Pa]         | 486.3     | 441.17 | 389.52 | 0.195 | 798.47     | 9.28    | 39.10     | 80.24      | p = 0.19 OUTSIDE [0.5, 2.5]: not in the asymptotic range, GCI is indicative only                             |

Procedure: Celik et al., J. Fluids Eng. 130(7):078001 (2008); Fs = 1.25 for three grids.
