# Grid convergence, impinging jet at Re_d = 1e4, H/d = 6 (`--family coarse`)

Three uniformly refined meshes, refinement ratios r21 = 2.000, r32 = 2.000.

| grid | refine | cells | n_t | C_P | tau_at_15mm_Pa | tau_peak_Pa |
|---|---|---|---|---|---|---|
| 1 (fine) | 2 | 122880 | 1.9119 | 1.5178 | 6.2448 | 454.4698 |
| 2  | 1 | 30720 | 1.9023 | 1.3518 | 5.8658 | 389.5219 |
| 3  | 0.5 | 7680 | 1.8871 | 1.3008 | 5.8893 | 341.6188 |

## Iterative convergence

A GCI measures discretisation error and assumes iterative error is far smaller. Residual target in `fvSolution` is 1e-05.

| grid | iterations | `SIMPLE converged` | worst final residual | fields |
|---|---|---|---|---|
| refine 2 | 3000 | **no (cap)** | **3e-05** ⚠ above target | Ux 7.4e-06, Uy 3e-05, Uz 2.1e-05, k 2e-05, omega 7.6e-07, p 2.9e-08 |
| refine 1 | 2261 | yes | **1e-05** | Ux 7.6e-07, Uy 1e-05, Uz 1.6e-06, k 2.8e-06, omega 3.7e-08, p 8e-09 |
| refine 0.5 | 823 | yes | **9.9e-06** | Ux 1e-06, Uy 9.9e-06, Uz 6e-06, k 2.6e-06, omega 3.4e-08, p 9.7e-09 |

## Discretisation error

| quantity                | phi1_fine | phi2   | phi3   | p     | phi_extrap | e_a21_% | e_ext21_% | GCI_fine_% | verdict                                                                                                                                                                                                                                                   |
|-------------------------|-----------|--------|--------|-------|------------|---------|-----------|------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| closure exponent n_t    | 1.9119    | 1.9023 | 1.8871 | 0.657 | 1.9286     | 0.50    | 0.86      | 1.09       | OK                                                                                                                                                                                                                                                        |
| closure prefactor C_P   | 1.5178    | 1.3518 | 1.3008 | --    | --         | 10.94   | --        | --         | NOT CONVERGENT: DIVERGING under refinement: |eps21/eps32| = 3.256 >= 1, i.e. the gap between the two FINEST grids is larger than between the two coarsest. Celik et al. exclude this case; the procedure would report |p| = 1.703 for a true p of -1.703. |
| tau_w at r = 15 mm [Pa] | 6.2448    | 5.8658 | 5.8893 | --    | --         | 6.07    | --        | --         | NOT CONVERGENT: non-monotone in the grid (sign of eps32/eps21 is negative); Richardson extrapolation invalid                                                                                                                                              |
| peak tau_w [Pa]         | 454.47    | 389.52 | 341.62 | --    | --         | 14.29   | --        | --         | NOT CONVERGENT: DIVERGING under refinement: |eps21/eps32| = 1.356 >= 1, i.e. the gap between the two FINEST grids is larger than between the two coarsest. Celik et al. exclude this case; the procedure would report |p| = 0.439 for a true p of -0.439. |

Procedure: Celik et al., J. Fluids Eng. 130(7):078001 (2008); Fs = 1.25 for three grids.
