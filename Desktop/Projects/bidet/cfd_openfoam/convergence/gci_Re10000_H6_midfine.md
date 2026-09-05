# Grid convergence, impinging jet at Re_d = 1e4, H/d = 6 (`--family midfine`)

Three uniformly refined meshes, refinement ratios r21 = 1.414, r32 = 1.414.

| grid | refine | cells | n_t | C_P | tau_at_15mm_Pa | tau_peak_Pa |
|---|---|---|---|---|---|---|
| 1 (fine) | 2 | 122880 | 1.8922 | 1.3928 | 6.2415 | 441.1746 |
| 2  | 1.414 | 61009 | 1.9237 | 1.4098 | 5.7726 | 417.1974 |
| 3  | 1 | 30720 | 1.9023 | 1.3518 | 5.8658 | 389.5219 |

## Iterative convergence

A GCI measures discretisation error and assumes iterative error is far smaller. Residual target in `fvSolution` is 1e-05.

| grid | iterations | `SIMPLE converged` | worst final residual | fields |
|---|---|---|---|---|
| refine 2 | 5105 | yes | **1e-05** | Ux 8.7e-07, Uy 1e-05, Uz 3.3e-06, k 5.2e-06, omega 1.3e-08, p 9.7e-09 |
| refine 1.414 | 20000 | **no (cap)** | **0.018** ⚠ above target | Ux 2.9e-09, Uy 2.4e-09, Uz 0.018, k 9.7e-09, omega 5.7e-09, p 8.3e-09 |
| refine 1 | 2261 | yes | **1e-05** | Ux 7.6e-07, Uy 1e-05, Uz 1.6e-06, k 2.8e-06, omega 3.7e-08, p 8e-09 |

## Discretisation error

| quantity                | phi1_fine | phi2   | phi3   | p     | phi_extrap | e_a21_% | e_ext21_% | GCI_fine_% | verdict                                                                                                      |
|-------------------------|-----------|--------|--------|-------|------------|---------|-----------|------------|--------------------------------------------------------------------------------------------------------------|
| closure exponent n_t    | 1.8922    | 1.9237 | 1.9023 | --    | --         | 1.66    | --        | --         | NOT CONVERGENT: non-monotone in the grid (sign of eps32/eps21 is negative); Richardson extrapolation invalid |
| closure prefactor C_P   | 1.3928    | 1.4098 | 1.3518 | --    | --         | 1.22    | --        | --         | NOT CONVERGENT: non-monotone in the grid (sign of eps32/eps21 is negative); Richardson extrapolation invalid |
| tau_w at r = 15 mm [Pa] | 6.2415    | 5.7726 | 5.8658 | --    | --         | 7.51    | --        | --         | NOT CONVERGENT: non-monotone in the grid (sign of eps32/eps21 is negative); Richardson extrapolation invalid |
| peak tau_w [Pa]         | 441.17    | 417.2  | 389.52 | 0.416 | 595.54     | 5.43    | 25.92     | 43.74      | p = 0.42 OUTSIDE [0.5, 2.5]: not in the asymptotic range, GCI is indicative only                             |

Procedure: Celik et al., J. Fluids Eng. 130(7):078001 (2008); Fs = 1.25 for three grids.
