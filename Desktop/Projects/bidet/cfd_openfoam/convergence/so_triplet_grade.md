# PREREG_so_triplet grade

## Health (Amendment 2)

| case | iter | worst non-Uz | Uz resid | H1 | max|Uz|/max|U| | H2 | orig. residualControl |
|---|---|---|---|---|---|---|---|
| case_Re10000_H6 | 2261 | 9.99e-06 | 1.60e-06 | PASS | 2.8e-11 | PASS | met |
| case_Re10000_H6_fo_r1p5 | 191 | 9.86e-06 | 1.33e-03 | PASS | 2.0e-11 | PASS | NOT met |
| case_Re10000_H6_fo_r1p5_cont | 591 | 7.81e-06 | 2.72e-04 | PASS | 2.0e-11 | PASS | NOT met |
| case_Re10000_H6_fo_r2p25 | 1290 | 5.31e-06 | 7.83e-04 | PASS | 2.0e-14 | PASS | NOT met |
| case_Re10000_H6_fo_r2p25_cont | 1690 | 5.76e-06 | 8.71e-04 | PASS | 2.0e-14 | PASS | NOT met |
| case_Re10000_H6_so_r1 | 1948 | 9.67e-06 | 1.22e-06 | PASS | 3.8e-11 | PASS | met |
| case_Re10000_H6_so_r1p5 | 3392 | 1.89e-05 | 6.88e-06 | FAIL | 3.6e-11 | PASS | NOT met |
| case_Re10000_H6_so_r2p25 | 765 | 3.19e-06 | 1.11e-03 | PASS | 1.2e-13 | PASS | NOT met |
| case_Re10000_H6_so_r2p25_cont | 1165 | 9.45e-06 | 9.09e-04 | PASS | 1.2e-13 | PASS | NOT met |
| case_Re10000_H6_so_r3p375 | 29 | 6.15e-06 | 4.95e-04 | PASS | 1.1e-13 | PASS | NOT met |

## C1 (Amendment 4): gate quantities against iteration, +2400 in steps of 400

- **case_Re10000_H6_fo_r1p5** (7 writes, t = 400..2591): last-two max 1.191% (bar 0.1%), +1600->+2400 max 2.249% (bar 0.2%), Uz ratio 2.4e-11 -> NOT CONVERGED
  | t | n_t | C_P | tau_at_15mm_Pa | tau_peak_Pa |
  |---|---|---|---|---|
  | 400 | 1.9109 | 1.3808 | 5.848 | 394.8 |
  | 800 | 1.9123 | 1.3779 | 5.8116 | 398 |
  | 1200 | 1.9158 | 1.3771 | 5.7641 | 402.95 |
  | 1600 | 1.9218 | 1.3863 | 5.709 | 407.89 |
  | 2000 | 1.9258 | 1.3944 | 5.6805 | 412.3 |
  | 2400 | 1.9291 | 1.4023 | 5.6634 | 415.92 |
  | 2591 | 1.9304 | 1.406 | 5.658 | 417.27 |
- case_Re10000_H6_fo_r2p25: NOT GRADED (0 trajectory points)
- **case_Re10000_H6_so_r1p5** (7 writes, t = 3600..5792): last-two max 0.025% (bar 0.1%), +1600->+2400 max 0.139% (bar 0.2%), Uz ratio 3.6e-11 -> CONVERGED
  | t | n_t | C_P | tau_at_15mm_Pa | tau_peak_Pa |
  |---|---|---|---|---|
  | 3600 | 1.9249 | 1.4427 | 5.9231 | 485.63 |
  | 4000 | 1.9313 | 1.4539 | 5.859 | 485.74 |
  | 4400 | 1.9352 | 1.4618 | 5.8218 | 485.83 |
  | 4800 | 1.9371 | 1.4661 | 5.8043 | 485.89 |
  | 5200 | 1.9377 | 1.4678 | 5.7987 | 485.93 |
  | 5600 | 1.9378 | 1.4682 | 5.7986 | 485.96 |
  | 5792 | 1.9377 | 1.4682 | 5.7994 | 485.96 |
- **case_Re10000_H6_so_r2p25** (7 writes, t = 800..3165): last-two max 0.277% (bar 0.1%), +1600->+2400 max 0.430% (bar 0.2%), Uz ratio 3.0e-15 -> NOT CONVERGED
  | t | n_t | C_P | tau_at_15mm_Pa | tau_peak_Pa |
  |---|---|---|---|---|
  | 800 | 1.9376 | 1.474 | 5.8388 | 480.44 |
  | 1200 | 1.9394 | 1.4762 | 5.8181 | 481.53 |
  | 1600 | 1.9407 | 1.4773 | 5.7975 | 482.63 |
  | 2000 | 1.9418 | 1.4786 | 5.7895 | 483.45 |
  | 2400 | 1.9428 | 1.4805 | 5.7787 | 484.2 |
  | 2800 | 1.9437 | 1.4823 | 5.7711 | 484.93 |
  | 3165 | 1.9442 | 1.4835 | 5.768 | 485.54 |
- case_Re10000_H6_so_r3p375: NOT GRADED (0 trajectory points)

## S6 (superseded by C1; reported only): +400 continuation of mapped-start cases

- case_Re10000_H6_fo_r2p25: NOT GRADED (base or continuation missing/unhealthy)
- case_Re10000_H6_so_r2p25: NOT GRADED (base or continuation missing/unhealthy)
- case_Re10000_H6_fo_r1p5: NOT GRADED (base or continuation missing/unhealthy)

## Triplets

**so_r15: INCOMPLETE** -- missing/unhealthy ['case_Re10000_H6_so_r2p25'], C1 not converged ['case_Re10000_H6_so_r2p25']

**fo_r15: INCOMPLETE** -- missing/unhealthy ['case_Re10000_H6_fo_r2p25', 'case_Re10000_H6_fo_r1p5'], C1 not converged ['case_Re10000_H6_fo_r2p25', 'case_Re10000_H6_fo_r1p5']

**so_r15_S7: INCOMPLETE** -- missing/unhealthy ['case_Re10000_H6_so_r3p375', 'case_Re10000_H6_so_r2p25'], C1 not converged ['case_Re10000_H6_so_r3p375', 'case_Re10000_H6_so_r2p25']

## Gates

- S1, S2: INCOMPLETE
- S3: INCOMPLETE
- S7: INCOMPLETE
- **S4** FO refine 1.5 health: H1 PASS (worst non-Uz 9.86e-06); original residualControl NOT met
- **S5** |dC_P|/C_P at refine 1, SO vs FO = 5.18% -> PASS (bar 1%)
