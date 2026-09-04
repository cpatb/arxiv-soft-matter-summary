#!/usr/bin/env bash
# uytrend.sh <log> -- Uy initial residual every 400th solve, to estimate time-to-criterion.
grep "Solving for Uy," "$1" | awk -F'Initial residual = ' '{split($2,a,","); print NR, a[1]}' | awk 'NR%400==0'
