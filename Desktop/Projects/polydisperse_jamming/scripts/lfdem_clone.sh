#!/bin/bash
# Clone LF_DEM (Seto/Mari/Morris) for the sigma* campaign capability check. 2026-09-24.
set -u
mkdir -p ~/lfdem && cd ~/lfdem
[ -d LF_DEM ] || git clone -q https://github.com/ryseto/LF_DEM
cd LF_DEM && git log -1 --format='%H %ci' && ls
ls ~/miniforge3/bin | grep -E "^(conda|mamba)$"
~/miniforge3/bin/conda env list
