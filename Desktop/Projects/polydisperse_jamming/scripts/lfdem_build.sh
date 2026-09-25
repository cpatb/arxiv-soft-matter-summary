#!/bin/bash
# Build LF_DEM against conda-forge SuiteSparse in a dedicated env (does not touch base/yade). 2026-09-24.
set -u
source ~/miniforge3/etc/profile.d/conda.sh
if ! conda env list | grep -q "^lfdem "; then
  mamba create -y -q -n lfdem -c conda-forge suitesparse openblas cxx-compiler make > ~/lfdem/conda_create.log 2>&1 || { tail -20 ~/lfdem/conda_create.log; exit 1; }
fi
conda activate lfdem
P=$CONDA_PREFIX
ls $P/include/suitesparse/cholmod.h $P/include/cholmod.h 2>/dev/null
ls $P/lib/libcholmod* | head
cd ~/lfdem/LF_DEM/LF_DEM
INC=$P/include; [ -f $P/include/suitesparse/cholmod.h ] && INC=$P/include/suitesparse
cat > config/Makefile_config.mk <<CFG
CXX = g++
override_default_cholmod = true
Cholmod_path = -I $INC
Cholmod_Linking_Flags = -L$P/lib -lcholmod -lamd -lcolamd -lsuitesparseconfig -Wl,-rpath=$P/lib
Blas_Linking_Flags = -lopenblas
Lapack_Linking_Flags =
install_dir = ~/lfdem/bin
CFG
make -j8 > ~/lfdem/make.log 2>&1; echo "make rc=$?"
tail -15 ~/lfdem/make.log
ls -la LF_DEM && ./LF_DEM 2>&1 | head -20
