#!/bin/bash
# Locations used by the run, scoring and status scripts. Every value can be overridden with an environment variable,
# e.g.  BENCH=/data/MICA-2026-benchmark RUNS_ROOT=/scratch/runs bash scripts/run_mica_batch.sh
#
# BENCH        benchmark root: this repository with the Zenodo data extracted into it (dataset/, af3/, runs/, eval/)
# METHODS_DIR  where the five methods are cloned (MICA, EModelX, model-angelo, CryoAtom, EMProt)
# RUNS_ROOT    per-method run folders; use a fast local Linux disk (not a Windows drive under WSL)
# PHENIX_ACT   Phenix 1.20.1 environment script
# USALIGN      US-align executable
# CONDA_SH     conda shell hook
BENCH=${BENCH:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
METHODS_DIR=${METHODS_DIR:-$HOME/bench/methods}
RUNS_ROOT=${RUNS_ROOT:-$HOME/bench/runs}
PHENIX_ACT=${PHENIX_ACT:-$HOME/phenix/phenix-1.20.1-4487/phenix_env.sh}
USALIGN=${USALIGN:-$HOME/bench/tools/USalign/USalign}
CONDA_SH=${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
