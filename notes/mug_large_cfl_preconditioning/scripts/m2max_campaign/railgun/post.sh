#!/bin/bash
# post.sh <case_dir>: plotting run (3 ranks) + plot_railgun.py sheet-speed analysis
set -e
G=${OFT_TREE:?set OFT_TREE to the build tree parent (holds mpich-*/ build_railgun/)}
cd "$1"
sed 's/^ plot_run=F / plot_run=T /' oft_coax.in > oft_coax_plot.in
cp oft_coax_plot.in oft_coax_plotrun.in
# plot_railgun.py reads oft_coax.in for parameters, plotting run uses its own copy
OMP_NUM_THREADS=1 GFORTRAN_UNBUFFERED_PRECONNECTED=y $G/mpich-4_2_3/bin/mpirun -np 3 $G/build_railgun/examples/MUG/railgun/MUG-railgun oft_coax_plot.in oft_coax.xml > plot_log.txt 2>&1
mkdir -p figs
MPLBACKEND=Agg $G/../../venv/bin/python plot_railgun.py --coax . --outdir figs > sheet_speed.txt 2>&1
cat sheet_speed.txt
