#!/bin/bash
#
#SBATCH --job-name=2DLandscape_Grid
#SBATCH --time=0-12:00:00 ##0-00:20:00 # 2-00:00:00 1-00:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32 ##128
#SBATCH --cpus-per-task=1 ##1
#SBATCH --account=e280-howlin
#SBATCH --partition=serial
#SBATCH --qos=serial ## short long standard
#SBATCH -o %j.%x.o 
#SBATCH -e %j.%x.e

cd $SLURM_SUBMIT_DIR
echo $SLURM_NODELIST

module load cray-python/3.9.7.1
module load tk

for ((i=0; i<1; i++))
do
    python3 gen_grid.py -d Na/REP${i} -s sf${i}.pdb -c 1.625 -r 5.125 -f sf_grid${i}.html -o grid >> grid_sod2.txt
    python3 gen_grid.py -d Na_HW/REP${i} -s sf${i}.pdb -c 1.625 -r 5.125 -f sf_grid${i}.html -o grid >> grid_sod_hw2.txt
    python3 gen_grid.py -d Li/REP${i} -s sf${i}.pdb -c 1.5 -r 5 -f sf_grid${i}.html -o grid >> grid_lit2.txt
    python3 gen_grid.py -d Li_HW/REP${i} -s sf${i}.pdb -c 1.5 -r 5 -f sf_grid${i}.html -o grid >> grid_lit_hw2.txt
    # python3 gen_grid.py -d K/REP${i} -s sf${i}.pdb -c 1.75 -r 5.25 -f sf_grid${i}.html -o scan >> grid_pot.txt
    # python3 gen_grid.py -d K_HW/REP${i} -s sf${i}.pdb -c 1.75 -r 5.25 -f sf_grid${i}.html -o scan >> grid_pot_hw.txt
done

exit;
