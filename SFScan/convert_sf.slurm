#!/bin/bash
#
#SBATCH --job-name=Convert_SF_GRO2PDB
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

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores

module load cray-python/3.9.4.1

# Set some environment variables 
LANDSCAPE_DIR=`pwd`
echo "2D Landscape home directory set to ${LANDSCAPE_DIR}"
DIR=Na

cd ${LANDSCAPE_DIR}

for ((i=0; i<6; i++))
do
    python3 convert_link.py -d ${DIR}/REP${i} -s sf_link${i}.gro -o sf${i}.pdb
done

DIR=Na_HW

for ((i=0; i<6; i++))
do
    python3 convert_link.py -d ${DIR}/REP${i} -s sf_link${i}.gro -o sf${i}.pdb
done

DIR=Li

for ((i=0; i<6; i++))
do
    python3 convert_link.py -d ${DIR}/REP${i} -s sf_link${i}.gro -o sf${i}.pdb
done

DIR=Li_HW

for ((i=0; i<6; i++))
do
    python3 convert_link.py -d ${DIR}/REP${i} -s sf_link${i}.gro -o sf${i}.pdb
done

DIR=K

for ((i=0; i<6; i++))
do
    python3 convert_link.py -d ${DIR}/REP${i} -s sf_link${i}.gro -o sf${i}.pdb
done

DIR=K_HW

for ((i=0; i<6; i++))
do
    python3 convert_link.py -d ${DIR}/REP${i} -s sf_link${i}.gro -o sf${i}.pdb
done

exit;
