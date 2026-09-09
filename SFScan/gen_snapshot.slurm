#!/bin/bash
#
#SBATCH --job-name=2DLandscape_Snapshot
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

module load gromacs/2021.3

# Set some environment variables 
LANDSCAPE_DIR=`pwd`
echo "2D Landscape home directory set to ${LANDSCAPE_DIR}"
FREE_ENERGY=/work/e280/e280-Howlin/cve280/FEP/
echo "Free energy home directory set to ${FREE_ENERGY}"

ION=Na

cd ${LANDSCAPE_DIR}/${ION}

echo 0 | gmx_mpi trjconv -s ${FREE_ENERGY}/${ION}/MD_REP2_FWD/Lambda_29/md29.gro -f ${FREE_ENERGY}/${ION}/MD_REP2_FWD/Lambda_29/md29.xtc -dt 100 -o snapshot.gro -sep

for ((i=0; i<6; i++))
do
    echo "QM_LINK" | gmx_mpi editconf -f snapshot${i}.gro -n ${LANDSCAPE_DIR}/index_qmmm.ndx -o REP${i}/sf_link${i}.gro
done

ION=Na_HW

cd ${LANDSCAPE_DIR}/${ION}

echo 0 | gmx_mpi trjconv -s ${FREE_ENERGY}/${ION}/MD_REP2_FWD/Lambda_29/md29.gro -f ${FREE_ENERGY}/${ION}/MD_REP2_FWD/Lambda_29/md29.xtc -dt 100 -o snapshot.gro -sep

for ((i=0; i<6; i++))
do
    echo "QM_LINK" | gmx_mpi editconf -f snapshot${i}.gro -n ${LANDSCAPE_DIR}/index_qmmm.ndx -o REP${i}/sf_link${i}.gro
done

ION=Li

cd ${LANDSCAPE_DIR}/${ION}

echo 0 | gmx_mpi trjconv -s ${FREE_ENERGY}/${ION}/MD_REP2_FWD/Lambda_29/md29.gro -f ${FREE_ENERGY}/${ION}/MD_REP2_FWD/Lambda_29/md29.xtc -dt 100 -o snapshot.gro -sep

for ((i=0; i<6; i++))
do
    echo "QM_LINK" | gmx_mpi editconf -f snapshot${i}.gro -n ${LANDSCAPE_DIR}/index_qmmm.ndx -o REP${i}/sf_link${i}.gro
done

ION=Li_HW

cd ${LANDSCAPE_DIR}/${ION}

echo 0 | gmx_mpi trjconv -s ${FREE_ENERGY}/${ION}/MD_REP2_FWD/Lambda_29/md29.gro -f ${FREE_ENERGY}/${ION}/MD_REP2_FWD/Lambda_29/md29.xtc -dt 100 -o snapshot.gro -sep

for ((i=0; i<6; i++))
do
    echo "QM_LINK" | gmx_mpi editconf -f snapshot${i}.gro -n ${LANDSCAPE_DIR}/index_qmmm.ndx -o REP${i}/sf_link${i}.gro
done

ION=K

cd ${LANDSCAPE_DIR}/${ION}

echo 0 | gmx_mpi trjconv -s ${FREE_ENERGY}/${ION}/MD_REP2_FWD/Lambda_29/md29.gro -f ${FREE_ENERGY}/${ION}/MD_REP2_FWD/Lambda_29/md29.xtc -dt 100 -o snapshot.gro -sep

for ((i=0; i<6; i++))
do
    echo "QM_LINK" | gmx_mpi editconf -f snapshot${i}.gro -n ${LANDSCAPE_DIR}/index_qmmm.ndx -o REP${i}/sf_link${i}.gro
done

ION=K_HW

cd ${LANDSCAPE_DIR}/${ION}

echo 0 | gmx_mpi trjconv -s ${FREE_ENERGY}/${ION}/MD_REP2_FWD/Lambda_29/md29.gro -f ${FREE_ENERGY}/${ION}/MD_REP2_FWD/Lambda_29/md29.xtc -dt 100 -o snapshot.gro -sep

for ((i=0; i<6; i++))
do
    echo "QM_LINK" | gmx_mpi editconf -f snapshot${i}.gro -n ${LANDSCAPE_DIR}/index_qmmm.ndx -o REP${i}/sf_link${i}.gro
done

exit;
