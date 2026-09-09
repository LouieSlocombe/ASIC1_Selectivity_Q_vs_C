#!/bin/bash

#SBATCH --time=00-00:20:00
#SBATCH --account=e280-Howlin
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH -o %j_%x.o
#SBATCH -e %j_%x.e

#SBATCH --distribution=block:block
#SBATCH --hint=nomultithread

#SBATCH --nodes=1
#SBATCH --tasks-per-node=64
#SBATCH --cpus-per-task=2

#SBATCH -J CP2K_BLYP_test_Na

cd $SLURM_SUBMIT_DIR
echo $SLURM_NODELIST

LANDSCAPE=`pwd`
echo "2D Landscape home directory set to ${LANDSCAPE}"
DIR=${LANDSCAPE}/Na

# define inputs & outputs
PROJECT='Na_SF_QM_Test2'
INPUT_FILE='test2_blyp_sod.inp'
OUTPUT_FILE='test2_blyp_sod.out'
COORD_FILE='test_sod.pdb'
TEMPLATE='template_BLYP_CP2K.inp'
ION='Na'

module load cp2k/cp2k-2022.1
# module load cray-python/3.9.7.1

# copy over the template input
cp ${TEMPLATE} ${INPUT_FILE}

echo -e "\033[1m""RunSinglePoint" $1 "\033[0m"

sed -i "s/__PROJECT__/$PROJECT/" $INPUT_FILE
sed -i "s|__COORD_FILE_NAME__|$COORD_FILE|" $INPUT_FILE
sed -i "s|__ION__|$ION|" $INPUT_FILE

# OMP setup
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores

# actual QM/MM SMD
echo "Running cp2k.psmp... ${INPUT_FILE} > ${OUTPUT_FILE}"
time srun cp2k.psmp -i $INPUT_FILE > $OUTPUT_FILE

exit 0
