#!/bin/bash

#SBATCH --time=00-00:20:00
#SBATCH --account=e280-Howlin
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH -o %j_%x.o
#SBATCH -e %j_%x.e

#SBATCH --distribution=block:block
#SBATCH --hint=nomultithread

#SBATCH --nodes=2
#SBATCH --tasks-per-node=64
#SBATCH --cpus-per-task=2

#SBATCH -J CP2K_QMMM_test_Na_Implicit

cd $SLURM_SUBMIT_DIR
echo $SLURM_NODELIST

# Define variables
LANDSCAPE_DIR=`pwd`
echo "2D Landscape home directory set to ${LANDSCAPE_DIR}"
DIR=${LANDSCAPE_DIR}/Na
MDP='gromacs_qmmm.mdp'
TOPOL='topol.top'
GRO_COORD='snapshot0.gro'
INDEX='index_qmmm.ndx'
GRO_OUTPUT='sod_gro_qmmm_test0'

module use /work/y07/shared/archer2-lmod/training/
module load gromacs/2021.1+cp2k

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores

cd ${DIR}

# Generate PDB and INP files using grompp
time gmx_cp2k grompp -f ../${MDP} -p ${TOPOL} -c ${GRO_COORD} -n ../${INDEX} -o ${GRO_OUTPUT}.tpr -maxwarn 1

sed -n '1,106794p;106795q' ${GRO_OUTPUT}.pdb > ${GRO_OUTPUT}_no_water.pdb 

# define inputs & outputs
TEMPLATE='template_B3LYP_QMMM_NW_CP2K.inp'
PROJECT='Na_SF_QMMM_Implicit_Test'
INPUT_FILE='sod_test0_qmmm_implicit.inp'
OUTPUT_FILE='test0_qmmm_implicit.out'
ION='Na'

cd ${LANDSCAPE_DIR}

# copy over the template input
cp ${TEMPLATE} ${DIR}/${INPUT_FILE}

cd ${DIR}

echo -e "\033[1m""RunSinglePoint" $1 "\033[0m"

sed -i "s/__PROJECT__/${PROJECT}/" ${INPUT_FILE}
sed -i "s|__COORD_FILE_NAME__|${GRO_OUTPUT}_no_water.pdb|" ${INPUT_FILE}
sed -i "s|__ION__|${ION}|" ${INPUT_FILE}

grep -B1000 "&END DFT" ${INPUT_FILE} > before.inp
grep -A1000 "&TOPOLOGY" ${INPUT_FILE} > after.inp

sed -n '47,145p;146q' ${GRO_OUTPUT}.inp > middle.inp 

cat before.inp > final.inp
cat middle.inp >> final.inp
cat after.inp >> final.inp

mv final.inp ${INPUT_FILE}

module load cp2k/cp2k-2022.1
# module load cray-python/3.9.7.1

# actual QM/MM SMD
echo "Running cp2k.psmp... ${INPUT_FILE} > ${OUTPUT_FILE}"
time srun cp2k.psmp -i ${INPUT_FILE} > ${OUTPUT_FILE}

exit 0
