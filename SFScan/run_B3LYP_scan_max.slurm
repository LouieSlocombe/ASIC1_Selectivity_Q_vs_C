#!/bin/bash

#SBATCH --time=01-00:00:00
#SBATCH --account=e280-Sacchi
#SBATCH --partition=standard
#SBATCH --qos=taskfarm
#SBATCH -o %j_%x.o
#SBATCH -e %j_%x.e

#SBATCH --distribution=block:block
#SBATCH --hint=nomultithread

#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --cpus-per-task=4

#SBATCH -J diab_B3LYP_2CM_mid_HL_N624A14

# define global paths
export WORK="/work/e280/e280-Sacchi/max"
export MWSHPATH=$WORK/MShTools
export GROMAX=$WORK/GromacsOnEUREKA
export ROOT=$WORK/DiabaticScan

cd $ROOT

# ../../ is DiabaticScan
# ../../../ArcherUmbrellaHelicase is ArcherUmbrellaHelicase

# # define inputs & outputs
# INPUT_DIRECTORY='../2CM_grid_HL06'
# OUTPUT_DIRECTORY='B3LYP_2CM_HL06'
# TEMPLATE='template_B3LYP_new.inp'

# TEMPLATE='template_B3LYP_new.inp'

# # define inputs & outputs
# INPUT_DIRECTORY='../2DM_grid_mid_L05'
# OUTPUT_DIRECTORY='B3LYP_2DM_mid_L05'
# TEMPLATE='template_B3LYP_2DM_new_c12.inp'

# define inputs & outputs
INPUT_DIRECTORY='../2CM_grid_HL_N624A14'
OUTPUT_DIRECTORY='B3LYP_2CM_mid_HL_N624A14'
TEMPLATE='template_B3LYP_N624A_new_c12.inp'

module load cp2k/cp2k-2022.1
# module load cray-python/3.9.7.1

function RunSinglePoint() {

	# module load cp2k/cp2k-2022.1

	cd $ROOT/$OUTPUT_DIRECTORY

	echo -e "\033[1m""RunSinglePoint" $1 "\033[0m"

	PDB_FILE=$1

	OUTKEY=$(basename "${PDB_FILE%.pdb}")

	DATA_FILE=$OUTKEY".energy"

	# check if the output datafile already exists and skip if so
	if test -f "$DATA_FILE"; then
		echo "Found existing data, skipping"
		return -1
	fi

	# make a subdirectory
	mkdir -pv $OUTKEY

	cd $OUTKEY

	# copy over the template input
	cp ../../$TEMPLATE $OUTKEY.inp

	sed -i "s/__PROJECT__/$OUTKEY/" $OUTKEY.inp
	sed -i "s|__COORD_FILE_NAME__|../$PDB_FILE|" $OUTKEY.inp

	# OMP setup
	export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
	export OMP_PLACES=cores

	# actual QM/MM SMD
	echo "Running cp2k.psmp... $OUTKEY.inp > $OUTKEY.out"
	time srun cp2k.psmp -i $OUTKEY.inp > $OUTKEY.out

	STATUS=$?
	if [ $STATUS -ne 0 ] ; then
		echo "srun cp2k.psmp failed"
		exit 2
	fi
	
	if [ $( grep -n "SCF run NOT converged" $OUTKEY.out | wc -l) -ne 0 ] ; then
		echo 'SCF did not converge!'

		# extract the QMMM energies (should be only one or two)
		grep -n "ENERGY|" $OUTKEY.out | awk '{print $1" "$10}' > "../"$DATA_FILE.non_converged

		STATUS=-2

	else

		# extract the QMMM energies (should be only one or two)
		grep -n "ENERGY|" $OUTKEY.out | awk '{print $1" "$10}' > "../"$DATA_FILE

	fi

	rm *.wfn*

	return $STATUS
}

# prepare output directory

mkdir -p $OUTPUT_DIRECTORY
cd $OUTPUT_DIRECTORY

# go through each gro configuration in the input folder

for PDB_INPUT in $INPUT_DIRECTORY/window_??_??.pdb; do

	# run a super short QM/MM for each one
	RunSinglePoint $PDB_INPUT

	echo -e "\033[1mSTATUS="$?"\033[0m"

done

# for JSON_INPUT in $INPUT_DIRECTORY/*.json; do

# 	module load cray-python/3.9.7.1

# 	export MOLPARSE=$WORK/MolParse
# 	export MWPYTPATH=$WORK/MPyTools
# 	PYTHONPATH=$PYTHONPATH:$MWPYTPATH
# 	PYTHONPATH=$PYTHONPATH:$MOLPARSE
# 	export PYTHONPATH

# 	echo $JSON_INPUT $INPUT_DIRECTORY/window_00_00.pdb ${JSON_INPUT%.json}.pdb

# 	# create the PDB input
# 	../shift.py $JSON_INPUT $INPUT_DIRECTORY/window_00_00.pdb ${JSON_INPUT%.json}.pdb

# 	# run a super short QM/MM for each one
# 	RunSinglePoint ${JSON_INPUT%.json}.pdb

# 	rm -v ${JSON_INPUT%.json}.pdb

# 	echo -e "\033[1mSTATUS="$?"\033[0m"

# done

exit $STATUS
