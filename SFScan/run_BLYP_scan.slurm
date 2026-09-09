#!/bin/bash

#SBATCH --time=01-00:00:00
#SBATCH --account=e280-Howlin
#SBATCH --partition=standard
#SBATCH --qos=taskfarm
#SBATCH -o %j_%x.o
#SBATCH -e %j_%x.e

#SBATCH --distribution=block:block
#SBATCH --hint=nomultithread

#SBATCH --nodes=1
#SBATCH --tasks-per-node=64
#SBATCH --cpus-per-task=2

#SBATCH -J 2DScan_CP2K_BLYP_Li_REP0
#SBATCH --exclusive

echo $SLURM_NODELIST

# define global paths
export HOME="/home/e280/e280/cve280"
export MWSHPATH=$HOME/MShTools
export WORK="/work/e280/e280-Howlin/cve280"
export GROMAX=$WORK/GromacsOnEUREKA
export ROOT=$WORK/2D_Landscape

cd $ROOT

# define inputs & outputs
DIRECTORY='Li/REP0'
TEMPLATE='template_BLYP_CP2K.inp'
ION='Li'
DIELECTRIC='7.65000000E+001' #7.65000000E+001 for H2O; 7.60000000E+001 for D2O


module load cp2k/cp2k-2022.1
# module load cray-python/3.9.7.1

function RunSinglePoint() {

	# module load cp2k/cp2k-2022.1

	cd $ROOT/$DIRECTORY

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
	cp $ROOT/$TEMPLATE $OUTKEY.inp

	sed -i "s/__PROJECT__/$OUTKEY/" $OUTKEY.inp
	sed -i "s|__COORD_FILE_NAME__|../$OUTKEY.pdb|" $OUTKEY.inp
	sed -i "s|__ION__|$ION|" $OUTKEY.inp
	sed -i "s|__DIELECTRIC__|$DIELECTRIC|" $OUTKEY.inp

	# OMP setup
	export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
	export OMP_PLACES=cores

	# actual QM/MM SMD
	echo "Running cp2k.psmp... $OUTKEY.inp > $OUTKEY.out"
	time srun cp2k.psmp -i $OUTKEY.inp > $OUTKEY.out

	STATUS=$?
	if [ $STATUS -ne 0 ] ; then
		echo "srun cp2k.psmp failed"
		return $STATUS
	fi
	
	if [ $( grep -n "SCF run NOT converged" $OUTKEY.out | wc -l) -ne 0 ] ; then
		echo 'SCF did not converge!'

		# extract the QMMM energies (should be only one or two)
		grep -n "ENERGY|" $OUTKEY.out | awk '{print $1" "$10}' > $DATA_FILE.non_converged

		STATUS=-2

	else

		# extract the QMMM energies (should be only one or two)
		grep -n "ENERGY|" $OUTKEY.out | awk '{print $1" "$10}' > $DATA_FILE

	fi

	rm *.wfn*

	return $STATUS
}

# go through each pdb configuration in the input folder

for PDB_INPUT in $DIRECTORY/grid_???.pdb; do

	# run a super short QM/MM for each one
	RunSinglePoint $PDB_INPUT

	echo -e "\033[1mSTATUS="$?"\033[0m"

done

exit $STATUS
