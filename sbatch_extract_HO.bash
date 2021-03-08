#!/bin/bash
#SBATCH -J 3droi
#SBATCH --time=0-01:00:00
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3G
#SBATCH -p ncf
#SBATCH --account=mclaughlin_lab
# Outputs ----------------------------------
#SBATCH -o log/%x-%A-%a.out

module load afni/18.3.05-fasrc01

set -aeuxo pipefail

DIRLIST=($( cat $1 ))
N=$SLURM_ARRAY_TASK_ID
SUBSESS=${DIRLIST[${N}]}
outfilepre=$(echo ${SUBSESS} | sed -e 's/.*\([0-9]\{4\}\)\/month\([0-9]\{2\}\).*/\1_\2/')

TASK="faceReactivity"
RUN="FaceReactivity1"
FUNC="${SUBSESS}/${TASK}/${RUN}_new_ssmooth.nii.gz"
HO=${SUBSESS}/${TASK}/HarvardOxford-sub-maxprob-thr50-2mm.nii.gz_${RUN}.nii.gz
MASK=${SUBSESS}/${TASK}/${RUN}_brainmask.nii.gz
HONEW=${SUBSESS}/${TASK}/HarvardOxford-sub-maxprob-thr50-2mm_${RUN}_masked.nii.gz

if [ ! -f ${FUNC} ]; then
	echo "File not found: ${FUNC}"
	exit 1
fi

3dcalc -a ${HO} -b ${MASK} -expr 'a*b' -prefix ${HONEW}

3dROIstats -mask ${HONEW} -1Dformat "${FUNC}" > HO/"${outfilepre}_TS.1D"
