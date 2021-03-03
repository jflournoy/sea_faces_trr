#!/bin/bash
#SBATCH -J ants
#SBATCH --time=0-01:00:00
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3G
#SBATCH -p ncf
#SBATCH --account=mclaughlin_lab
# Outputs ----------------------------------
#SBATCH -o log/%x-%A-%a.out

source /users/jflournoy/code/FSL-6.0.4.txt
source /users/jflournoy/code/ANTS-2.3.1.txt

############ Sch400 >> custom >> T1 >> func
set -aeuxo pipefail
DIRLIST=($( cat $1 ))
N=$SLURM_ARRAY_TASK_ID
SUBSESS=${DIRLIST[${N}]}

echo "Transforming Schaefer for ${SUBSESS}..."

SCH400="/ncf/mclaughlin/stressdevlab/stress_pipeline/Group/FaceReactivity/3dICC/bicc/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz"
TASK="faceReactivity"
RUN="FaceReactivity1"
FUNCMID="${SUBSESS}/${TASK}/${RUN}_FinalMidVol.nii.gz"
OUTFN=$( basename ${SCH400} )
OUTFN=${OUTFN%.nii.gz}_${RUN}.nii.gz
OUTIMAGE="${SUBSESS}/${TASK}/${OUTFN}"

if [ ! -f ${FUNCMID} ]; then
	echo "File not found: $FUNCMID"
	exit 1
fi

# Brain templates
PROJECT_DIR=/ncf/mclaughlin/stressdevlab/stress_pipeline
STANDARD_DIR=${PROJECT_DIR}/Template
STD_BRAIN=${FSLDIR}/data/standard/MNI152_T1_2mm_brain.nii.gz
MNI_BRAIN=${FSLDIR}/data/standard/MNI152_T1_2mm_brain.nii.gz
ST_BRAIN=${STANDARD_DIR}/ST_brain.nii.gz
ST_REG_PREFIX=${STANDARD_DIR}/ST_brain_to_MNI_brain

pushd $SUBSESS

antsApplyTransforms -i ${SCH400} -d 3 -e 0 \
       -n MultiLabel \
	-r ${FUNCMID} \
	-t [xfm_dir/${TASK}/${RUN}_to_T1_s_0GenericAffine.mat,1] \
	xfm_dir/${TASK}/${RUN}_to_T1_s_1InverseWarp.nii.gz \
	[xfm_dir/T1_to_custom_s_0GenericAffine.mat,1] \
	xfm_dir/T1_to_custom_s_1InverseWarp.nii.gz \
	[${ST_REG_PREFIX}_0GenericAffine.mat,1] \
	${ST_REG_PREFIX}_1InverseWarp.nii.gz \
	-v \
	-o ${OUTIMAGE}

popd
