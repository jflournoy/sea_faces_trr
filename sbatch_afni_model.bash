#!/bin/bash
#SBATCH -J 3ddeconfit 
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
#N=0
SUBSESS=${DIRLIST[${N}]}
TASK="faceReactivity"
RUN="FaceReactivity1"
FUNC="${SUBSESS}/${TASK}/${RUN}_new_ssmooth.nii.gz"
MASK=${SUBSESS}/${TASK}/${RUN}_brainmask.nii.gz
ICCBASE=/ncf/mclaughlin/stressdevlab/stress_pipeline/Group/FaceReactivity/3dICC/bicc
AFNIOUT=${ICCBASE}/afni_out

outfilepre=$(echo ${SUBSESS} | sed -e 's/.*\([0-9]\{4\}\)\/month\([0-9]\{2\}\).*/\1_\2/')

if [ ! -f ${FUNC} ]; then
	echo "File not found: ${FUNC}"
	exit 1
fi

EVS=( C F H )
for EV in ${EVS[@]}; do
	awk '{print $1}' ${SUBSESS}/${TASK}/evfiles/${RUN}_${EV}.txt | tr '\n' ' ' > ${SUBSESS}/${TASK}/evfiles/${RUN}_${EV}.1D
done


# generate model matrix using 3dDeconvolve for each of the subject
if [ ! -f ${AFNIOUT}/${outfilepre}_IM.x1D ]; then
    3dDeconvolve -input ${FUNC} -mask ${MASK} \
    	-polort A -nobout -local_times -jobs 1 \
    	-num_stimts 3 \
    	-stim_times_IM 1 ${SUBSESS}/${TASK}/evfiles/${RUN}_C.1D 'SPMG2(18)' -stim_label 1 Calm \
    	-stim_times_IM 2 ${SUBSESS}/${TASK}/evfiles/${RUN}_F.1D 'SPMG2(18)' -stim_label 2 Fear \
    	-stim_times_IM 3 ${SUBSESS}/${TASK}/evfiles/${RUN}_H.1D 'SPMG2(18)' -stim_label 3 Happy \
    	-ortvec ${SUBSESS}/${TASK}/${RUN}_nuisance_regressors.txt 'nuisance' \
    	-cbucket ${AFNIOUT}/${outfilepre}_betas_IM.nii.gz \
    	-x1D ${AFNIOUT}/${outfilepre}_IM.x1D \
    	-x1D_stop -xsave
else
    echo "${AFNIOUT}/${outfilepre}_IM.x1D exsits"
fi


## Run GLS+ARMA(1,1) using 3dREMLfit by looping through 57 subjects and 11 ROIs
#sch400 1-400
ROIs=($(seq 0 399))
for ROI in ${ROIs[@]}; do
    if [ ! -f ${AFNIOUT}/${outfilepre}_${ROI}_REMLfit.1D ]; then
        3dREMLfit -input ${ICCBASE}/sch400/${outfilepre}_TS.1D[${ROI}]\' \
        -matrix ${AFNIOUT}/${outfilepre}_IM.x1D -noFDR\
        -Grid 5 -tout -Rbuck ${AFNIOUT}/${outfilepre}_${ROI}_REMLfit.1D
    else
        echo "${AFNIOUT}/${outfilepre}_${ROI}_REMLfit.1D exsits"
    fi
done

#HO 1-21
ROIs=($(seq 0 20))
for ROI in ${ROIs[@]}; do
    if [ ! -f ${AFNIOUT}/HO_${outfilepre}_${ROI}_REMLfit.1D ]; then
        3dREMLfit -input ${ICCBASE}/HO/${outfilepre}_TS.1D[${ROI}]\' \
        -matrix ${AFNIOUT}/${outfilepre}_IM.x1D -noFDR\
        -Grid 5 -tout -Rbuck ${AFNIOUT}/HO_${outfilepre}_${ROI}_REMLfit.1D
    else
        echo "${AFNIOUT}/HO_${outfilepre}_${ROI}_REMLfit.1D exists"
    fi
done