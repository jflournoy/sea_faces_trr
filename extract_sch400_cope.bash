#!/bin/bash
#/ContrastName1  Calm > Fear 
#/ContrastName2  Calm > Happy 
#/ContrastName3  Calm > Scramble 
#/ContrastName4  Fear > Calm 
#/ContrastName5  Fear > Happy 
#/ContrastName6  Fear > Scramble 
#/ContrastName7  Happy > Calm 
#/ContrastName8  Happy > Fear 
#/ContrastName9  Happy > Scramble 
#/ContrastName10 Calm 
#/ContrastName11 Fear 
#/ContrastName12 Happy 
#/ContrastName13 Scramble 

files=($(cat /mnt/stressdevlab/stress_pipeline/Group/FaceReactivity/NewNeuropoint/datafiles/setfilenames_fearGTcalm.txt))

for file in ${files[@]}; do
	outfilepre=$(echo $file | sed -e 's/.*\([0-9]\{4\}\)\/month\([0-9]\{2\}\).*/\1_\2/')
	3dROIstats -mask Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz -numROI 400 -1DRformat "${file}" > sch400_cope/"${outfilepre}_fearGTcalm.tsv"
done
