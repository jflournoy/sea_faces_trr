#!/bin/bash
3dUndump -ijk -datum 'float' -prefix TRR_meta_est.nii.gz -master Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz TRR_meta_est.txt
3dUndump -ijk -datum 'float' -prefix TRR_meta_ll.nii.gz -master Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz TRR_meta_ll.txt
3dUndump -ijk -datum 'float' -prefix TRR_meta_uu.nii.gz -master Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii.gz TRR_meta_uu.txt
