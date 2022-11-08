#!/bin/tcsh -xef
# user specifications
set apth = $1 #path to individual subject BOLD files
set anat_pth = $2 #path to individual subject T2 anatomical files
set corrpth = $3 #directory for correlation (& regressor) files
set mskpth = $4 #path to template masks, and anatomical files.
set output_dir = $5 #path out output files
set temp_dir = $6 #temp directory, iterim files later deleted to save space
set p = $7
set r = $8 

if(! -f ${temp_dir}/errts.${p}_bold_${r}.tproject_to_template.nii.gz)then
antsApplyTransforms -e 3 -i ${temp_dir}/errts.${p}_bold_${r}.tproject_to_t2.nii.gz -r ${temp_dir}/t2_to_template.5iso_Warped.nii.gz -o ${temp_dir}/errts.${p}_bold_${r}.tproject_to_template.nii.gz -t ${temp_dir}/t2_to_template.5iso_0GenericAffine.mat -t ${temp_dir}/t2_to_template.5iso_1Warp.nii.gz
endif
3dcalc -a ${temp_dir}/errts.${p}_bold_${r}.tproject_to_template.nii.gz -b ${mskpth}/gray_template.5iso.nii.gz -expr '(a*b)' -prefix ${corrpth}/${p}_bold_${r}_to_template_errts.5iso_masked_gm.nii.gz

#create nuisance regressors
3dmaskave -quiet -mask ${mskpth}/white_matter.5iso.nii ${temp_dir}/errts.${p}_bold_${r}.tproject_to_template.nii.gz > ${corrpth}/${p}_bold_${r}_wm_nui_regressor.1D
3dmaskave -quiet -mask ${mskpth}/csf.5iso.nii ${temp_dir}/errts.${p}_bold_${r}.tproject_to_template.nii.gz > ${corrpth}/${p}_bold_${r}_csf_nui_regressor.1D
paste ${corrpth}/${p}_bold_${r}_wm_nui_regressor.1D ${corrpth}/${p}_bold_${r}_csf_nui_regressor.1D > ${corrpth}/${p}_bold_${r}_nui_regressors.1D

rm ${corrpth}/${p}_bold_${r}_wm_nui_regressor.1D
rm ${corrpth}/${p}_bold_${r}_csf_nui_regressor.1D