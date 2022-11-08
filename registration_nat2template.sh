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

# flirt mean functional to t2 with skull
flirt -searchrx -360 360 -searchry -360 360 -searchrz -360 360 -in  ${temp_dir}/SEEPI_${p}_mean_reg_bold_${r}.nii.gz -ref ${anat_pth}/InplaneT2.nii.gz -out ${temp_dir}/${p}_bold_${r}.mean_to_t2.nii.gz -omat ${temp_dir}/${p}_bold_${r}.mean_to_t2.mat
flirt -in ${temp_dir}/errts.${p}_bold_${r}.tproject.nii.gz -ref ${anat_pth}/InplaneT2.nii.gz -applyxfm -init ${temp_dir}/${p}_bold_${r}.mean_to_t2.mat -out ${temp_dir}/errts.${p}_bold_${r}.tproject_to_t2.nii.gz -interp trilinear
if(! -f ${temp_dir}/mask_to_t2.nii.gz)then
flirt -searchrx -360 360 -searchry -360 360 -searchrz -360 360 -in ${anat_pth}/mask.nii.gz -ref ${anat_pth}/InplaneT2.nii.gz -out ${temp_dir}/mask_to_t2.nii.gz
endif
if(! -f ${temp_dir}/mask_to_t2_binary.nii.gz)then
3dcalc -a ${temp_dir}/mask_to_t2.nii.gz -expr 'ispositive(a)' -prefix ${temp_dir}/mask_to_t2_binary.nii.gz
endif
if(! -f ${temp_dir}/InplaneT2_masked.nii.gz)then
3dcalc -a ${anat_pth}/InplaneT2.nii.gz -float -b ${temp_dir}/mask_to_t2_binary.nii.gz -expr '(a*b)' -prefix ${temp_dir}/InplaneT2_masked.nii.gz
endif
if(! -f ${temp_dir}/t2_to_template.5iso_)then
antsRegistrationSyNQuick.sh -d 3 -f ${mskpth}/template_T2w_brain.5iso.nii.gz -m ${temp_dir}/InplaneT2_masked.nii.gz -o ${temp_dir}/t2_to_template.5iso_
endif