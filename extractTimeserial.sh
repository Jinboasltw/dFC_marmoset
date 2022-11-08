#!/bin/tcsh -xef

# user specifications
set p = $1
set r = $2
set temp_dir = $3
set corrpth = $4 #correlation path
set mskpth = $5 #path to mask

if (-f ${corrpth}/${p}_bold_${r}_to_template_errts.5iso_masked_gm.nii.gz)then
#create output directory if it doesn't already exist
if (! -d ${corrpth}/${p}_${r}_split/)then
mkdir ${corrpth}/${p}_${r}_split/
endif
endif

#calculate the correlation between each voxel and every other voxel, output is each voxel (seed) and resultant connectivity as a volume (ordered by coordinate)
if (-f ${corrpth}/${p}_bold_${r}_nui_regressors.1D)then
if(! -f ${corrpth}/${p}_${r}_corrmap.nii.gz)then
3dTcorrMap -ort ${corrpth}/${p}_bold_${r}_nui_regressors.1D -input ${corrpth}/${p}_bold_${r}_to_template_errts.5iso_masked_gm.nii.gz -mask ${mskpth} -Corrmap ${corrpth}/${p}_${r}_corrmap.nii.gz
endif
endif

endif
endif