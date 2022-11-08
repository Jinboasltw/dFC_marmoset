#!/bin/tcsh -xef
set wd = $1
set p = $2
set r = $3
set tr_counts = 502
set tr = 2
set n_cpu = 6
# compute de-meaned motion parameters (for use in regression)
1d_tool.py -infile "${wd}/dfile.${p}_bold_${r}_est.1D" -set_nruns 1 -demean -write "${wd}/${p}_bold_${r}_motion_demean.1D"

# compute motion parameter derivatives (for use in regression)
1d_tool.py -infile "${wd}/dfile.${p}_bold_${r}_est.1D" -set_nruns 1 -derivative -demean -write "${wd}/motion_deriv_${p}_bold_${r}.1D"

#  create censor file, for censoring motion
1d_tool.py -infile "${wd}/dfile.${p}_bold_${r}_est.1D" -set_nruns 1 -show_censor_count -censor_prev_TR -censor_motion .5 "${wd}/motion_${p}_bold_${r}"

# create bandpass regressors
1dBport -nodata ${tr_counts} ${tr} -band 0.01 0.1 -invert -nozero > "${wd}/bandpass_rall_${p}_bold_${r}.1D"

set ktrs = `1d_tool.py -infile ${wd}\/motion_${p}_bold_${r}_censor.1D -show_trs_uncensored encoded`

# run the regression analysis
3dDeconvolve -input ${wd}/pb04.bold_${r}_${p}.blur+orig.HEAD                          \
-censor ${wd}/motion_${p}_bold_${r}_censor.1D                                       \
-ortvec ${wd}/bandpass_rall_${p}_bold_${r}.1D bandpass                                      \
-polort 5                                                             \
-num_stimts 12                                                         \
-stim_file 1 ${wd}/${p}_bold_${r}_motion_demean.1D'[0]' -stim_base 1 -stim_label 1 roll_01  \
-stim_file 2 ${wd}/${p}_bold_${r}_motion_demean.1D'[1]' -stim_base 2 -stim_label 2 pitch_01 \
-stim_file 3 ${wd}/${p}_bold_${r}_motion_demean.1D'[2]' -stim_base 3 -stim_label 3 yaw_01   \
-stim_file 4 ${wd}/${p}_bold_${r}_motion_demean.1D'[3]' -stim_base 4 -stim_label 4 dS_01    \
-stim_file 5 ${wd}/${p}_bold_${r}_motion_demean.1D'[4]' -stim_base 5 -stim_label 5 dL_01    \
-stim_file 6 ${wd}/${p}_bold_${r}_motion_demean.1D'[5]' -stim_base 6 -stim_label 6 dP_01    \
-stim_file 7 ${wd}/motion_deriv_${p}_bold_${r}.1D'[0]' -stim_base 7 -stim_label 7 roll_02   \
-stim_file 8 ${wd}/motion_deriv_${p}_bold_${r}.1D'[1]' -stim_base 8 -stim_label 8 pitch_02  \
-stim_file 9 ${wd}/motion_deriv_${p}_bold_${r}.1D'[2]' -stim_base 9 -stim_label 9 yaw_02    \
-stim_file 10 ${wd}/motion_deriv_${p}_bold_${r}.1D'[3]' -stim_base 10 -stim_label 10 dS_02  \
-stim_file 11 ${wd}/motion_deriv_${p}_bold_${r}.1D'[4]' -stim_base 11 -stim_label 11 dL_02  \
-stim_file 12 ${wd}/motion_deriv_${p}_bold_${r}.1D'[5]' -stim_base 12 -stim_label 12 dP_02  \
-fout -tout -x1D ${wd}/${p}_bold_${r}.X.xmat.1D -xjpeg ${wd}/${p}_bold_${r}.X.jpg                                \
-x1D_uncensored ${wd}/${p}_bold_${r}.X.nocensor.xmat.1D                                     \
-fitts ${wd}/fitts.${p}_bold_${r}                                                    \
-errts ${wd}/errts.${p}_bold_${r}                                                  \
-x1D_stop                                                              \
-bucket ${wd}/stats.${p}_bold_${r} \
-jobs ${n_cpu}

# 3dTproject to project out regression matrix
3dTproject -polort 0 -input ${wd}/pb04.bold_${r}_${p}.blur+orig.HEAD                   \
-censor ${wd}/motion_${p}_bold_${r}_censor.1D -cenmode ZERO                  \
-ort ${wd}/${p}_bold_${r}.X.nocensor.xmat.1D -prefix ${wd}/errts.${p}_bold_${r}.tproject

3dcalc -a ${wd}/errts.${p}_bold_${r}.tproject+orig. -expr '(a)' -prefix ${wd}/errts.${p}_bold_${r}.tproject.nii.gz

3dTstat -prefix ${wd}/${p}_bold_${r}.mean.nii.gz ${wd}/bold_${r}_rest-pec_${p}.nii.gz
