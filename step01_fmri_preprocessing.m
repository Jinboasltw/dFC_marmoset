clear,clc
%% Jinbo Zhang
% 2022-08-18, zhangjinbo@cibr.ac.cn
%% software version
% All code run on a laptop with Ubuntu MATE 22.04.1 LTS x86_64
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AFNI: Version AFNI_22.2.05 'Marcus Aurelius'
% FSL: FSL 6.0.5
% ANTS: 2.4
% Connectome Workbench: Version 1.5.0
% Python: 3.9.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pyton -> run in terminal: conda create --name marmosetENV python=3.9.2
% fsl -> sudo python fslinstaller.py
% ANTS -> unzip and add to path
% Connectome Workbench -> unzip and add to path
% AFNI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% META INFO
wd = which('step01_fmri_preprocessing.m');
[wd_pth,~,~] = fileparts(wd);
addpath(fullfile(wd_pth,'code','lib'));
mList = f_ls_super(fullfile(wd_pth,'marmoset_rsfMRI_data','sub-*'));
callIt = @(n) eval(n);
% ANALYSIS START
%% PREPROCESSING
for iSubj=6:numel(mList)
    [~,mID,~] = fileparts(mList{iSubj});
    temp_dir = fullfile(wd_pth,'temp',mID);
    if ~exist(temp_dir)
        mkdir(temp_dir)
    end
    output_dir = fullfile(wd_pth,'output',mID);
    if ~exist(output_dir)
        mkdir(output_dir)
    end
    %% Key job: phase encoding correction
    % step01: average SEEPI image
    cmd=['!fslmaths ' mList{iSubj} filesep 'SEEPI_up.nii.gz -Tmean ' temp_dir filesep 'SEEPI_up_mean.nii.gz'];
    callIt(cmd);
    cmd=['!fslmaths ' mList{iSubj} filesep 'SEEPI_down.nii.gz -Tmean ' temp_dir filesep 'SEEPI_down_mean.nii.gz'];
    callIt(cmd);
    
    % step02: extract middle volume and registration to it
    runs_ups = f_ls_super(fullfile(wd_pth,'marmoset_rsfMRI_data',mID,'BOLD_up_*.nii.gz'));
    runs_downs = f_ls_super(fullfile(wd_pth,'marmoset_rsfMRI_data',mID,'BOLD_down_*.nii.gz'));
    for irun=1:1%numel(runs_ups)
        %calc middle mean which picture
        targetBOLD_up = runs_ups{irun};
        targetBOLD_down = runs_downs{irun};
        cmd=['!fslinfo ' targetBOLD_up '> temp_' num2str(irun) '.txt'];
        callIt(cmd)
        content=importdata(['temp_' num2str(irun) '.txt']);
        timePointsMiddle = round(content.data(4)/2)-1; % afni sub-brick start at 0
        %extract middle picture
        cmd=['!3dcalc -a ' runs_ups{irun} '[' num2str(timePointsMiddle) '] -expr ''(a)'' -prefix ' [temp_dir filesep 'rest-up_bold_' num2str(irun) '_middle_vol.nii.gz']];
        callIt(cmd)
        
        %concatenate for registration
        cmd=['!fslmerge -t ' [temp_dir filesep 'merged_bold_', num2str(irun),'.nii.gz '] [temp_dir filesep 'rest-up_bold_' num2str(irun) '_middle_vol.nii.gz '] [temp_dir filesep 'SEEPI_up_mean.nii.gz '] temp_dir filesep 'SEEPI_down_mean.nii.gz'];
        callIt(cmd)
        
        %register "up" and "down" runs to chosen registration volume
        cmd=['!3dvolreg -prefix ' temp_dir filesep 'reg_up_down_to_middle_bold_' num2str(irun) '.nii.gz ' [temp_dir filesep 'merged_bold_', num2str(irun),'.nii.gz']];
        callIt(cmd)
        
        %extract the registred volumes from "within-volume" registration
        cmd=['!3dcalc -a ' fullfile(temp_dir,['reg_up_down_to_middle_bold_' num2str(irun) '.nii.gz']) '[1]' ' -expr ' '''(a)''' ' -prefix ' fullfile(temp_dir,['SEEPI_up_mean_reg_bold_' num2str(irun) '.nii.gz'])];
        callIt(cmd)
        cmd=['!3dcalc -a ' fullfile(temp_dir,['reg_up_down_to_middle_bold_' num2str(irun) '.nii.gz']) '[2]' ' -expr ' '''(a)''' ' -prefix ' fullfile(temp_dir,['SEEPI_down_mean_reg_bold_' num2str(irun) '.nii.gz'])];
        callIt(cmd)
        
        %merge up and down volumes
        cmd=['!fslmerge -t ' fullfile(temp_dir,['SEEPI_merged_reg_bold_',num2str(irun),'.nii.gz ']) fullfile(temp_dir,['SEEPI_up_mean_reg_bold_' num2str(irun) '.nii.gz ']) fullfile(temp_dir,['SEEPI_down_mean_reg_bold_' num2str(irun) '.nii.gz'])];
        callIt(cmd)
        
        %run topup
        cmd = ['!topup --imain=' fullfile(temp_dir,['SEEPI_merged_reg_bold_',num2str(irun),'.nii.gz']) ' --datain=' fullfile(wd_pth,'distortion_correction_parameters','acqparams_1_0_0.txt')  ' --config=' fullfile(wd_pth,'distortion_correction_parameters','updown.cnf') ' --out=' fullfile(output_dir,['bold_' num2str(irun) '_topup'])];
        callIt(cmd)
        
        %volume removal for steady-state stabilization, default is the first 10 volumes
        cmd = ['!3dTcat -prefix ' [temp_dir,filesep,'pb00.up_',num2str(irun),'.tcat ',targetBOLD_up,'''[10..$]''']];
        callIt(cmd)
        cmd = ['!3dTcat -prefix ' [temp_dir,filesep,'pb00.down_',num2str(irun),'.tcat ',targetBOLD_down,'''[10..$]''']];
        callIt(cmd)
        
        %data check: compute outlier fraction for each volume
        cmd = ['!touch ' temp_dir filesep 'out.pre_ss_warn.txt'];
        callIt(cmd)
        cmd = ['!3dToutcount -automask -fraction -polort 5 -legendre ' [temp_dir,filesep,'pb00.up_',num2str(irun),'.tcat+orig > '] fullfile(temp_dir,['outcount.up_bold_',num2str(irun),'.1D'])];
        callIt(cmd)
        cmd = ['!3dToutcount -automask -fraction -polort 5 -legendre ' [temp_dir,filesep,'pb00.down_',num2str(irun),'.tcat+orig > '] fullfile(temp_dir,['outcount.down_bold_',num2str(irun),'.1D'])];
        callIt(cmd)
        
        %catenate outlier counts into a single time series
        cmd = ['!cat ',fullfile(temp_dir,['outcount.up_bold_',num2str(irun),'.1D']),'>',fullfile(temp_dir,['outcount.up_bold_',num2str(irun),'_all.1D'])];
        callIt(cmd);
        cmd = ['!cat ',fullfile(temp_dir,['outcount.down_bold_',num2str(irun),'.1D']),'>',fullfile(temp_dir,['outcount.down_bold_',num2str(irun),'_all.1D'])];
        callIt(cmd);
        
        %apply 3dDespike to each run
        cmd = ['!3dDespike -NEW -nomask -prefix ' [temp_dir,filesep,'pb01.up_',num2str(irun),'.despike '] [temp_dir,filesep,'pb00.up_',num2str(irun),'.tcat+orig']];
        callIt(cmd)
        cmd = ['!3dDespike -NEW -nomask -prefix ' [temp_dir,filesep,'pb01.down_',num2str(irun),'.despike '] [temp_dir,filesep,'pb00.down_',num2str(irun),'.tcat+orig']];
        callIt(cmd)
        
        %time shift data so all slice timing is the same
        cmd = ['!3dTshift -tzero 0 -quintic -prefix ', fullfile(temp_dir,['pb02.up_bold_',num2str(irun),'.tshift ']),[temp_dir,filesep,'pb01.up_',num2str(irun),'.despike+orig']];
        callIt(cmd);
        cmd = ['!3dTshift -tzero 0 -quintic -prefix ', fullfile(temp_dir,['pb02.down_bold_',num2str(irun),'.tshift ']),[temp_dir,filesep,'pb01.down_',num2str(irun),'.despike+orig']];
        callIt(cmd);
        
        %align each dset to base volume, base volume calculated in pe correction above, timePointsMiddle volume of "UP" before stabilzation volume removal
        cmd = ['!3dvolreg -verbose -zpad 1 -base ' fullfile(temp_dir,'rest-up_bold_1_middle_vol.nii.gz') ' -1Dfile ' fullfile(temp_dir,['dfile.up_bold_' num2str(irun) '.1D']) ' -prefix ' fullfile(temp_dir,['pb03.up_bold_' num2str(irun) '.volreg.nii.gz']) ' -cubic ' fullfile(temp_dir,['pb02.up_bold_',num2str(irun),'.tshift+orig'])];
        callIt(cmd)
        cmd = ['!3dvolreg -verbose -zpad 1 -base ' fullfile(temp_dir,'rest-up_bold_1_middle_vol.nii.gz') ' -1Dfile ' fullfile(temp_dir,['dfile.down_bold_' num2str(irun) '.1D']) ' -prefix ' fullfile(temp_dir,['pb03.down_bold_' num2str(irun) '.volreg.nii.gz']) ' -cubic ' fullfile(temp_dir,['pb02.down_bold_',num2str(irun),'.tshift+orig'])];
        callIt(cmd)
        
        %make a single file of registration params
        cmd = ['!cat ' fullfile(temp_dir,['dfile.up_bold_' num2str(irun) '.1D']) ' > ' fullfile(temp_dir,['dfile.up_bold_' num2str(irun) '_all.1D'])];
        callIt(cmd)
        cmd = ['!cat ' fullfile(temp_dir,['dfile.down_bold_' num2str(irun) '.1D']) ' > ' fullfile(temp_dir,['dfile.down_bold_' num2str(irun) '_all.1D'])];
        callIt(cmd)
        
        %apply PE correction to registered images
        cmd = ['!applytopup --imain=' fullfile(temp_dir,['pb03.up_bold_' num2str(irun) '.volreg.nii.gz ']) '--topup=' fullfile(output_dir,['bold_' num2str(irun) '_topup ']) '--datain=' fullfile(wd_pth,'distortion_correction_parameters','acqparams_1_0_0.txt') ' --method=jac --inindex=1 --out=' fullfile(temp_dir,['bold_' num2str(irun) '_rest-pec_up.nii.gz'])]
        callIt(cmd)
        cmd = ['!applytopup --imain=' fullfile(temp_dir,['pb03.down_bold_' num2str(irun) '.volreg.nii.gz ']) '--topup=' fullfile(output_dir,['bold_' num2str(irun) '_topup ']) '--datain=' fullfile(wd_pth,'distortion_correction_parameters','acqparams_1_0_0.txt') ' --method=jac --inindex=2 --out=' fullfile(temp_dir,['bold_' num2str(irun) '_rest-pec_down.nii.gz'])]
        callIt(cmd)
        
        %estimate motion
        cmd = ['!3dvolreg -verbose -zpad 1 -based ' num2str(timePointsMiddle) ' -1Dfile ' fullfile(temp_dir,['dfile.up_bold_' num2str(irun) '_est.1D']) ' -prefix ' fullfile(temp_dir,['pb03.up_bold_' num2str(irun) '.volreg.est.nii.gz ' ' -cubic ' fullfile(temp_dir,['pb02.up_bold_',num2str(irun),'.tshift+orig'])])];
        callIt(cmd)
        
        cmd = ['!3dvolreg -verbose -zpad 1 -based ' num2str(timePointsMiddle) ' -1Dfile ' fullfile(temp_dir,['dfile.down_bold_' num2str(irun) '_est.1D']) ' -prefix ' fullfile(temp_dir,['pb03.down_bold_' num2str(irun) '.volreg.est.nii.gz ' ' -cubic ' fullfile(temp_dir,['pb02.down_bold_',num2str(irun),'.tshift+orig'])])];
        callIt(cmd)
        
        %blur each volume of each run
        cmd = ['!3dmerge -1blur_fwhm 1.5 -doall -prefix ' fullfile(temp_dir,['pb04.bold_',num2str(irun),'_up.blur ']) fullfile(temp_dir,['bold_' num2str(irun) '_rest-pec_up.nii.gz'])];
        callIt(cmd)
        cmd = ['!3dmerge -1blur_fwhm 1.5 -doall -prefix ' fullfile(temp_dir,['pb04.bold_',num2str(irun),'_down.blur ']) fullfile(temp_dir,['bold_' num2str(irun) '_rest-pec_down.nii.gz'])];
        callIt(cmd)
        
        %motion regressin
        cmd = ['!./motionHandle.sh ' temp_dir ' up ', num2str(irun)];callIt(cmd);
        cmd = ['!./motionHandle.sh ' temp_dir ' down ', num2str(irun)];callIt(cmd);
    end
end
%% From Nativ to Standard Space
%downsampling MBM template https://marmosetbrainmapping.org/atlas.html#v3 to .5mm isotropic
cmd = ['!3dresample -input ' fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1','MBM_v3.0.1','template_T2w_brain.nii.gz ') '-prefix ' fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1','MBM_v3.0.1','template_T2w_brain.5iso.nii.gz -dxyz 0.5 0.5 0.5')];
callIt(cmd)
cmd = ['!3dresample -input ' fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1','MBM_v3.0.1','segmentation_three_types_prob_1_gray.nii.gz ') '-prefix ' fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1','MBM_v3.0.1','gray_prob.5iso.nii.gz -dxyz 0.5 0.5 0.5')];
callIt(cmd)
cmd = ['!3dcalc -a ' fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1','MBM_v3.0.1','gray_prob.5iso.nii.gz ') '-expr ''step(a-0.5)'' -prefix ' fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1','MBM_v3.0.1','gray_template.5iso.nii.gz')];
callIt(cmd)

cmd = ['!3dresample -input ' fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1','MBM_v3.0.1','segmentation_three_types_prob_2_white.nii.gz ') '-prefix ' fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1','MBM_v3.0.1','white_prob.5iso.nii.gz -dxyz 0.5 0.5 0.5')];
callIt(cmd)
cmd = ['!3dcalc -a ' fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1','MBM_v3.0.1','white_prob.5iso.nii.gz ') '-expr ''step(a-0.5)'' -prefix ' fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1','MBM_v3.0.1','white_matter.5iso.nii.gz')];
callIt(cmd)

cmd = ['!3dresample -input ' fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1','MBM_v3.0.1','segmentation_three_types_prob_3_csf.nii.gz ') '-prefix ' fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1','MBM_v3.0.1','csf_prob.5iso.nii.gz -dxyz 0.5 0.5 0.5')];
callIt(cmd)
cmd = ['!3dcalc -a ' fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1','MBM_v3.0.1','csf_prob.5iso.nii.gz ') '-expr ''step(a-0.5)'' -prefix ' fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1','MBM_v3.0.1','csf.5iso.nii.gz')];
callIt(cmd)

for iSubj=5:numel(mList)
    [~,mID,~] = fileparts(mList{iSubj});
    temp_dir = fullfile(wd_pth,'temp',mID);
    if ~exist(temp_dir)
        mkdir(temp_dir)
    end
    output_dir = fullfile(wd_pth,'output',mID);
    if ~exist(output_dir)
        mkdir(output_dir)
    end
    corr_dir = fullfile(wd_pth,'corrpth',mID);
    if ~exist(corr_dir)
        mkdir(corr_dir)
    end
    % if killed by system, clean sys mem first
    % sudo swapoff -a
    % sudo sync && echo 3 | sudo tee /proc/sys/vm/drop_caches***
    % sudo swapon -a
    for irun=1:1%numel(runs_ups)
        cmd=['!./registration_nat2template.sh ',mList{iSubj},' ',mList{iSubj},' ',corr_dir,' ',fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1','MBM_v3.0.1'),' ', output_dir, ' ', temp_dir, ' up ',num2str(irun)];
        callIt(cmd);
        cmd=['!./applyReg_createRegressor.sh ',mList{iSubj},' ',mList{iSubj},' ',corr_dir,' ',fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1','MBM_v3.0.1'),' ', output_dir, ' ', temp_dir, ' up ',num2str(irun)];
        callIt(cmd);
        pause(15)
        cmd=['!./registration_nat2template.sh ',mList{iSubj},' ',mList{iSubj},' ',corr_dir,' ',fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1','MBM_v3.0.1'),' ', output_dir, ' ', temp_dir, ' down ',num2str(irun)];
        callIt(cmd);
        cmd=['!./applyReg_createRegressor.sh ',mList{iSubj},' ',mList{iSubj},' ',corr_dir,' ',fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1','MBM_v3.0.1'),' ', output_dir, ' ', temp_dir, ' down ',num2str(irun)];
        callIt(cmd);
    end
end

%% EXTRACT SIGNAL
addpath('/data/home/jinbo/Project/JojiRotation/code/toolbox/spm12');
addpath(genpath('/data/home/jinbo/Project/JojiRotation/code/toolbox/DPABI'));
addpath(genpath('/data/home/jinbo/Project/JojiRotation/GRETNA'));
atlas = fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1/MBM_v3.0.1','atlas_MBM_cortex_vM.nii.gz');
atlasRegion = [1:54];
for iSubj=1:numel(mList)
    [~,mID,~] = fileparts(mList{iSubj});
    temp_dir = fullfile(wd_pth,'temp',mID);
    if ~exist(temp_dir)
        mkdir(temp_dir)
    end
    output_dir = fullfile(wd_pth,'output',mID);
    if ~exist(output_dir)
        mkdir(output_dir)
    end
    corr_dir = fullfile(wd_pth,'corrpth',mID);
    if ~exist(corr_dir)
        mkdir(corr_dir)
    end
    % if killed by system, clean sys mem first
    % sudo swapoff -a
    % sudo sync && echo 3 | sudo tee /proc/sys/vm/drop_caches***
    % sudo swapon -a
    for irun=1:1%numel(runs_ups)
        covFile = fullfile(corr_dir,'down_bold_1_nui_regressors.1D');
        niiFile = fullfile(corr_dir,'down_bold_1_to_template_errts.5iso_masked_gm.nii.gz');
        TSmat{iSubj,1} = f_extractTSbyAtlas(niiFile,atlas,covFile,atlasRegion)';
        %         cmd=['!./extractTimeserial.sh ','up ',num2str(irun),' ',temp_dir,' ',corr_dir,' ',fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1','MBM_v3.0.1','gray_template.5iso.nii.gz')];
        %         callIt(cmd);
        %         pause(30)
        %         cmd=['!./extractTimeserial.sh ','down ',num2str(irun),' ',temp_dir,' ',corr_dir,' ',fullfile(wd_pth,'Marmoset_Brain_Mappping_v3.0.1','MBM_v3.0.1','gray_template.5iso.nii.gz')];
        %         callIt(cmd);
    end
end
save TSmat TSmat