function [TSmat] = f_extractTSbyAtlas(nii,atlas,covFile,atlasRegion)
covData = importdata(covFile);
[DataNii, VoxelSizeNii, ~, ~] = y_ReadAll(nii);
[DataAtlas, VoxelSizeAtlas, ~, ~] = y_ReadAll(atlas);
if unique(VoxelSizeNii) ~= unique(VoxelSizeAtlas)
    cmd = ['!3dresample -input ' atlas ' -prefix ' fullfile(pwd,'temp.nii.gz -dxyz 0.5 0.5 0.5')];
    eval(cmd)
    [DataAtlas, VoxelSizeAtlas, ~, ~] = y_ReadAll(fullfile(pwd,'temp.nii.gz'));
    !rm temp.nii.gz
end
%idReg = unique(DataAtlas(:));
idReg = atlasRegion;
k=1;
clc
for i=1:size(idReg,2)
    for j=1:size(DataNii,4)
        tempTS_package = (DataNii(:,:,:,j).*(DataAtlas==atlasRegion(i)));
        TS_singleRegion(:,j) = tempTS_package(find(DataAtlas==atlasRegion(i)));
    end
    avgTS_singleReg{k,1} = sum(TS_singleRegion)./size(TS_singleRegion,1);
    [b,r,SSE,SSR] = gretna_regress_ss(avgTS_singleReg{k,1}',[ones(size(covData,1),1),covData]);
    avgTS_singleReg_rg{k,1} = r';
    k=k+1;
    disp(sprintf('region#%02d',i))
    clear TS_singleRegion
end
TSmat = cell2mat(avgTS_singleReg_rg);
disp('ExtractTimeSerials')
