function [dwelltime,switchProbs] = f_calc_dwell_time(data)
elementNum = numel(data);
for i=1:elementNum
    temp = data{i};
    stateInfo(i).elementInfo = [temp(diff(temp)~=0);temp(end)];
    if isrow(temp)
        stateInfo(i).countInfo = diff([0, find(diff(temp)), numel(temp)]);
    else
        stateInfo(i).countInfo = diff([0; find(diff(temp)); numel(temp)]);
    end
    stat = unique(temp);
    for j=1:numel(stat)
        dwelltime(i).dwelltime{j} = (stateInfo(i).countInfo((stateInfo(i).elementInfo == j)));
        dwelltime(i).avgdwelltime(j) = mean(dwelltime(i).dwelltime{j});
    end
    dwelltime(i).avgdwelltime_all = mean(dwelltime(i).avgdwelltime);
    dwelltime(i).switchFreq = sum(diff(temp)~=0)/(2*502);
    
    switchProb_temp = stateInfo(i).elementInfo;
    switchProb = zeros(numel(stat),numel(stat));
    for k=1:size(switchProb_temp,1)-1
        switchProb(switchProb_temp(k),switchProb_temp(k+1)) = switchProb(switchProb_temp(k),switchProb_temp(k+1)) + 1;
    end
    switchProbs{i} = switchProb./sum(switchProb,2);
end
end

