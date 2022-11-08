restoredefaultpath()
load TSmat

addpath(genpath(fullfile(pwd,'LEiDA')))
%%
% Get the BOLD data Subjects names
N_sub=length(TSmat);

% Set parameters and preallocate variables
N=54; % Number of AAL brain areas considered
Tmax=502; % TR=2s, so total time is 1004s
Phases=zeros(N,Tmax);
Leading_Eig=zeros(N_sub,Tmax,N);
Var_Eig=zeros(N_sub,Tmax);
FCD_eig =zeros(N_sub,Tmax-2,Tmax-2);

for s=1:N_sub
    % Load BOLD data and save it as "signaldata'
    disp(['  Subject ' num2str(s) ' from ' num2str(N_sub)])
    signaldata=TSmat{s};
    
    % Get the Phase of BOLD data
    for seed=1:N
        timeseriedata=demean(detrend(signaldata(:,seed)));
        Phases(seed,:) = angle(hilbert(timeseriedata));
    end
    
    % Get the iFC leading eigenvector at each time point 
    iFC=zeros(N);
    
    for t=1:Tmax
        for n=1:N
            for p=1:N
                iFC(n,p)=cos(adif(Phases(n,t),Phases(p,t)));
            end
        end
        [eVec,eigVal]=eig(iFC);
        eVal=diag(eigVal);
        [val1, i_vec_1] = max(eVal);
        Leading_Eig(s,t,:)=eVec(:,i_vec_1);
        % Calculate the variance explained by the leading eigenvector
        Var_Eig(s,t)=val1/sum(eVal);
    end
    
    % Compute the FCD for each subject
    for t=1:Tmax
        eig1=squeeze(Leading_Eig(s,t,:));
        for t2=1:Tmax
            eig2=squeeze(Leading_Eig(s,t2,:));
            FCD_eig(s,t,t2)=dot(eig1,eig2)/(norm(eig1)*norm(eig2));
            %FCD_eig(s,t,t2)=corr(eig1,eig2);
        end
    end
end
plot(signaldata)
ylabel('BOLD Signal')
xlabel('Time (TR)')
xlim([0,502])

save('LEiDA_data','Leading_Eig','FCD_eig','Var_Eig')

iFD_temp = squeeze(FCD_eig(8,:,:));
imagesc(iFD_temp)
colormap(hot)
axis square
caxis([-1 1])
%%
% Generate vector X concatenating all eigenvectors from all subjects and
% time points, where rows are observations and collumns are variables.
clear,clc
load LEiDA_data Leading_Eig
load TSmat
N_sub=length(TSmat);

X=[];
for s=1:N_sub
    X=cat(1,X,squeeze(Leading_Eig(s,:,:))); 
    corrLeadingEig(:,:,s) = corr(squeeze(Leading_Eig(s,:,:))');
end
clear Leading_Eig

%Kmeans clustering
maxk=10;
opt= statset('UseParallel',0); %,'UseSubstreams',1);
% The options may vary according to the Matlab version
Kmeans_results=cell(1,maxk);
rng(35) % for reproducibility
for k=2:maxk  
    disp(['Calculating for ' num2str(k) 'clusters'])
    [IDX, C, SUMD, D]=kmeans(X,k,'Distance','cityblock','Replicates',20,'Display','final','Options',opt); %,'Options',opt);  
    %[IDX, C, SUMD, D]=kmeans(X,k,'Distance','cityblock','Replicates',20,'Display','final'); %,'Options',opt);   
    Kmeans_results{k}.IDX=IDX;
    Kmeans_results{k}.C=C; 
    Kmeans_results{k}.SUMD=SUMD; 
    Kmeans_results{k}.D=D; 
end

%Evaluate Clustering performance
distM_fcd=squareform(pdist(X,'cityblock'));
dunn_score=zeros(maxk,1);
for j=2:maxk
    dunn_score(j)=dunns(j,distM_fcd,Kmeans_results{j}.IDX);
    disp(['Performance for ' num2str(j) 'clusters'])
end
plot(dunn_score)
ylabel('Dunn''s index')
xlabel('k')

[~,ind_max]=max(dunn_score);
disp(['Best clustering solution: ' num2str(ind_max) ' clusters']);

Clusters= Kmeans_results{ind_max};

save('LEiDA_Clusters','Clusters')

%%
load LEiDA_Clusters Clusters
[N_Cl, N_ba]=size(Clusters.C);
h=hist(Clusters.IDX,N_Cl);
[y, ind]=sort(h,'descend');
V=Clusters.C(ind,:);
VVT_mean=zeros(N_ba);

% To reorder matrix plots
Order=[1:2:N_ba N_ba:-2:2];

figure
colormap(jet)    
% Pannel A
% Plot the centroids outerproduct

for c=1:N_Cl
%     subplot(2,N_Cl+1,c)
%     plot_nodes_in_cortex(V(c,:))
%     title(['#' num2str(c)])
%     axis off
    subplot(1,N_Cl+1,c)
    VVT=V(c,:)'*V(c,:);   
    imagesc(VVT(Order,Order))   
    caxis([-0.02 0.02])
    axis square
    title('Outer product') 
    ylabel('Brain area #')
    xlabel('Brain area #') 
    VVT_mean=VVT_mean+squeeze(VVT)*y(c);
end

VVT_mean=VVT_mean/sum(y);

% Pannel B
% Plot the probabilities of each cluster
% The Colormap is adjusted for 5 clusters 
%mymap=[1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1];

x=y/sum(y);
%subplot(2,N_Cl,N_Cl+1)
hold on
for c=1:N_Cl
    %bar(c,x(c),'Facecolor',mymap(c,:))
    bar(c,x(c))
end
set(gca,'xtick',1:N_Cl)
box off
xlabel('State #')
ylabel('Probability')
xlim([0 8])

% Pannel C plot the weighted sum of outer products
%subplot(2,N_Cl+1,N_Cl+3)
imagesc(VVT_mean(Order,Order))
title({'Weighted sum of VV^T'})
ylabel('Brain area')
xlabel('Brain area')
colorbar
axis square

load TSmat
for i=1:numel(TSmat)
    temp(:,:,i) = corr(TSmat{i});
end
Static_FC_mean = nanmean(temp,3);
I_sup_diag=find(triu(ones(N_ba),3));
[cc p]=corrcoef(Static_FC_mean(I_sup_diag),VVT_mean(I_sup_diag))
%subplot(2,N_Cl+1,N_Cl+4)
scatter(Static_FC_mean(I_sup_diag),VVT_mean(I_sup_diag))
%imagesc(atanh(Static_FC_mean))
imagesc(Static_FC_mean)
%caxis([-0.5 0.8])
axis square
ylabel('Brain area #')
xlabel('Brain area #')
h=colorbar
h.Label.String = 'Pearson r';

%% extract state sequence
TRcount = 502;
for i=1:size(Clusters.IDX,1)/TRcount
    ss{1,i} = Clusters.IDX((i-1)*TRcount+1:i*TRcount,:);
end
% lifetimes
[results.dwelltime switchProbs] = f_calc_dwell_time(ss);
tblDwell = struct2table(results.dwelltime);
ageSplit = [1,0,0,0,0,1,1,1,1]';
sexMarmoset = [1,1,0,1,0,1,0,1,0]';
age = [71,24,41,40,22,63,75,78,77]';
statData = [tblDwell(:,3),array2table(sexMarmoset),array2table(age),array2table(tblDwell.avgdwelltime),array2table(ageSplit),tblDwell(:,end)];
statData.Properties.VariableNames = {'avgLifetime','sex','age','s1','s2','s3','s4','s5','s6','s7','ageSplit','switchFreq'};
[p,h,stats] = ranksum(statData.switchFreq(statData.ageSplit == 1),statData.switchFreq(statData.ageSplit == 0));
grpstats(statData,'ageSplit')
corr(statData.age,statData.switchFreq)
scatter(statData.age,statData.switchFreq)

writetable(statData,'statData.csv');
%% plot on surface
addpath('/home/jinbo/Project/JojiRotation/gifti');
% load V vector
load LEiDA_Clusters Clusters
[N_Cl, N_ba]=size(Clusters.C);
h=hist(Clusters.IDX,N_Cl);
[y, ind]=sort(h,'descend');
V=Clusters.C(ind,:);
for i=1:size(V,1)
    singleState = V(i,:)';
    f_plotSurface(singleState,sprintf('state#%02d',i));
end
%% stat analysis