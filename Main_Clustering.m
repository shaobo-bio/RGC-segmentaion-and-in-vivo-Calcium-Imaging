close all;
clear all;


[file,path]  = uigetfile('*.csv','where is the csv file of the traces?');
FItenAll = csvread(fullfile(path,file));
FItenAll = FItenAll-1;

cutoffs = input('Do you want to filter out the traces with less than 0.10 variation?  (yes) or (no):  ','s');

if cutoffs(1) == 'y'
    disp('Removing the cells with less than 10% variation')
    varia=[min(FItenAll,[],2),max(FItenAll,[],2)];
    %varia = max(varia,[],2);
    NoresList = find(varia(:,1) > -0.10 & varia(:,2) < 0.10);
    Noresponse = FItenAll(NoresList,:);
    disp(['Removed  ' num2str(size(NoresList,1)) '  Traces']);
    % Noresponse correspond to the filtered traces (considered as No responding cells)
    FItenAll(NoresList,:)=[];
    clear varia;
end

clear cutoffs;

% savedata = input('Save the data as excel file? input (yes) or (no) ','s');
%
% if savedata(1) == 'y'
%     filename = input('what would be the file name?:    ','s');
%     filename = append(filename,".csv");
%     csvwrite(filename,ItenAll);
%     disp('Job done! Save the result as result.xlsx');
% end

%% Code for clustering
%%% Hierachy Tree Plot
%
% Z = linkage(FItenAll,'ward');
% figure,dendrogram(Z,0);
% title('Hierachy Tree');

%Cluster = input('How many Groups you want to group the cells into for separated analysis?   ');
Cluster =2;
% [B SD] = spca(FNMIntens, [], 20, inf, -[250 125 100], 3000, true);
% PCAInten=B(:,1:11)'* FNMIntens';
%IdxCLU = cluster(Z,'maxclust',2)
disp('Separating the traces into two groups: On-Activation/On-suppression....')

ONOFF=sum(FItenAll(:,45:138),2);
IdxCLU = ONOFF>0;
IdxCLU = IdxCLU+1;

figure

for i=1:Cluster
    subplot(ceil(Cluster/2),2,i);
    plot(FItenAll(IdxCLU==i,:)','g');
    hold on
    MeanCInten(:,i)= mean(FItenAll(IdxCLU(:,1)==i,:));
    plot(mean(FItenAll(IdxCLU==i,:)),'r');
    hold off;
    title(['Group ' num2str(i) ' with ' num2str(size(FItenAll(IdxCLU==i,:),1)) ' Traces']);
end
sgtitle('Waveforms of clustered Calcium waves');

disp('Performing tSNE analysis----');
Y = tsne(FItenAll);
figure, gscatter(Y(:,1),Y(:,2),IdxCLU(:,1),[],'.',12);
title('tSNE plot of the identified clusters');


%% Clustering Analysis for separated groups

SeparateFlag = input('Are you sure to perform clustering analysis separately?(yes) Or (no)   ','s');
%SeparateFlag='yes';
idx=[];
MeanInten = [];

if SeparateFlag(1) == 'y'
    ItenAllGroup = cell(Cluster,1);
    idxSub = cell(Cluster,1);
    MeanIntenSub = cell(Cluster,1);
    RItenAll=[];
    for i = 1:2
        k=35;
        ItenAllGroup{i} = FItenAll(IdxCLU(:,1)==i,:);
        GroupName = ['Group' num2str(i)];
        [idxSub{i},MeanIntenSub{i},~,~] = ClusteringAnalysis(ItenAllGroup{i},GroupName,k,-1,1);
        if i == 1
            idx=idxSub{i};
        else
            idx=[idx;max(idx)+idxSub{i}];
        end
        MeanInten=[MeanInten,MeanIntenSub{i}];
        RItenAll = [RItenAll; ItenAllGroup{i}];
        pause(2);
    end
    
    disp('Combining two taces for further regrouping.');
    GroupName='All Traces';
    FItenAll = RItenAll;
    clear RItenAll;
    [MeanInten,outperm,RegroupID]=ClusterHiearchyTree(max(idx),idx,FItenAll,GroupName,-1,1);
    
    
else
    k=8;
    GroupName='All Traces';
    [idx,MeanInten,outperm,RegroupID] = ClusteringAnalysis(FItenAll,GroupName,k,-1,1);

end

saveFlag = input('Do you want to save the traces and group ID? (yes) or (no):  ');
if saveFlag(1) == 'y'
    csvwrite('Reordered-Trace-File.csv',FItenAll);
    csvwrite('Cluster-IDs-Two-Group-Analysis.csv',idx);
    disp('Results saved as csv files!');
end

figure
for i=1:max(idx)
    subplot(ceil(max(idx)/4),4,i);
    colormap(redbluecmap);
    imagesc(FItenAll(idx==outperm(end+1-i),:));
    xticks([1,45,92,139,186,233,280]);
    xticklabels({'-10','0','10','20','30','40','50'});
    caxis([-1,1]);
    title(['Cluster ' num2str(outperm(end+1-i))]);
end

%% Re-grouping
groupflag = input('Do you want to regroup based on the hierarchy tree cutoff?(yes) or (no)   ','s');

ReID = idx;

if any(groupflag == 'y')
    GroupName='Regrouping';
    for i=1:max(idx)
        ReID(idx==i) = RegroupID(i);
    end
    
else
    Flag = 1;
    while Flag
        CombineC = input('Which clusters do you want to combine? (input like [1,3,4] to end type in [])  ');
        
        if isempty(CombineC)
            Flag = 0;
        else
            MinIndex = min(CombineC);
            for i = 1:numel(CombineC)
                ReID(idx==CombineC(i)) = MinIndex;
            end
        end
    end
    
    UID = unique(ReID);
    for j=1:size(UID)
        ReID(ReID==UID(j)) = j;
    end
    
  
end

[MeanInten,outperm,~]=ClusterHiearchyTree(max(RegroupID),ReID,FItenAll,GroupName,0,2);

figure
for i=1:max(ReID)
    subplot(ceil(max(ReID)/4),4,i);
    colormap(redbluecmap);
    imagesc(FItenAll(ReID==outperm(end+1-i),:));
    xticks([1,45,92,139,186,233,280]);
    xticklabels({'-10','0','10','20','30','40','50'});
    caxis([0,2]);
    title(['Cluster ' num2str(outperm(end+1-i))]);
end

csvwrite('Cluster-IDs-Regrouping-Group.csv',ReID);
%% Measure the properties of each group

viewPlots = 1;
saveFlag = 1;
approxStart = 10;
approxEnd = 30;
tvals = linspace(0,60,283);
tvals = tvals(4:279);

kSigs = input('Enter signals to keep (use brackets): ');

sigsToMeasure = zeros(1,size(FItenAll,2));
rawSigsToMeasure = zeros(1,size(FItenAll,2));
cellIDsToMeasure = {};

for kk = kSigs
    tmpSigs = FItenAll(ReID(:,1)==kk,:);
    sigsToMeasure(end+1:end+size(tmpSigs,1),:) = tmpSigs;
    tmpSigsRaw = FItenAll(ReID(:,1)==kk,:);
    rawSigsToMeasure(end+1:end+size(tmpSigs,1),:) = tmpSigsRaw;
    tmpIDs = num2cell(find(ReID(:,1)==kk));
    cellIDsToMeasure(end+1:end+length(tmpIDs)) = tmpIDs;
end
sigsToMeasure(1,:) = [];
rawSigsToMeasure(1,:) = [];

if saveFlag
    svCheck = 1;
    while exist([path file(1:end-4) '_' 'Analysis_' num2str(svCheck)],'dir')
        svCheck = 1 + svCheck;
    end
    
    mkdir([path file(1:end-4) '_' 'Analysis_' num2str(svCheck)]);
    
    saveDir = [path file(1:end-4) '_' 'Analysis_' num2str(svCheck)];
end

for kk = kSigs
    figure
    plot(FItenAll(ReID==kk,:)','g');
    xticks([1,45,92,139,186,233,280]);
    xticklabels({'-10','0','10','20','30','40','50'});
    %axis([0,283,Ymin,Ymax]);
    hold on
    plot(MeanInten(:,kk),'r');
    title(['Cluster' num2str(kk)])
    hold off;
    
    cellType = input('Please Specify Cell Type: ','s');
    
    if strcmp(cellType,'OFF-TRANS')
        off_trans_Analysis(sigsToMeasure,rawSigsToMeasure,tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
    elseif strcmp(cellType, 'ON-SUST1')
        on_sust1_Analysis(sigsToMeasure,rawSigsToMeasure,tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
    elseif strcmp(cellType, 'ON-SUST2')
        on_sust2_Analysis(sigsToMeasure,rawSigsToMeasure,tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
    elseif strcmp(cellType, 'ON-SUST3')
        on_sust3_Analysis(sigsToMeasure,rawSigsToMeasure,tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
    elseif strcmp(cellType, 'ON-OFF')
        on_off_Analysis(sigsToMeasure,rawSigsToMeasure,tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
    elseif strcmp(cellType, 'ON-TRANS')
        on_trans_Analysis(sigsToMeasure,rawSigsToMeasure,tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
    elseif strcmp(cellType, 'SUPP2')
        supp2_Analysis(sigsToMeasure,rawSigsToMeasure,tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
    elseif strcmp(cellType, 'OFF-SUST')
        off_sust_Analysis(sigsToMeasure,rawSigsToMeasure,tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
    elseif strcmp(cellType, 'SUPP1')
        supp1_Analysis(sigsToMeasure,rawSigsToMeasure,tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
    else
        disp('Cell Type Not Recognized');
    end

end