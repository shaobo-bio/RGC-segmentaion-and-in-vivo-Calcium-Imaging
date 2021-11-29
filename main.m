clear all;
close all;


Folder = dir('/Users/CaVideos');

minR = 5; %5
maxR = 16; %20 for GCamp

% the distance threshold to merge circles
DistanceThresh1 = 12; % 12 for GCamp
DistanceThresh2 = 10; % 10 for GCamp

% sensitivity of the circle detection
Sense = 0.87;
% Sense = 0.92;

% sigma for the gaussian filter blur to extract background intensity
sigma= 20; % 20 for GCamp

% Threshold for subtracting intensity to generate BW images
Thresh =1.25; % 1 for GCamp

FFC=[];
FFR=[];
FFMI=[];

for k = 4:size(Folder)
    k
    Filenames = Folder(k).name
    v = VideoReader(Filenames);
    frames = read(v,[1 Inf]);
    Img =squeeze(frames);
    Planes = size(Img,3);
    clear frames;
    clear v;
    
    %% Detect Neurons and Intensity measurement
    
    Correction = 1; % need intensity correction with blood vessel intensity changes
    
    [FinalC,FinalR,MIntens,MBInten] = DetectRGC(Img,Thresh,sigma,minR,maxR,Sense,DistanceThresh1,DistanceThresh2,Correction); % last parameter represents whether needs correction
    
    %MIntens: Mean Intensity of individual RGCs
    %MBIten: Mean Blood vessel Intensity;
    
    if isempty(MIntens)
        NMIntens = [];
    else
        NMIntens = MIntens./mean(MIntens(:,1:45),2);
        
        if Correction == 1 & isempty(MBInten) == 0
            NMBInten = MBInten./mean(MBInten(:,1:45),2);
            NMIntens = NMIntens./NMBInten;
        end
    end
    
    KeepFlag = input('Do you want to keep the data from this movie?(yes) or (no):  ','s');
    
    if KeepFlag(1) == 'y'
        FFC=[FFC;FinalC];
        FFR=[FFR;FinalR];
        FFMI=[FFMI;NMIntens];
    end
    
end


FNMIntens = FFMI;
figure, plot(FNMIntens');
title('Plot of subtracted Calcium waves');

% lowpass filtering
FilterFlag=input('Do you want to use lowpass filter to smooth the data?(yes) or (no) ','s');

if any(FilterFlag ==  'y')
    disp('Now processding the trace with Lowpass to remove noises... (May take some time ...)')
    for i=1:size(FNMIntens,1)
        FS=lowpass(FNMIntens(i,:),0.3);
        FS(1:3)=[];
        FS(end-3:end)=[];
        FilteredFNMItens(i,:) = FS;
    end
    
    FNMIntens=FilteredFNMItens;
end

%% Remove cells with variation less than 0.15
% cutoffs = input('Do you want to filter out the traces with less than 0.15 variation?  (yes) or (no):  ','s');
%
% if cutoffs(1) == 'y'
%     disp('Removing the cells with less than 10% variation')
%     varia=[abs(1-min(FNMIntens,[],2)),abs(1-max(FNMIntens,[],2))];
%     varia = max(varia,[],2);
%     FNMIntens(find(varia<0.10),:)=[];
% end


savedata = input('Save the data as excel file? input (yes) or (no) ','s');

if savedata(1) == 'y'
    filename = input('what would be the file name?:    ','s');
    filename = append(filename,".csv");
    csvwrite(filename,FNMIntens);
    disp('Job done! Save the result as result.xlsx');
end

%% Code for clustering
%%% Hierachy Tree Plot
Z = linkage(FNMIntens,'ward');
figure,dendrogram(Z,0);
title('Hierachy Tree');

Cluster = input('How many Groups you want to group the cells into for separated analysis?   ');

% [B SD] = spca(FNMIntens, [], 20, inf, -[250 125 100], 3000, true);
% PCAInten=B(:,1:11)'* FNMIntens';


IdxCLU = cluster(Z,'maxclust',2)
%s = mdwtcluster(FNMIntens,'maxclust',Cluster);

ONOFF=sum(FNMIntens(:,48:142),2)-94;
IdxCLU = ONOFF>0;
IdxCLU = IdxCLU+1;

%IdxCLU = s.IdxCLU;
figure

for i=1:Cluster
    subplot(ceil(Cluster/2),2,i);
    plot(FNMIntens(IdxCLU(:,1)==i,:)','g');
    hold on
    MeanCInten(:,i)= mean(FNMIntens(IdxCLU(:,1)==i,:));
    plot(mean(FNMIntens(IdxCLU(:,1)==i,:)),'r');
    hold off;
    title(['Group ' num2str(i) ' with ' num2str(size(FNMIntens(IdxCLU(:,1)==i,:),1)) ' Traces']);
end
sgtitle('Waveforms of clustered Calcium waves');

disp('Performing tSNE analysis----');
Y = tsne(FNMIntens);
figure, gscatter(Y(:,1),Y(:,2),IdxCLU(:,1),[],'.',12);
title('tSNE plot of the identified clusters');

%% Clustering Analysis for separated groups

SeparateFlag = input('Are you sure to perform clustering analysis separately?(yes) Or (no)   ','s');

idx=[];
MeanInten = [];

if SeparateFlag(1) == 'y'
    FNMIntensGroup = cell(Cluster,1);
    idxSub = cell(Cluster,1);
    MeanIntenSub = cell(Cluster,1);
    RFNMIntens=[];
    for i = 1: Cluster
        k=35;
        FNMIntensGroup{i} = FNMIntens(IdxCLU(:,1)==i,:);
        GroupName = ['Group' num2str(i)];
        [idxSub{i},MeanIntenSub{i},~] = ClusteringAnalysis(FNMIntensGroup{i},GroupName,k,0,2);
        if i == 1
            idx=idxSub{i};
        else
            idx=[idx;max(idx)+idxSub{i}];
        end
        MeanInten=[MeanInten,MeanIntenSub{i}];
        RFNMIntens = [RFNMIntens; FNMIntensGroup{i}];
        pause(2);
    end
    
    disp('Combining two taces for further regrouping.');
    GroupName='All Traces';
    [MeanInten,outperm]=ClusterHiearchyTree(max(idx),idx,RFNMIntens(:,10:273),GroupName,0,2);
    FNMIntens = RFNMIntens;
    
else
    k=40;
    GroupName='All Traces';
    [idx,MeanInten,outperm] = ClusteringAnalysis(FNMIntens,GroupName,k,0,2);
end


figure
for i=1:max(idx)
    subplot(ceil(max(idx)/4),4,i);
    imagesc(FNMIntens(idx==outperm(end+1-i),:));
    caxis([0.5,2]);
    title(['Cluster ' num2str(outperm(end+1-i))]);
end
