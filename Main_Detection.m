clear all;
close all;

F = uigetdir('','Select Input-folder');

saveFlag =1;
if saveFlag
    svCheck = 1;
    while exist([F '_' 'Segementation_' num2str(svCheck)],'dir')
        svCheck = 1 + svCheck;
    end
    
    mkdir([F '_' 'Segementation_' num2str(svCheck)]);
    
    saveDir = [F '_' 'Segementation_' num2str(svCheck)];
end

%Folder = dir(F);
FileList = dir(fullfile(F, '**', '*.avi'));
clear F;
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
IndexImage = []; % index for the movie

for iFile = 1:numel(FileList)
    file = fullfile(FileList(iFile).folder, FileList(iFile).name);
%     Filenames = Folder(k).name
%     Filenames = fullfile(F,Filenames);
    v = VideoReader(Filenames);
    frames = read(v,[1 Inf]);
    Img =squeeze(frames);
    Planes = size(Img,3);
    clear frames;
    clear v;
    
    %% Detect Neurons and Intensity measurement
    
    Correction = 1; % need intensity correction with blood vessel intensity changes
    
    [FinalC,FinalR,MIntens,MBInten] = DetectRGC(Img,Thresh,sigma,minR,maxR,Sense,DistanceThresh1,DistanceThresh2,Correction,saveDir,FileList(iFile).name(1:end-4)); % last parameter represents whether needs correction
    
    
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
        IndexImage=[IndexImage; iFile*ones(size(FinalR),1)];
    end
    
end

FNMIntens = FFMI;
figure, plot(FNMIntens');
title('Plot of subtracted Calcium waves');

savedata = input('Save the original data as excel file? input (yes) or (no) ','s');

if savedata(1) == 'y'
    filename = input('what would be the file name?:    ','s');
    filename = append(filename,".csv");
    csvwrite(filename,FNMIntens);
    disp('Job done! Saved!');
end