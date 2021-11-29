function [idx,MeanInten,outperm] = ClusteringAnalysis(FNMIntens,GroupName,k,Ymin,Ymax)

%% function for RGC Calcium signal traces clustering analysis

disp(['Now running analysis with ' GroupName]);

[B SD] = spca(FNMIntens, [], 20, inf, -[250 125 100], 3000, true);

figure,imagesc(B');

PCACom = input('How many Principle Components to Process?   ');

% HEATMAP of the PCA components
figure,imagesc(B(:,1:PCACom)');
colorbar;
title(['Heatmap of Sparse PCA components - ' GroupName]);

% Convert the orignal signal profiles into PCA projected profiles
PCAInten=B(:,1:PCACom)'* FNMIntens';

%% options for mdwtcluster
%
% %% plot subcluster signals
% Cluster = input('How many clusters you want to group the cells into?');
% s = mdwtcluster(PCAInten','maxclust',Cluster);
% IdxCLU = s.IdxCLU;
% figure
% for i=1:Cluster
%     subplot(ceil(Cluster/5),5,i);
%     plot(FNMIntens(IdxCLU(:,1)==i,:)','g');
%     hold on
%     MeanCInten(:,i)= mean(FNMIntens(IdxCLU(:,1)==i,:));
%     plot(mean(FNMIntens(IdxCLU(:,1)==i,:)),'r');
%     hold off;
%     title(['Cluster' num2str(i)]);
% end
% sgtitle('Waveforms of clustered Calcium waves');

%% Plot Inten in PCA dimension
% figure, hold on;
% CM = jet(101);
% for i= 1:10
%     plot3(PCAInten(1,IdxCLU(:,1)==i),PCAInten(2,IdxCLU(:,1)==i),PCAInten(3,IdxCLU(:,1)==i),'color',CM(110-10*i,:),'marker','.', 'MarkerSize',12,'LineStyle', 'none' )
%     text(-10,10-0.5*i,['Cluster' num2str(i)],'Color',CM(110-10*i,:));
% end
% hold off;
% title('Plot of PCA-dimension reducted data');

%% tsne plot of the clustered signals


%% Automatic Gaussian Mixture Model Clustering
% automatic caculate the cluster numbers with the lowest BIC index

disp('Now running Clustering analysis and BIC indexing');
BIC=[];
GMModels=cell(40,1);
options = statset('MaxIter',2000);

for k = 1:k
    k
    GMModels{k} = fitgmdist(PCAInten',k,'Options',options,'Replicates',20,'CovarianceType','diagonal');
    BIC(k)=GMModels{k}.BIC;
end

figure,plot(BIC);
title(['BIC index and Cluster Number  ' GroupName]);
[minBIC,index]=min(BIC);
hold on;
plot(index,BIC(index),'ro');
hold off;


% set range for plotting
if BIC(index)>0
    PY=1.05;
else
    PY=0.95;
end
text(index*0.8,BIC(index)*PY,['Index:' num2str(index) num2str(BIC(index))]);

disp(['Based on the BIC Index, ' num2str(index) ' Clusters were found!']);
%GMModels{index} = fitgmdist(PCAInten',index,'Options',options,'Replicates',20,'CovarianceType','diagonal');
IndexFlag= input(['Proceed with ' num2str(index) ' Clusters? (yes) o r(no)  '],'s');
if IndexFlag(1) == 'n'
    index = input('How Many Clusters? You want to cluster the traces into?  '); 
end

idx = cluster(GMModels{index},PCAInten');

%% Plot the clustering results

[MeanInten,outperm]=ClusterHiearchyTree(index,idx,FNMIntens,GroupName,Ymin,Ymax);

end