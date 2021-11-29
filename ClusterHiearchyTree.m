function [MeanInten,outperm]=ClusterHiearchyTree(index,idx,FNMIntens,GroupName,Ymin,Ymax)

MeanInten=[];
figure
for i=1:index
    subplot(ceil(index/4),4,i);
    plot(FNMIntens(idx==i,:)','g');
    axis([0,283,Ymin,Ymax]);
    hold on
    MeanInten(:,i)= mean(FNMIntens(idx==i,:));
    plot(MeanInten(:,i),'r');
    title(['Cluster' num2str(i)])
end

sgtitle(['Waveforms of clustered Calcium waves(guassian mixture model) -- ' GroupName]);

Z = linkage(MeanInten','average','correlation');
D = pdist(MeanInten');
leafOrder = optimalleaforder(Z,D);
figure
[H,T,outperm]=dendrogram(Z,0,'reorder',leafOrder,'ColorThreshold',0.05,'Orientation','left');
title(['Hierchay Plot of the Clusters for Further Grouping --' GroupName]);
% outperm contains the ordered leaves #

figure
imagesc(MeanInten(:,flip(outperm))');
caxis([Ymin Ymax]);
colormap(redbluecmap);
set(gca,'YTick',[1:size(outperm,2)],'YTickLabel',flip(outperm),'TickDir','out','TickLength',[0 0],'fontsize', 14);
title(['Heatmap of the mean intensity changes of each cluster in --' GroupName]);



% cgo = clustergram(MeanInten','Colormap',redbluecmap,'Symmetric',false);
% cgo.Dendrogram = 0.5;
% 
% cgAxes =plot(cgo);
% 
% set(cgAxes, 'Clim', [min(MeanInten,[],'all')*0.9,max(MeanInten,[],'all')*0.7]);
% title(['Clustergram Plot of the Clusters for Further Grouping --' GroupName]);
% disp(['Clustering analysis done with ' GroupName]);



tSNEFLAG=input('Do you want tSNE plot?(yes) or (no)  ','s');
if tSNEFLAG(1) == 'y'
    Y = tsne(FNMIntens);
    figure, gscatter(Y(:,1),Y(:,2),idx,[],'.',12);
    title(['tSNE plot of the identified clusters in ' GroupName]);
end

end