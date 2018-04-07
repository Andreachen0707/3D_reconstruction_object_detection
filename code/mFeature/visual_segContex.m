clear all;
load(['/n/fs/modelnet/NYUdataSet/NYUdatafeatureNew/' 'feature_' num2str(1) '.mat'],'XYZworldframeTest','SpaceTest','featureTest','Rorg','Rtilt','imgDepth','rgbTest');
load('featurecenter_tsdf.mat');
s =0.1;
[pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,Space,removedpts]=init_gridInd(XYZworldframeTest,s,Qnum);
load('/n/fs/modelnet/rgbd-master/cachedir/release/output/ucm/img_5001.mat', 'ucm2');
sp = bwlabel(ucm2 < 0.35); 
segMaskorg =sp([2,2:2:end,end], [2:2:end,end]); 
% vis_point_cloud(XYZworldframeTest,double([segMaskorg(:),segMaskorg(:),segMaskorg(:)])/max(segMaskorg(:)),30,5000);
% feature = mFeature_segContext_raw(xyzShifted,removedpts,segMaskorg,cellIdx,bincellDim);

%% when pause run following code to visualize 
load(['/n/fs/modelnet/NYUdataSet/NYUdatafeatureNew/' 'feature_' num2str(1) '.mat'],'rgbTest');
rgbTest = rgbTest(~removedpts,:);

ptsofThisSeg =  xyzShifted(segMask ==uniquesegIncell(sIdx),:);
rgbofThisSeg = rgbTest(segMask ==uniquesegIncell(sIdx),:);
ptsIncell = xyzShifted(ptsIncellidx,:);

figure,
vis_point_cloud(ptsofThisSeg,rgbofThisSeg,30,7000);
hold on;
vis_point_cloud(ptsIncell,[],30,5000);
plotcube([0.1,0.1,0.1],cell_center([1,2,3])-0.05,0.1,[1 1 0]);
axis equal;


shapeContext = shapeContext + weight*getShapeContext(cell_center,ptsofThisSeg,distance_threshold,disBinNum,directionsBin);
figure,imagesc(shapeContext)
% for each bin draw sphere
[angleIdx,distanceIdx]=ndgrid(1:size(shapeContext,1),1:size(shapeContext,2));
angleIdx=angleIdx(:);
distanceIdx =distanceIdx(:);
linearSep = (distance_threshold(2)-distance_threshold(1))/disBinNum;    
distance = (distanceIdx -1)*linearSep + distance_threshold(1);
vector = bsxfun(@times,directionsBin(angleIdx,:), distance);
ceters = bsxfun(@plus,vector,cell_center);

scatter3(ceters(:,1),ceters(:,2),ceters(:,3),(shapeContext(:)+mean(shapeContext(:)))/mean(shapeContext(:))*20,'fill');
axis tight;