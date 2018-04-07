% visual_planeFea
clear all;
load(['/n/fs/modelnet/NYUdataSet/NYUdatafeatureNew/' 'feature_' num2str(1) '.mat'],'XYZworldframeTest','SpaceTest','featureTest','Rorg','Rtilt','imgDepth','rgbTest');
load('featurecenter_tsdf.mat');
s =0.1;
KNN =50;
[pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,Space,removedpts]=init_gridInd(XYZworldframeTest,s,Qnum);
[~,normhistMatrix,normals] = mFeature_normal(pts,cellIdx,bincellDim,KNN);

% feature = mFeature_plane_raw(normhistMatrix,xyzShifted,removedpts,cellIdx,bincellDim,size(imgDepth),normals,KNN);
close all;
load(['/n/fs/modelnet/NYUdataSet/NYUdatafeatureNew/' 'feature_' num2str(1) '.mat'],'rgbTest');
rgbTest = rgbTest(~removedpts,:);

figure,
vis_point_cloud(xyzShifted,rgbTest,30,7000);
%vis_point_cloud(xyzShifted,[],30,7000);
hold on;
vis_point_cloud(xyzShifted(onPlaneIdx,:),repmat([0,1,0],[sum(onPlaneIdx),1]),30,5000);
vis_point_cloud(xyzShifted(onPlaneIdx&conectedIdx,:),repmat([1,0,0],[sum(onPlaneIdx&conectedIdx),1]),30,5000);
vis_point_cloud(ptsIncell,[],30,5000);
plotcube([0.1,0.1,0.1],cell_center([1,2,3])-0.05,0.1,[1 1 0]);


quiver3(cell_center(1),cell_center(2),cell_center(3),startLineA(1),startLineA(2),startLineA(3),'-r','Linewidth',1);
quiver3(cell_center(1),cell_center(2),cell_center(3),startLineB(1),startLineB(2),startLineB(3),'-g','Linewidth',1);
quiver3(cell_center(1),cell_center(2),cell_center(3),planeNormal(1),planeNormal(2),planeNormal(3),'-b','Linewidth',1);
axis equal;

figure,
subplot(1,3,1)
imagesc(Mask)
axis image
subplot(1,3,2)
imagesc(labelMask)
axis image
subplot(1,3,3)
imagesc(labelMask == labelMask(IJcenter(1),IJcenter(2)));
axis image
%%
projectionOnPlane = bsxfun(@minus,ptsOnConectedplane,cell_center) - disToPlane(onPlaneIdx&conectedIdx,:)*planeNormal;
vector = bsxfun(@plus,cell_center',lineDir.*repmat(planeShape(:)',[3,1]))';
vectorOnplane = bsxfun(@minus,vector,cell_center) - bsxfun(@minus,vector,cell_center)*planeNormal'*planeNormal;
figure,
plot(projectionOnPlane(:,1),projectionOnPlane(:,2),'.','Linewidth',10);
hold on;
plot(0,0,'.r','Linewidth',10)
plot(vectorOnplane(:,1),vectorOnplane(:,2),'rx','Linewidth',10);

figure,
vis_point_cloud(xyzShifted(onPlaneIdx&conectedIdx,:),repmat([0,0,1],[sum(onPlaneIdx&conectedIdx),1]),30,5000);
hold on;
quiver3(cell_center(1)*ones(size(lineDir(1,:))),cell_center(2)*ones(size(lineDir(1,:))),cell_center(3)*ones(size(lineDir(1,:))),...
        lineDir(1,:),lineDir(2,:),lineDir(3,:))
plot3(vector(:,1),vector(:,2),vector(:,3),'xg','Linewidth',10)
plot3(cell_center(:,1),cell_center(:,2),cell_center(:,3),'xr','Linewidth',10)
axis equal;
axis tight;