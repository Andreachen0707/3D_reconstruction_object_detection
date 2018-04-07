clear all;
GraphicFea=[];
GraphicSource =[];
KNN =50;
raw =1;
%%
Qnum =6;
filesAll ={'/n/fs/modelnet/warehouseQuery/Data/chair/chair_000000232/chair_000000232.off', ...
           '/n/fs/modelnet/warehouseQuery/Data/chair/chair_000001184/chair_000001184.off',...
           '/n/fs/modelnet/warehouseQuery/Data/chair/chair_000001688/chair_000001688.off',...
           '/n/fs/modelnet/warehouseQuery/Data/chair/chair_000002237/chair_000002237.off',...
           '/n/fs/modelnet/warehouseQuery/Data/night_stand//night_stand_000000207/night_stand_000000207.off',...
           '/n/fs/modelnet/warehouseQuery/Data/night_stand//night_stand_000000028/night_stand_000000028.off',...
           '/n/fs/modelnet/warehouseQuery/Data/toilet/toilet_000000062/toilet_000000062.off',...
           '/n/fs/modelnet/warehouseQuery/Data/toilet//toilet_000000059/toilet_000000059.off',...
           };
Scale = [0.85,0.85,0.8,0.9,0.6,0.65,0.8,0.8];
Cdepth =[1,2,3,4,5,6];
Rotation =0.1:0.3:2;
[fileid,d,r]=ndgrid([1:length(filesAll)],[1:length(Cdepth)],[1:length(Rotation)]);
[svmpara,featurepara]=setParameters();
s = featurepara.s;
cnt = 1;
%%
for k=1:length(fileid(:))
    [depth,K,crop] = off2im(filesAll{fileid(k)},1,Rotation(r(k))*pi,-0.25,Cdepth(d(k)),-0.9,-0.9,[1;1;1]*Scale(fileid(k)),mod(k,5)==1,0.1,0);
    %                off2im(offfile,        ratio,xzRot ,tilt , objz,objx,objy,modelsize,addfloor,enlargefloor)
    [~,XYZworldframe]=read_3d_pts_general(depth,K,size(depth),[],crop);
    XYZworldframe(sum(isnan(XYZworldframe),2)>0,:)=[];
    if ~isempty(XYZworldframe)
    [pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,Space]=init_gridInd(XYZworldframe,s,Qnum,[0,0,0]);
    feature = mFeature_shape(bincellDim,cellIdx,xyzShifted,Space,Qnum,KNN,featurepara,raw);
    featureDim =size(feature,4);
    linearIndPts =sub2ind(bincellDim,cellIdx(:,2),cellIdx(:,1),cellIdx(:,3));
    
    [c,~]  = shiftdim(feature,3);
    featureReshaped=reshape(c,featureDim,[]);
    sampleind =randi(size(featureReshaped,2),[1,3000]);    
    sampleind = unique(sampleind);
    GraphicFea =[ GraphicFea,featureReshaped(:,sampleind)];
    GraphicSource  = [GraphicSource;k*ones(size(sampleind')),sampleind'];
   
    for i =1:length(sampleind),
         ptsIncell = pts(linearIndPts==sampleind(i),:);
         offset = (cellIdx(linearIndPts==sampleind(i),:) -1)*s;
         if ~isempty(xyzShifted(linearIndPts==sampleind(i),:))
            ptsAll(cnt) ={[xyzShifted(linearIndPts==sampleind(i),:)-offset]'};
         else
            ptsAll(cnt) ={[]};                   
         end
         cnt =cnt+1;
    end
    end
end
a=[GraphicFea]';
removedind = sum(a,2)==0;
a(removedind,:)=[];
GraphicFea(:,removedind)=[];
GraphicSource(removedind,:)=[];
ptsAll(removedind) =[];
[IDX,shape_center] = kmeans(a,49,'EmptyAction', 'singleton');% N-by-P data matrix X
save('./featurecenters/learnfeature_shape_2','GraphicFea','GraphicSource','shape_center','IDX','ptsAll','-v7.3')

%% visualize kmeans
clear all;
load('./featurecenters/learnfeature_shape_2.mat')
feacenter = shape_center;
IDX = IDX;
Fea = GraphicFea;
clear centerind;
for gid = 1:size(feacenter)
    thisgroupIDX = find(IDX == gid);
    distance = bsxfun(@minus,Fea(:,thisgroupIDX)',feacenter(gid,:));
    distance = sum(distance.*distance,2);
    [~,I]=sort(distance);
    centerind(gid,:) = thisgroupIDX(I([1:10,end-5,end])); 
end

figure,
for i = 1:size(centerind,1)
    for k = 1
        pts = ptsAll{centerind(i,k)};
        subplot(5,10,i)
        scatter3(pts(1,:),pts(2,:),pts(3,:),'+');
        axis equal;
        axis([0,.1,0,.1,0,.1])
    end
end

figure,
for i =11
    for k = 1:size(centerind,2)
        pts = ptsAll{centerind(i,k)};
        subplot(1,size(centerind,2),k)
        scatter3(pts(1,:),pts(2,:),pts(3,:),'+');
        axis equal;
        axis([0,.1,0,.1,0,.1])
    end
end
figure,
plot(GraphicFea(:,1:10));

%{
shape_center =[zeros(1,size(shape_center,2));shape_center];
save('featurecenter_shape','shape_center','IDX');
%}
%{
KinectFea =[];
sqc =[2,5,7,20,67,300,50,37,1003,898,455,1001,560,454,333,222,111,367,456,234,877,577,677,344,233,643,700:800];
parfor i=1:120
    n =sqc(i)
    tmp = load(['/n/fs/modelnet/NYUdataSet/NYUdatafeatureNew/' 'feature_' num2str(n) '.mat'],'XYZworldframeTest','SpaceTest','Rtilt');   
    %[featureNeg,~]=compute_3d_fea_tsdf_raw(tmp.imgDepth,tmp.Rtilt,tmp.SpaceTest,startLoc,endLoc,Qnum,[]);
SpaceTest.s = s;
    [featureNeg,bincellDim]=compute_shape_fea_raw(XYZworldframe,Space,Qnum,[]);
    [b,~]  = shiftdim(featureNeg,3);
    b=reshape(b,84,[]);
    b=b(:,randi(size(b,2),[1,4000]));
    KinectFea =[KinectFea,b];
    %{
    load(['/n/fs/modelnet/NYUdataSet/NYUdatafeatureRaw/' 'feature_' num2str(n) '.mat'],'XYZworldframeTest','SpaceTest');   
    featureNeg = compute_3d_fea7_raw(XYZworldframeTest,SpaceTest,Qnum,startLoc,endLoc);
    [b,~]  = shiftdim(featureNeg,3);
    b=reshap(b,feamun+36,[]);
    b=b(:,randi(size(b,2),[1,4000]));
    KinectFea =[KinectFea,b];
    %}
end
a=[GraphicFea,KinectFea]';
save('learnfeature_shape','a','-v7.3')
%}