clear all;
Qnum =6;
s =0.1;
GraphicFea=[];
%%
Modelpath =  '/n/fs/modelnet/TrainedModels/toilet_Feb13_all/toilet/toilet_000000062/toilet_000000062_r0.3_s0.8_d2.6_t-0.3_x0.5_h-0.7_sft0.0500_asp111.mat';
load(Modelpath);
offfile ='/n/fs/modelnet/warehouseQuery/Data/chair/chair_000000358/chair_000000358.off';
[depth,K,crop,bb3d,vmat,meshObj,BBTight,segMask] = off2im_sys(offfile,1,Model.Rotation*pi,Model.Tilt,Model.Cx,Model.Cdepth+1,-0.95,[1;1;1]*Model.Scale,1,1,-1);
[rgb,XYZworldframe] = read_3d_pts_general(depth,K,size(depth),[],crop);
Rtilt =  getRotationMatrix('tilt',Model.Tilt);
XYZworldframe = (Rtilt*XYZworldframe')';
[pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,SpacePre,removedpts]=init_gridInd(XYZworldframe,s,Qnum,[0,0,0]);
feature = mFeature_curvature_raw(xyzShifted,cellIdx,bincellDim);
fdim = size(feature,4);
enlargeBB =0.1;
[Space,pts,bb3dfea,enlarge]= initModelBB(XYZworldframe,s,enlargeBB,bb3d,SpacePre);
feature  =getTemplate_3d(feature,bb3dfea);
[c,~]  = shiftdim(feature,3);
c=reshape(c,fdim,[]);
c=c(:,randi(size(c,2),[1,300]));
GraphicFea=[GraphicFea,c];

%%
offfile =['/n/fs/modelnet/warehouseQuery/Data/sofa/sofa_000000026/sofa_000000026.off'];
[depth,K,crop,bb3d,vmat,meshObj,BBTight,segMask] = off2im_sys(offfile,1,Model.Rotation*pi,Model.Tilt,Model.Cx,Model.Cdepth+1,-0.95,[1;1;1]*Model.Scale,1,1,-1);
[rgb,XYZworldframe] = read_3d_pts_general(depth,K,size(depth),[],crop);
Rtilt =  getRotationMatrix('tilt',Model.Tilt);
XYZworldframe = (Rtilt*XYZworldframe')';
[pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,SpacePre,removedpts]=init_gridInd(XYZworldframe,s,Qnum,[0,0,0]);
feature = mFeature_curvature_raw(xyzShifted,cellIdx,bincellDim);

enlargeBB =0.1;
[Space,pts,bb3dfea,enlarge]= initModelBB(XYZworldframe,s,enlargeBB,bb3d,SpacePre);
feature  =getTemplate_3d(feature,bb3dfea);

[c,~]  = shiftdim(feature,3);
c=reshape(c,fdim,[]);
c=c(:,randi(size(c,2),[1,300]));
GraphicFea=[GraphicFea,c];

filename = '/n/fs/modelnet/warehouseQuery/Data/toilet//toilet_000000059/toilet_000000059.off';
[depth,K,crop,bb3d,vmat,meshObj,BBTight,segMask] = off2im_sys(offfile,1,Model.Rotation*pi,Model.Tilt,Model.Cx,Model.Cdepth+1,-0.95,[1;1;1]*Model.Scale,1,1,1);
[rgb,XYZworldframe] = read_3d_pts_general(depth,K,size(depth),[],crop);
Rtilt =  getRotationMatrix('tilt',Model.Tilt);
XYZworldframe = (Rtilt*XYZworldframe')';
[pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,SpacePre,removedpts]=init_gridInd(XYZworldframe,s,Qnum,[0,0,0]);
feature = mFeature_curvature_raw(xyzShifted,cellIdx,bincellDim);

enlargeBB =0.1;
[Space,pts,bb3dfea,enlarge]= initModelBB(XYZworldframe,s,enlargeBB,bb3d,SpacePre);
feature  =getTemplate_3d(feature,bb3dfea);

[c,~]  = shiftdim(feature,3);
c=reshape(c,fdim,[]);
c=c(:,randi(size(c,2),[1,300]));
GraphicFea=[GraphicFea,c];

%%

KinectFea =[];
sqc =randi(1074,1,300);
 load('featurecenter_tsdf.mat');
 s =0.1
for i=1:length(sqc)
    i/length(sqc)
    n =sqc(i);
    tmp = load(['/n/fs/modelnet/NYUdataSet/NYUdatafeatureNew/' 'feature_' num2str(n) '.mat'],'XYZworldframeTest','SpaceTest','Rtilt'); 
    Rtilt = tmp.Rtilt;
  
    XYZworldframeTest =tmp.XYZworldframeTest;
    segMask = getSegMaskNYU(n);
    KNN =50;
    [pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,Space,removedpts]=init_gridInd(XYZworldframeTest,s,Qnum,[0,0,0]);
    featureNeg = mFeature_curvature_raw(xyzShifted,cellIdx,bincellDim);
    [b,~]  = shiftdim(featureNeg,3);
    b=reshape(b,fdim,[]);
    b=b(:,randi(size(b,2),[1,1000]));
    KinectFea =[KinectFea,b];

end
a=[GraphicFea,KinectFea]';
save('learnfeature_curv','a','-v7.3')
a(sum(a,2)==0,:)=[];
normscale = max(a,[],2)+1;
a_norm = bsxfun(@rdivide, a, normscale);
[IDX,curv_center] = kmeans(a,49,'EmptyAction', 'singleton','MaxIter',500);% N-by-P data matrix X
curv_center =[zeros(1,size(curv_center,2));curv_center];
save('featurecenter_curv','curv_center','IDX');

for i =1:10:1000,
    fearaw = a(randi(size(a,1),1),:);
    fea1 =pdist2(fearaw(1:1000),stick_center); 
    stdf(i) = std(fea1);
end
save('featurecenter_curv','curv_center','stdf','IDX1','normscale');
