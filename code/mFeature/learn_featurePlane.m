clear all;
Qnum =6;
s =0.1;
GraphicFea=[];
KNN =50;

%%
Modelpath =  '/n/fs/modelnet/TrainedModels/Plan0.15_chair_000000358pointtsdfnormalshape_out_pft0_addout_1_missing_0_addmissing_1_occ_1addempty0/chair/chair_000000358/chair_000000358_r0.4_s1_d2_t-0.3201_x-0.3_h-0.95.mat';
load(Modelpath);

offfile ='/n/fs/modelnet/warehouseQuery/Data/chair/chair_000000358/chair_000000358.off';
[depth,K,crop,bb3d,vmat,meshObj,BBTight,segMask] = off2im_sys(offfile,1,Model.Rotation*pi,Model.Tilt,Model.Cx,Model.Cdepth+1,-0.95,[1;1;1]*Model.Scale,1,1,-1);
[rgb,XYZworldframe] = read_3d_pts_general(depth,K,size(depth),[],crop);
Rtilt =  getRotationMatrix('tilt',Model.Tilt);
XYZworldframe = (Rtilt*XYZworldframe')';

[pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,SpacePre,removedpts]=init_gridInd(XYZworldframe,s,Qnum);
[~,normhistMatrix,normals] = mFeature_normal(pts,cellIdx,bincellDim,KNN);
feature = mFeature_plane_raw(normhistMatrix,xyzShifted,removedpts,cellIdx,bincellDim,size(depth),normals);
feaDim =  size(feature,4)
enlargeBB =0.1;
[Space,pts,bb3dfea,enlarge]= initModelBB(XYZworldframe,s,enlargeBB,bb3d,SpacePre);
feature  =getTemplate_3d(feature,bb3dfea);

[c,~]  = shiftdim(feature,3);
c=reshape(c,feaDim,[]);
c=c(:,randi(size(c,2),[1,300]));
GraphicFea=[GraphicFea,c];

%%
offfile =['/n/fs/modelnet/warehouseQuery/Data/sofa/sofa_000000026/sofa_000000026.off'];
[depth,K,crop,bb3d,vmat,meshObj,BBTight,segMask] = off2im_sys(offfile,1,Model.Rotation*pi,Model.Tilt,Model.Cx,Model.Cdepth+1,-0.95,[1;1;1]*Model.Scale,1,1,-1);
[rgb,XYZworldframe] = read_3d_pts_general(depth,K,size(depth),[],crop);
Rtilt =  getRotationMatrix('tilt',Model.Tilt);
XYZworldframe = (Rtilt*XYZworldframe')';
[pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,SpacePre,removedpts]=init_gridInd(XYZworldframe,s,Qnum);
[~,normhistMatrix,normals] = mFeature_normal(pts,cellIdx,bincellDim,KNN);
feature = mFeature_plane_raw(normhistMatrix,xyzShifted,removedpts,cellIdx,bincellDim,size(depth),normals);

enlargeBB =0.1;
[Space,pts,bb3dfea,enlarge]= initModelBB(XYZworldframe,s,enlargeBB,bb3d,SpacePre);
feature  =getTemplate_3d(feature,bb3dfea);

[c,~]  = shiftdim(feature,3);
c=reshape(c,feaDim,[]);
c=c(:,randi(size(c,2),[1,300]));
GraphicFea=[GraphicFea,c];
%%
filename = '/n/fs/modelnet/warehouseQuery/Data/toilet//toilet_000000059/toilet_000000059.off';
[depth,K,crop,bb3d,vmat,meshObj,BBTight,segMask] = off2im_sys(offfile,1,Model.Rotation*pi,Model.Tilt,Model.Cx,Model.Cdepth+1,-0.95,[1;1;1]*Model.Scale,1,1,1);
[rgb,XYZworldframe] = read_3d_pts_general(depth,K,size(depth),[],crop);
Rtilt =  getRotationMatrix('tilt',Model.Tilt);
XYZworldframe = (Rtilt*XYZworldframe')';

[pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,SpacePre,removedpts]=init_gridInd(XYZworldframe,s,Qnum);
[~,normhistMatrix,normals] = mFeature_normal(pts,cellIdx,bincellDim,KNN);
feature = mFeature_plane_raw(normhistMatrix,xyzShifted,removedpts,cellIdx,bincellDim,size(depth),normals);

enlargeBB =0.1;
[Space,pts,bb3dfea,enlarge]= initModelBB(XYZworldframe,s,enlargeBB,bb3d,SpacePre);
feature  =getTemplate_3d(feature,bb3dfea);

[c,~]  = shiftdim(feature,3);
c=reshape(c,feaDim,[]);
c=c(:,randi(size(c,2),[1,300]));
GraphicFea=[GraphicFea,c];

%%

KinectFea =[];
sqc =randi(1074,1,1074);
 load('featurecenter_tsdf.mat');
 s =0.1
%for i=1:length(sqc)
%     i/length(sqc)
%     n =sqc(i);
%     tmp = load(['/n/fs/modelnet/NYUdataSet/NYUdatafeatureNew/' 'feature_' num2str(n) '.mat'],'XYZworldframeTest','SpaceTest','Rtilt','imgDepth'); 
%     Rtilt = tmp.Rtilt;
%     XYZworldframeTest =tmp.XYZworldframeTest;
%     KNN =50;
%     [pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,Space,removedpts]=init_gridInd(XYZworldframeTest,s,Qnum);
%     [~,normhistMatrix,normals] = mFeature_normal(pts,cellIdx,bincellDim,KNN);
%     featureNeg = mFeature_plane_raw(normhistMatrix,xyzShifted,removedpts,cellIdx,bincellDim,size(tmp.imgDepth),normals,KNN);
cnt =0; 
for i =1:1074
    if cnt>5000, break; end
    try
        load(['/n/fs/modelnet/NYUdataSet/NYUdatafeatureNewCombinedWithPlaneRaw_pointcount/' 'feature_' num2str(i) '.mat'])
        featureNeg = featureTest(:,:,:,301:336);    
        [b,~]  = shiftdim(featureNeg,3);
        b=reshape(b,feaDim,[]);
        b=b(:,randi(size(b,2),[1,500]));
        KinectFea =[KinectFea,b];
        cnt =cnt+1
    catch ee
        display(ee)
    end
    
end
a=[GraphicFea,KinectFea]';
save('learnfeature_plane','a','-v7.3')
a(sum(a,2)==0,:)=[];

[IDX,plane_center] = kmeans(a,49,'EmptyAction', 'singleton','MaxIter',500);% N-by-P data matrix X
plane_center =[zeros(1,size(plane_center,2));plane_center];
save('featurecenter_plane','plane_center','IDX');
 

for i =1:1000,
    fearaw = a(randi(size(a,1),1),:);
    fea1 =pdist2(fearaw,plane_center); 
    stdf(i) = std(fea1);
end
