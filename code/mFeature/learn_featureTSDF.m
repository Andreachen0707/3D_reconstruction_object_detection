clear all;
Qnum =6;
feamun =1000;
GraphicFea=[];
[startLoc,endLoc]=get_start_end_ind(Qnum,feamun);
s =0.1
%{
%%
Modelpath =  '/n/fs/modelnet/TrainedModels/Plan0.15_chair_000000358pointtsdfnormalshape_out_pft0_addout_1_missing_0_addmissing_1_occ_1addempty0/chair/chair_000000358/chair_000000358_r0.4_s1_d2_t-0.3201_x-0.3_h-0.95.mat';
load(Modelpath);
offfile ='/n/fs/modelnet/warehouseQuery/Data/chair/chair_000000358/chair_000000358.off';
[depth,K,crop,bb3d,vmat,meshObj,BBTight,segMask] = off2im_sys(offfile,1,Model.Rotation*pi,Model.Tilt,Model.Cx,Model.Cdepth+1,-0.95,[1;1;1]*Model.Scale,1,1,-1);
[rgb,XYZworldframe] = read_3d_pts_general(depth,K,size(depth),[],crop);
Rtilt =  getRotationMatrix('tilt',Model.Tilt);
XYZworldframe = (Rtilt*XYZworldframe')';
[pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,Space,removedpts]=init_gridInd(XYZworldframe,s,Qnum);
[feature,gridDists] = mFeature_tsdf_raw(depth,histSize,bincellDim,Rtilt,crop,Space);
feaDim =  size(feature,4)
enlargeBB =0.1;
[Space,pts,bb3dfea,enlarge]= initModelBB(XYZworldframe,s,enlargeBB,bb3d,Space);
feature  =getTemplate_3d(feature,bb3dfea);

[c,~]  = shiftdim(feature,3);
c=reshape(c,feaDim,[]);
c=c(:,randi(size(c,2),[1,500]));
GraphicFea=[GraphicFea,c];

%%
offfile =['/n/fs/modelnet/warehouseQuery/Data/sofa/sofa_000000026/sofa_000000026.off'];
[depth,K,crop,bb3d,vmat,meshObj,BBTight,segMask] = off2im_sys(offfile,1,Model.Rotation*pi,Model.Tilt,Model.Cx,Model.Cdepth+1,-0.95,[1;1;1]*Model.Scale,1,1,-1);
[rgb,XYZworldframe] = read_3d_pts_general(depth,K,size(depth),[],crop);
Rtilt =  getRotationMatrix('tilt',Model.Tilt);
XYZworldframe = (Rtilt*XYZworldframe')';
[pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,Space,removedpts]=init_gridInd(XYZworldframe,s,Qnum);
[feature,gridDists] = mFeature_tsdf_raw(depth,histSize,bincellDim,Rtilt,crop,Space);
feaDim =  size(feature,4)
enlargeBB =0.1;
[Space,pts,bb3dfea,enlarge]= initModelBB(XYZworldframe,s,enlargeBB,bb3d,Space);
feature  =getTemplate_3d(feature,bb3dfea);

[c,~]  = shiftdim(feature,3);
c=reshape(c,feaDim,[]);
c=c(:,randi(size(c,2),[1,500]));
GraphicFea=[GraphicFea,c];

%%

KinectFea =[];
sqc =randi(1074,1,200);
parfor i=1:length(sqc)
    n =sqc(i)
    tmp =load(['/n/fs/modelnet/NYUdataSet/NYUdatafeatureNew/' 'feature_' num2str(n) '.mat'],'XYZworldframeTest','SpaceTest','featureTest','Rorg','Rtilt','imgDepth','rgbTest');
    Rtilt = tmp.Rtilt;
    Rorg = tmp.Rorg;
    XYZworldframeTest =tmp.XYZworldframeTest;
    [pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,Space,removedpts]=init_gridInd(XYZworldframeTest,s,Qnum);
    [featureNeg,gridDists] = mFeature_tsdf_raw(tmp.imgDepth,histSize,bincellDim,Rtilt,[],Space);
        
    [b,~]  = shiftdim(featureNeg,3);
    b=reshape(b,feaDim,[]);
    b=b(:,randi(size(b,2),[1,2000]));
    KinectFea =[KinectFea,b];
end
a=[GraphicFea,KinectFea]';
save('learnfeature_tsdf','a','-v7.3')
%}
load('learnfeature_tsdf.mat')
stickfea =a(:,1:size(startLoc,2));
fprintf('number of stick features1 =%d\n',size(stickfea,1));
stickfea(sum(stickfea,2)==0,:)=[];
fprintf('number of stick features2 =%d\n',size(stickfea,1));
[IDX1,stick_center] = kmeans(stickfea,49,'EmptyAction', 'singleton','MaxIter',500);% N-by-P data matrix X
stick_center = [zeros(1,feamun);stick_center];


countfea =a(:,size(endLoc,2)+1:end);
fprintf('number of countfea features1 =%d\n',size(countfea,1));
countfea(sum(countfea-1,2)==0,:)=[];
fprintf('number of countfea features2 =%d\n',size(countfea,1));
countfea(sum(countfea+1,2)==0,:)=[];
fprintf('number of countfea features3 =%d\n',size(countfea,1));
[IDX2,count_center] = kmeans(countfea,48,'EmptyAction', 'singleton','MaxIter',500);% N-by-P data matrix X
count_center = [ones(1,216);-1*ones(1,216);count_center];
%}
% Save centers
for i =1:10:1000,
    fearaw = a(randi(size(a,1),1),:);
    fea1 =pdist2(fearaw(1:1000),stick_center); 
    stdf_stick(i) = std(fea1);
    fea2 =pdist2(fearaw(1001:end),count_center); 
    stdf_count(i) = std(fea2);
end
save('featurecenter_tsdf2','stick_center','count_center','Qnum','feamun','stdf_count','stdf_stick','IDX1','IDX2');
