clear all;
Qnum =6;
feamun =1000;
GraphicFea=[];
for m=1
    m
    [depth,K,crop] = off2im('chair_000009.off',1,0.1*pi,-0.25,2.5,-0.9,-0.9,[1;1;1],1,0.1);
    %off2im(offfile,       ratio,xzRot ,tilt , objz,objx,objy,modelsize,addfloor,enlargefloor)
    [rgb,XYZworldframe]=read_3d_pts_general(depth,K,size(depth),[],crop);
    ptsInBB =XYZworldframe';
    s =0.1;
    Space= initSpace(XYZworldframe,s);

    [feature,bincellDim]= compute_normal_fea_raw(XYZworldframe,Space,Qnum,[]);
    [c,~]  = shiftdim(feature,3);
    c=reshape(c,25,[]);
    c=c(:,randi(size(c,2),[1,9000]));
    GraphicFea=[GraphicFea,c];
end

KinectFea =[];
sqc =[2,5,7,20,67,300,50,37,1003,898,455,1001,560,454,333,222,111,367,456,234,877,577,677,344,233,643];


parfor i=1:25
    n =sqc(i)
    tmp = load(['/n/fs/modelnet/NYUdataSet/NYUdatafeatureNew/' 'feature_' num2str(n) '.mat'],'XYZworldframeTest','SpaceTest','Rtilt');   
    %[featureNeg,~]=compute_3d_fea_tsdf_raw(tmp.imgDepth,tmp.Rtilt,tmp.SpaceTest,startLoc,endLoc,Qnum,[]);
    [featureNeg,~]=compute_normal_fea_raw(XYZworldframe,Space,Qnum,[]);
    [b,~]  = shiftdim(featureNeg,3);
    b=reshape(b,25,[]);
    b=b(:,randi(size(b,2),[1,8000]));
    KinectFea =[KinectFea,b];
end
a=[GraphicFea,KinectFea]';
save('learnfeature_normal','a','-v7.3')
a(sum(a,2)==0,:)=[];
[IDX,normal_center] = kmeans(a,49,'EmptyAction', 'singleton');% N-by-P data matrix X
normal_center = [zeros(1,25);normal_center];
save('featurecenter_normal','normal_center','IDX');

