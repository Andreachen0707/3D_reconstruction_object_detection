clear all
[svmpara,featurepara]=setParameters();
model =0;
enlargewall =0;enlargefloor =0.6;
featurepara.maxdistance =2;
modelname = 'night_stand_000000028_r0_s0.62_d2.8_t-0.2_x-0.6_h-0.7_sft000_asp111';
modelname ='night_stand_000000028_r-0.2_s0.62_d2.8_t-0.2_x-0.6_h-0.7_sft000_asp111';
Qnum =10; s=0.1;KNN =50;
enlargeBB =0.1;
experimetnName = 'newfea2tril-1tree0kernel0selfoccweigth1occ_0occNonlinear_0care_1en0.1fea_PlaneCodeSegmaxplane_0.5minplane_0floatCode1';
if model
    imageNum =0;
    load(['/n/fs/modelnet/TrainedModels/' experimetnName '/night_stand/night_stand_000000028/' modelname '.mat'])
    offfile ='/n/fs/modelnet/warehouseQuery/Data_done/night_stand/night_stand_000000028/night_stand_000000028.off';

    [depth,K,crop,bb3d,vmat,meshObj,BBTight,segMask] = off2im_sys(offfile,1,Model.Rotation*pi,Model.Tilt,Model.Cx,Model.Cdepth,Model.Cheight,[1;1;1]*Model.Scale,1,enlargefloor,enlargewall);
    [rgb,XYZworldframe] = read_3d_pts_general(depth,K,size(depth),[],crop);
    Rtilt =  getRotationMatrix('tilt',Model.Tilt);
    XYZworldframe = (Rtilt*XYZworldframe')';
    
    
    [pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,SpacePre,removedpts]=init_gridInd(XYZworldframe,s,Qnum,[0,0,0]);
    [ImgplaneSeg]= getPlanSeg_complete(XYZworldframe,SpacePre,size(depth),imageNum);
    feature = mFeature_planeSeg_raw(xyzShifted,removedpts,cellIdx,bincellDim,size(depth),ImgplaneSeg,featurepara);
    
    result.bbw_3d_f =[1,1,1,1,1,1];
    rgbTest =rgb;
    imageSize = size(depth);
    
    %{
    [Space,XYZworldframe, bb3dfea,Modelenlarge]= initModelBB(XYZworldframe,s,enlargeBB,bb3d,SpacePre);
    feature = getTemplate_3d(feature,bb3dfea);
    feature_coded = linear_coding(feature,featurepara);
    result.bbw_3d_f =bb3dfea;
    
        weight = shiftdim(Model.svm.w,3);
        weight = reshape(weight,size(Model.svm.w,4),[]);
        weight = weight(1:396,:);
    %}
else
    imageNum = 706;
    load(['/n/fs/modelnet/NYUdataSet/NYUdatafeatureNew/' 'feature_' num2str(imageNum) '.mat'],'XYZworldframeTest','SpaceTest','Rorg','Rtilt','imgDepth','rgbTest');
    [pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,Space,removedpts]=init_gridInd(XYZworldframeTest,s,Qnum,[0,0,0]);
    imageSize = size(imgDepth);
    [ImgplaneSeg]= getPlanSeg_complete(XYZworldframeTest,Space,imageSize,imageNum);
    
    featureTest = mFeature_planeSeg_raw(xyzShifted,removedpts,cellIdx,bincellDim,imageSize,ImgplaneSeg,featurepara);
    %load(['/n/fs/modelnet/NYUdataSet/NYUdatafeature_s0.1Seg2/' 'feature_' num2str(imageNum) '.mat'],'featureTest');
    
    resultFolder = ['/n/fs/modelnet/Output/MatResult/' experimetnName '0/'];
    load([resultFolder 'night_stand/night_stand_000000028/' modelname sprintf('/Test%05d-',imageNum) modelname '.mat'])
    feature  = getTemplate_3d( featureTest,result.bbw_3d_f(1,:));
    isInBB = ptsInBB(XYZworldframeTest,result.bb3d(1,:));
    XYZworldframe = XYZworldframeTest(isInBB&~removedpts,:);
    rgbTest = rgbTest(isInBB&~removedpts,:);
   
end
%}
%%
directions = icosahedron2sphere(1); 
directions = directions(directions(:,2) <= 0,:);
numbin = 36;
angle = 0:2*pi/numbin:2*pi;
angle =angle(1:numbin);
[ind1,ind2,ind3]=ndgrid(1:size(feature,1),1:size(feature,2),1:size(feature,3));
ind1 =ind1(:);ind2 =ind2(:);ind3=ind3(:);

vectorall =[];
rgball=[];
for i =1:size(ind1(:),1)
    planeSeg = ImgplaneSeg(~removedpts);
    planeShape = feature(ind1(i),ind2(i),ind3(i),end-35:end);
    ptsIncellidx = find(cellIdx(:,2)==ind1(i)+result.bbw_3d_f(1,2)-1&cellIdx(:,1)== ind2(i)+result.bbw_3d_f(1,1)-1&cellIdx(:,3)==ind3(i)+result.bbw_3d_f(1,3)-1);
    ptsIncell = xyzShifted(ptsIncellidx,:);
    planIdIncell = planeSeg(ptsIncellidx);
    planIdIncell = planIdIncell(planIdIncell>0);
    numOfsegIncell = length(unique(planIdIncell));
    planIdIncell_pick = mode(planIdIncell(planIdIncell>0));
	
    if sum(planeShape(:)) >0.360001,
       onPlaneIdx = planeSeg==planIdIncell_pick;
       ptsOnConectedplane = xyzShifted(onPlaneIdx,:);
       [C,V]= PCA(ptsOnConectedplane',3);
       planeNormal = C(:,3);
       if [0,1,0]*planeNormal>0,
           planeNormal =-1*planeNormal;
       end
       planeNormal =planeNormal/norm(planeNormal);
        spaceface =[0,0,1;0 1 0];
        [~,facechoose] = min(abs(spaceface*planeNormal));
        startLineA = cross(planeNormal,spaceface(facechoose,:));
        startLineA = startLineA/norm(startLineA);
        startLineB = cross(planeNormal,startLineA);
         
        if facechoose ==1,
            % z up
            startLineA = startLineB(3)/abs(startLineB(3))*startLineA;
            startLineB = startLineB(3)/abs(startLineB(3))*startLineB;
        else
            % x right   
            startLineA = startLineA(1)/abs(startLineA(1))*startLineA;
            startLineB = startLineA(1)/abs(startLineA(1))*startLineB;
        end
             
    lineDir = [startLineA(:),startLineB(:),planeNormal(:)]*[cos(angle);sin(angle);zeros(size(angle))];
    cell_center = nanmean(ptsIncell(planIdIncell == planIdIncell_pick,:),1); 
    %cell_center = nanmean(ptsIncell(planIdIncell == planIdIncell_pick,:),1); 
    planeShape(planeShape(:)>featurepara.maxdistance) = featurepara.maxdistance;
    vector = bsxfun(@plus,cell_center',lineDir.*repmat(planeShape(:)'-0.01,[3,1])*s*0.5)';   
    vectorall =[vectorall;vector;cell_center];
    rgball = [rgball;repmat(rand(1,3),[36,1]);zeros(1,3)];
    %{
    %% visulize cell and high light 
    figure,
    vis_point_cloud(xyzShifted,rgbTest,30,5000);
    hold on;
    vis_point_cloud(xyzShifted(onPlaneIdx,:),repmat([1,0,0],[sum(onPlaneIdx),1]),30,5000);
    vis_point_cloud(ptsIncell,[],30,3000);
    plotcube([0.1,0.1,0.1],cell_center([1,2,3])-0.05,0.1,[1 1 0]);
    quiver3(cell_center(1),cell_center(2),cell_center(3),startLineA(1),startLineA(2),startLineA(3),'-r','Linewidth',1);
    quiver3(cell_center(1),cell_center(2),cell_center(3),startLineB(1),startLineB(2),startLineB(3),'-g','Linewidth',1);
    quiver3(cell_center(1),cell_center(2),cell_center(3),planeNormal(1),planeNormal(2),planeNormal(3),'-b','Linewidth',1);
    
    figure,
    scatter3(ptsOnConectedplane(:,1),ptsOnConectedplane(:,2),ptsOnConectedplane(:,3)); axis equal;
    hold on;
    quiver3(repmat(cell_center(1),[1,36]),repmat(cell_center(2),[1,36]),repmat(cell_center(3),[1,36]),lineDir(1,:),lineDir(2,:),lineDir(3,:),'-r','Linewidth',1);
    boundary = bsxfun(@plus,lineDir.*repmat(planeShape(:)'-0.01,[3,1]),cell_center');
    scatter3(boundary(1,:),boundary(2,:),boundary(3,:),'fill','r'); axis equal;
    axis equal;
    pause;
    
    close all;
    %}
    end
end
points2ply( ['iFeature/PlaneN' num2str(imageNum) '.ply'],  vectorall', uint8(255*rgball));
points2ply( ['iFeature/orgwholeN' num2str(imageNum) '.ply'],  bsxfun(@plus,XYZworldframe',-1*[Space.Rx(1);Space.Ry(1);Space.Rz(1)]), uint8(255*rgbTest));
