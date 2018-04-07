clear all
[svmpara,featurepara]=setParameters();
model =0;
enlargewall =0;enlargefloor =0.6;
featurepara.maxdistance =2;
modelname = 'night_stand_000000028_r0.2_s0.62_d2.8_t-0.2_x0.6_h-0.7_sft000_asp111';
Qnum =10; s=0.1;KNN =50;
if model
    imageNum =0;
    % for each cell
    load(['/n/fs/modelnet/TrainedModels/tree0kernel0selfoccweigth1occ_0occNonlinear_0care_1en0.1fea_PlaneCode/night_stand/night_stand_000000028/' modelname '.mat'])
    offfile ='/n/fs/modelnet/warehouseQuery/Data_done/night_stand/night_stand_000000028/night_stand_000000028.off';

    [depth,K,crop,bb3d,vmat,meshObj,BBTight,segMask] = off2im_sys(offfile,1,Model.Rotation*pi,Model.Tilt,Model.Cx,Model.Cdepth,Model.Cheight,[1;1;1]*Model.Scale,1,enlargefloor,enlargewall);
    [rgb,XYZworldframe] = read_3d_pts_general(depth,K,size(depth),[],crop);
    Rtilt =  getRotationMatrix('tilt',Model.Tilt);
    XYZworldframe = (Rtilt*XYZworldframe')';
    
    
    [pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,Space,removedpts]=init_gridInd(XYZworldframe,s,Qnum,[0,0,0]);
    [planeSeg]= getPlanSeg_complete(XYZworldframe,Space,size(depth),imageNum);
    [~,normhistMatrix,normals] = mFeature_normal(pts,cellIdx,bincellDim,KNN,featurepara);
    feature = mFeature_planeSeg_raw(xyzShifted,cellIdx,bincellDim,normals,planeSeg,featurepara);
    feature = mFeature_plane_raw(normhistMatrix,xyzShifted,removedpts,cellIdx,bincellDim,size(depth),normals);
    rgbTest =rgb;
else
    imageNum = 770;
    load(['/n/fs/modelnet/NYUdataSet/NYUdatafeatureNew/' 'feature_' num2str(imageNum) '.mat'],'XYZworldframeTest','SpaceTest','Rorg','Rtilt','imgDepth','rgbTest');
    %load('/n/fs/modelnet/Output/MatResult/chair_000000358pointtsdfnormalshapePlaneRaw_out_pft0_addout_1_missing_0_addmissing_1_occ_1enlarge0/chair/chair_000000358/chair_000000358_r0.7_s1_d2_t-0.3201_x-0.3_h-0.95/Test00384-chair_000000358_r0.7_s1_d2_t-0.3201_x-0.3_h-0.95.mat')
    resultFolder = '/n/fs/modelnet/Output/MatResult/tree0kernel0selfoccweigth1occ_0occNonlinear_0care_1en0.1fea_PlaneCode0/';
    load([resultFolder 'night_stand/night_stand_000000028/' modelname sprintf('/Test%05d-',imageNum) modelname '.mat'])
    
    [pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,Space,removedpts]=init_gridInd(XYZworldframeTest,s,Qnum,[0,0,0]);
    isInBB = ptsInBB(XYZworldframeTest,result.bb3d(1,:));
    XYZworldframe = XYZworldframeTest(isInBB&~removedpts,:);
    rgbTest = rgbTest(isInBB&~removedpts,:);
    [~,normhistMatrix,normals] = mFeature_normal(pts,cellIdx,bincellDim,KNN,featurepara); 
    imageSize = size(imgDepth);
    featureTest = mFeature_plane_raw(normhistMatrix,xyzShifted,removedpts,cellIdx,bincellDim,size(imgDepth),normals,featurepara);
    %load(['/n/fs/modelnet/NYUdataSet/NYUdatafeature_s0.1_1/' 'feature_' num2str(imageNum) '.mat'],'featureTest');
    feature  = getTemplate_3d( featureTest,result.bbw_3d_f(1,:));
    normhistMatrix= getTemplate_3d( normhistMatrix,result.bbw_3d_f(1,:));
end
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
for i =5:size(ind1(:))
    planeShape = feature(ind1(i),ind2(i),ind3(i),end-35:end);
    if sum(planeShape(:)) >0.36,
       [~,normalIdx] = max(normhistMatrix(ind1(i),ind2(i),ind3(i),:));
         planeNormal = directions(normalIdx,:);
         spaceface =[0,0,1;0 1 0];
         [~,facechoose] = min(abs(spaceface*planeNormal'));
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
    [~,normalIdx] = max(normhistMatrix(ind1(i),ind2(i),ind3(i),1:end));
             planeNormal = directions(normalIdx,:);

   
    cell_center =s*([ind2(i),ind1(i),ind3(i)]-0.05);
    planeShape(planeShape(:)>featurepara.maxdistance) = featurepara.maxdistance;
    vector = bsxfun(@plus,cell_center',lineDir.*repmat(planeShape(:)'-0.01,[3,1])*s*0.5)';   
    vectorall =[vectorall;vector;cell_center];
    rgball = [rgball;repmat(rand(1,3),[36,1]);zeros(1,3)];
    
    %% visulize cell and high light 
    figure,
    onPlaneThreshold =0.05;
    disThreshold = Inf;
    vis_point_cloud(XYZworldframe,rgbTest,30,7000);
    hold on;
    ptsIncellidx = find(cellIdx(:,2)==ind1(i)+result.bbw_3d_f(1,2)-1&cellIdx(:,1)== ind2(i)+result.bbw_3d_f(1,1)-1&cellIdx(:,3)==ind3(i)+result.bbw_3d_f(1,3)-1);
    ptsIncell = xyzShifted(ptsIncellidx,:);
    cell_center = nanmean(ptsIncell,1); 
       
    [~,pointcloseTocenterIdx] = pdist2(ptsIncell,cell_center,'euclidean','Smallest',1);
    planCenter = ptsIncell(pointcloseTocenterIdx,:);
    [I,J] = ndgrid(1:imageSize(1),1:imageSize(2));
    I =I(~removedpts);
    J =J(~removedpts);
    IJcenter = [I(ptsIncellidx(pointcloseTocenterIdx)),J(ptsIncellidx(pointcloseTocenterIdx))];
    disToPlane = abs(bsxfun(@minus,xyzShifted,planCenter)*planeNormal');
    normalAgree = abs(planeNormal*normals);
    onPlaneIdx = (disToPlane< onPlaneThreshold)&(normalAgree'>0.8);
    % conected component 
    Mask = zeros(imageSize);
    Mask(IJcenter(1),IJcenter(2)) =1;
    Mask(sub2ind(size(Mask),I(onPlaneIdx),J(onPlaneIdx))) =1; 

    labelMask = bwlabel(Mask,8);
    conectedIdx = labelMask == labelMask(IJcenter(1),IJcenter(2));
    conectedIdx = conectedIdx(~removedpts);

    
    vis_point_cloud(xyzShifted(disToPlane< onPlaneThreshold,:),repmat([0,0,1],[sum(disToPlane< onPlaneThreshold),1]),30,5000);
     vis_point_cloud(xyzShifted(normalAgree'>0.8,:),repmat([0,1,1],[sum(normalAgree'>0.8),1]),30,5000);
   
    vis_point_cloud(xyzShifted(onPlaneIdx,:),repmat([0,1,0],[sum(onPlaneIdx),1]),30,5000);
    vis_point_cloud(xyzShifted(onPlaneIdx&conectedIdx,:),repmat([1,0,0],[sum(onPlaneIdx&conectedIdx),1]),30,5000);
    vis_point_cloud(ptsIncell,[],30,5000);
    plotcube([0.1,0.1,0.1],cell_center([1,2,3])-0.05,0.1,[1 1 0]);
        quiver3(cell_center(1),cell_center(2),cell_center(3),startLineA(1),startLineA(2),startLineA(3),'-r','Linewidth',1);
    quiver3(cell_center(1),cell_center(2),cell_center(3),startLineB(1),startLineB(2),startLineB(3),'-g','Linewidth',1);
    quiver3(cell_center(1),cell_center(2),cell_center(3),planeNormal(1),planeNormal(2),planeNormal(3),'-b','Linewidth',1);
    axis equal;
    pause;
    
    close all;
    end
end
points2ply( ['iFeature/PlaneN' num2str(imageNum) '.ply'],  vectorall', uint8(255*rgball));
points2ply( ['iFeature/orgwholeN' num2str(imageNum) '.ply'],  bsxfun(@plus,XYZworldframe',-1*[Space.Rx(1);Space.Ry(1);Space.Rz(1)]), uint8(255*rgbTest));
