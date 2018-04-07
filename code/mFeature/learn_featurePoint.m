clear all;
GraphicFea=[];
GraphicSource =[];
KinectFea =[];
KinectFeaSource =[];
KNN =50;
raw =1;
%%
Qnum =6;
feamun =1000;
[startLoc,endLoc]=get_start_end_ind(Qnum,feamun);
[svmpara,featurepara]=setParameters();
s = featurepara.s;
cnt = 1;
%{
filesAll ={'/n/fs/modelnet/warehouseQuery/Data/chair/chair_000000232/chair_000000232.off', ...
           '/n/fs/modelnet/warehouseQuery/Data/chair/chair_000001184/chair_000001184.off',...
           '/n/fs/modelnet/warehouseQuery/Data/chair/chair_000001688/chair_000001688.off',...
           '/n/fs/modelnet/warehouseQuery/Data/chair/chair_000002237/chair_000002237.off',...
           '/n/fs/modelnet/warehouseQuery/Data/night_stand//night_stand_000000207/night_stand_000000207.off',...
           '/n/fs/modelnet/warehouseQuery/Data/night_stand//night_stand_000000028/night_stand_000000028.off',...
           '/n/fs/modelnet/warehouseQuery/Data/night_stand//night_stand_000000639/night_stand_000000639.off',...
           '/n/fs/modelnet/warehouseQuery/Data/night_stand//night_stand_000000164/night_stand_000000164.off',...
           '/n/fs/modelnet/warehouseQuery/Data/toilet/toilet_000000062/toilet_000000062.off',...
           '/n/fs/modelnet/warehouseQuery/Data/toilet//toilet_000000059/toilet_000000059.off',...
           };
Scale = [0.85,0.85,0.8,0.9,0.6,0.6,0.6,0.65,0.8,0.8];
Cdepth =[1.6,2,  2,2.8,2.8,2.8,  4,4,4,5,6];
Cx =    [0,-0.3,0.3,-0.6,0,0.6,-1,0,1,0,0];
Tilt =-0.25;
Rotation =[0:0.25:2];
[fileid,d,r]=ndgrid([1:length(filesAll)],[1:length(Cdepth)],[1:length(Rotation)]);
%%
for k=1:length(fileid(:))
    [depth,K,crop] = off2im(filesAll{fileid(k)},1,Rotation(r(k))*pi,Tilt,Cdepth(d(k)),Cx(d(k)),-0.8,[1;1;1]*Scale(fileid(k)),mod(k,5)==1,0.1,0);
    %off2im(offfile,       ratio,xzRot ,tilt , objz,objx,objy,modelsize,addfloor,enlargefloor)
    [rgb,XYZworldframe]=read_3d_pts_general(depth,K,size(depth),[],crop);
    XYZworldframe(sum(isnan(XYZworldframe),2)>0,:)=[];
    Rtilt =  getRotationMatrix('tilt',Tilt);
    XYZworldframe = (Rtilt*XYZworldframe')';
    if ~isempty(XYZworldframe)
        [pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,Space]=init_gridInd(XYZworldframe,s,Qnum,[0,0,0]);
        [feature,~]= mFeature_points(gridIdx,remainder,bincellDim,histSize,KNN,KNN,featurepara,raw); 
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

a=[GraphicFea,KinectFea]';

stickfea =a(:,1:size(startLoc,2));
countfea =a(:,size(endLoc,2)+1:end);
removedind = sum(a,2)==0;
stickfea(removedind,:)=[];
countfea(removedind,:)=[];
GraphicSource(removedind,:)=[];
ptsAll(removedind) =[];

[IDX1,stick_center] = kmeans(stickfea,49,'EmptyAction', 'singleton');
[IDX2,count_center] = kmeans(countfea,49,'EmptyAction', 'singleton');
save('./featurecenters/learnfeature_point_2','stickfea','countfea','GraphicSource','stick_center','count_center','IDX1','IDX2','ptsAll','-v7.3')
%}
%%
for k=1:100
    n =k
    tmp = load(['/n/fs/modelnet/NYUdataSet/NYUdatafeatureNew/' 'feature_' num2str(n) '.mat'],'XYZworldframeTest','SpaceTest','Rtilt');     
    [pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,Space]=init_gridInd(tmp.XYZworldframeTest,s,Qnum,[0,0,0]);
    [feature,~]= mFeature_points(gridIdx,remainder,bincellDim,histSize,KNN,KNN,featurepara,raw); 
    featureDim =size(feature,4);
    linearIndPts =sub2ind(bincellDim,cellIdx(:,2),cellIdx(:,1),cellIdx(:,3));
    
    [c,~]  = shiftdim(feature,3);
    featureReshaped=reshape(c,featureDim,[]);
    sampleind =randi(size(featureReshaped,2),[1,3000]);    
    sampleind = datasample(unique(linearIndPts),min(300,length(unique(linearIndPts))),'Replace',false);
    
   
    for i =1:length(sampleind),
        if sum(linearIndPts==sampleind(i))>10
         KinectFea =[KinectFea,featureReshaped(:,sampleind(i))];
         KinectFeaSource  = [KinectFeaSource;k,sampleind(i)];
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
a=[GraphicFea,KinectFea]';
save('./featurecenters/learnfeature_point_Knect','a','KinectFeaSource','ptsAll','-v7.3')

stickfea =a(:,1:size(startLoc,2));
countfea =a(:,size(endLoc,2)+1:end);
% removedind = sum(a,2)==0;
% stickfea(removedind,:)=[];
% countfea(removedind,:)=[];
% KinectFeaSource(removedind,:)=[];
% ptsAll(removedind) =[];
[IDX1,stick_center] = kmeans(stickfea,49,'EmptyAction', 'singleton');
[IDX2,count_center] = kmeans(countfea,49,'EmptyAction', 'singleton');
save('./featurecenters/learnfeature_point_Knect','stickfea','countfea','KinectFeaSource','stick_center','count_center','IDX1','IDX2','ptsAll','-v7.3')
%%
%{
% test different size of center 
load('./featurecenters/learnfeature_point_Knect')
numofcenter = 1000;
[IDX1,stick_center] = kmeans(stickfea,numofcenter,'EmptyAction', 'singleton');
[IDX2,count_center] = kmeans(countfea,numofcenter,'EmptyAction', 'singleton');
save(['./featurecenters/learnfeature_point_Knect' num2str(numofcenter)],'stick_center','count_center','IDX1','IDX2','-v7.3')

%}
%%

clear all
load('./featurecenters/learnfeature_point_Knect')
load('./featurecenters/learnfeature_point_Knect1000')
feacenter = count_center;IDX = IDX2;Fea = countfea';
feacenter = stick_center;IDX = IDX1;Fea = stickfea';
clear centerind distanceAll;
for gid = 1:size(feacenter)
    thisgroupIDX = find(IDX == gid);
%   distance = bsxfun(@minus,Fea(:,thisgroupIDX)',feacenter(gid,:));
%   distance = sum(distance.*distance,2);
    distance = pdist2(Fea(:,thisgroupIDX)',feacenter(gid,:));
    [D,I]=sort(distance);
    centerind(gid) = {thisgroupIDX(I)}; 
    distanceAll(gid) = {D};
end


for i = 1:length(centerind)
        pts = ptsAll{centerind{i}(1)};
        %subplot(round(size(feacenter,1)/10),10,i)
        if ~isempty(pts)
            f = figure('visible','off');
            scatter3(pts(1,:),pts(2,:),pts(3,:),'+');
        end
        axis equal;
        axis([0,.1,0,.1,0,.1])
        saveas(f, ['./debug/centers2/' num2str(i),'.jpg'])
        close(f)
end
%%
figure,
cnt =1;
numcol =10;
numrow = 2;
for i =10
    for k = 1:min(numcol*numrow,size(centerind{i},1))
        pts = ptsAll{centerind{i}(k)};
        subplot(2*numrow,numcol,cnt)
        if ~isempty(pts)
            scatter3(pts(1,:),pts(2,:),pts(3,:),'+');
        end
        axis equal;
        axis([0,.1,0,.1,0,.1])
        cnt= cnt+1;
        title(sprintf('%d,%.0f',k,distanceAll{i}(k)))
    end
    cnt =numcol*numrow+1;
    for k = max(1,size(centerind{i},1)-numcol*numrow)+1:size(centerind{i},1)
        pts = ptsAll{centerind{i}(k)};
        subplot(2*numrow,numcol,cnt)
        if ~isempty(pts)
            scatter3(pts(1,:),pts(2,:),pts(3,:),'+');
        end
        axis equal;
        axis([0,.1,0,.1,0,.1])
        cnt= cnt+1;
        title(sprintf('%d,%.0f',k,distanceAll{i}(k)))
    end
end
%%
figure,
plot(Fea(:,[1,265,261]));

%}
%{
count_center = [zeros(1,36);count_center];
stick_center = [zeros(1,feamun);stick_center];
% Save centers
save('featurecenter_point2','stick_center','count_center','Qnum','feamun','IDX1','IDX2');
%}

%}
