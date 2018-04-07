function [feature,occcell,missDepthcell,occsource] = mFeature_tsdf(depthMap,histSize,bincellDim,Rtilt,crop,Space,KNN_stick,KNN_count,featurepara)
load('./featurecenters/noNomalizedcenter/featurecenter_tsdf.mat'); 

stickfeamun =size(stick_center,2);
[startInd,endInd]=get_start_end_ind(Qnum,stickfeamun);

stepSize =Space.s/Qnum;
mu =0.04;

[X,Y,Z] = ndgrid(1:histSize(1),1:histSize(2),1:histSize(3));
gridIdx = sub2ind(size(X),X,Y,Z); 
gridIdx = gridIdx(:);
gridCoord = bsxfun(@times,[X(:) Y(:) Z(:)]-1,[stepSize stepSize stepSize]);
gridCoord = bsxfun(@plus,gridCoord,[Space.Rx(1),Space.Ry(1),Space.Rz(1)]+0.5*[stepSize stepSize stepSize]);
gridDists = mu * ones(size(X));
gridMissDepth = zeros(size(X));
gridoccSource = zeros(size(X));
clear X Y Z;

% projection
[gridProj,gridProjDepth] = project3dPtsTo2d(gridCoord,Rtilt,crop);   
gridProj = round(gridProj);

% remove invalid projections
isValid = gridProj(:,1) >= 1 & gridProj(:,1) <= size(depthMap,2) & ...
          gridProj(:,2) >= 1 & gridProj(:,2) <= size(depthMap,1) & ...
          gridProjDepth > 0;
gridIdx = gridIdx(isValid);
gridCoord = gridCoord(isValid,:);
gridProj = gridProj(isValid,:);
gridProjDepth = gridProjDepth(isValid);

 % remove point that has no depth
gridProjIdx = sub2ind(size(depthMap),gridProj(:,2),gridProj(:,1));
isValid = ~isnan(depthMap(gridProjIdx)) & depthMap(gridProjIdx)~=0;
gridMissDepth(gridIdx(~isValid)) =1;
gridIdx = gridIdx(isValid);
gridCoord = gridCoord(isValid,:);
gridProjIdx = gridProjIdx(isValid,:);
gridProjDepth = gridProjDepth(isValid);

% project depth map to 3d point cloud
try 
    points3d = rgb_plane2rgb_world(depthMap);
catch
    [~,points3d] = read_3d_pts_general(depthMap,[],size(depthMap),[],crop);
end
points3d = (Rtilt * points3d')';

% compute distance
dists = sqrt(sum((gridCoord-points3d(gridProjIdx,:)).^2,2));
dists(gridProjDepth > depthMap(gridProjIdx)) = -dists(gridProjDepth > depthMap(gridProjIdx)); % wheter is grid is before or after the points
gridDists(gridIdx) = dists;
gridDists = gridDists / mu;
gridDists = max(-1,gridDists);
gridDists = min(1,gridDists);

% record the world cordinate depth value 
gridoccSource(gridIdx) = points3d(gridProjIdx,2);

% to maske it consistent with other part of code.
gridDists = permute(gridDists,[2 1 3]);
gridMissDepth =permute(gridMissDepth,[2,1,3]);
gridoccSource =permute(gridoccSource,[2,1,3]);

ptxMatrix =gridDists;
feature = zeros([bincellDim,size(stick_center,1)+size(count_center,1)]);
occcell  =zeros(bincellDim);
occsource = zeros(bincellDim);
missDepthcell = zeros(bincellDim);
[ind1,ind2,ind3]=ndgrid(1:size(feature,1),1:size(feature,2),1:size(feature,3));
%cell = shrinkHist(ptxMatrix,6);
for i=1:length(ind1(:))      
    cellmatrix = ptxMatrix((ind1(i)-1)*Qnum+1:min(ind1(i)*Qnum,size(ptxMatrix,1)),...
                              (ind2(i)-1)*Qnum+1:min(ind2(i)*Qnum,size(ptxMatrix,2)),...
                              (ind3(i)-1)*Qnum+1:min(ind3(i)*Qnum,size(ptxMatrix,3)));
                        
    misscellmatrix = gridMissDepth((ind1(i)-1)*Qnum+1:min(ind1(i)*Qnum,size(ptxMatrix,1)),...
                              (ind2(i)-1)*Qnum+1:min(ind2(i)*Qnum,size(ptxMatrix,2)),...
                              (ind3(i)-1)*Qnum+1:min(ind3(i)*Qnum,size(ptxMatrix,3)));
    occSourcecellmatrix = gridoccSource((ind1(i)-1)*Qnum+1:min(ind1(i)*Qnum,size(ptxMatrix,1)),...
                              (ind2(i)-1)*Qnum+1:min(ind2(i)*Qnum,size(ptxMatrix,2)),...
                              (ind3(i)-1)*Qnum+1:min(ind3(i)*Qnum,size(ptxMatrix,3)));                  
    %cellmatrix(1:size(cellmatrix_tmp,1),1:size(cellmatrix_tmp,2),1:size(cellmatrix_tmp,3))= cellmatrix_tmp;
    
    % second order feature
    stickfea_raw=cellmatrix(startInd)-cellmatrix(endInd);
    if featurepara.normalize
        stickfea_norm =(stickfea_raw-stick_norm_mean)./stick_norm_scale;
    else
        stickfea_norm = stickfea_raw;
        stdf_stick =10;
    end
    stickfea_code = sparse_coding(stickfea_norm,stick_center,stdf_stick,KNN_stick);
    % first order feature
    countfea_raw  = cellmatrix(:)';
    if featurepara.normalize||featurepara.highdim
        count_norm = (countfea_raw-count_norm_mean)./count_norm_scale;
    else
        count_norm = countfea_raw;
        stdf_count = 7;
    end
    countfea_code = sparse_coding(count_norm,count_center,stdf_count,KNN_count);
    feature(ind1(i),ind2(i),ind3(i),1:end)=[stickfea_code(:);countfea_code(:)]; 
    
    % whether this cell is occluded (>.9 raw tsdf value is -1 )
    occcell(ind1(i),ind2(i),ind3(i)) = sum(countfea_raw(:)+1 <1.0e-15)/size(countfea_raw(:),1)>0.9;
    % where is the source of occlusion ?
    occsource(ind1(i),ind2(i),ind3(i)) = max(occSourcecellmatrix(:));
    % cell is missing when >.9  grids are projectectd to a missing depth
    missDepthcell(ind1(i),ind2(i),ind3(i)) = sum(misscellmatrix(:))/size(misscellmatrix(:),1) >0.9;
end

end