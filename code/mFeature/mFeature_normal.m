function [feature,normhistMatrix,normals] = mFeature_normal(pts,cellIdx,bincellDim,KNN,featurepara)

load('./featurecenters/noNomalizedcenter/featurecenter_normal.mat');

% caculte normal
directions = icosahedron2sphere(1); 
directions = directions(directions(:,2) <= 0,:);
normals = points2normals(pts);
[~, idx] = max(abs(directions * normals),[],1);

% build up the histogram
normhistMatrix = zeros([bincellDim,size(directions,1)]);
tmpIdx =[cellIdx(:,2),cellIdx(:,1),cellIdx(:,3),idx'];
weight =ones(size(tmpIdx,1),1);
normhistMatrix = accumarray(tmpIdx,weight,size(normhistMatrix));

% compute feature in cell
feature = zeros([bincellDim,size(normal_center,1)]);
[ind1,ind2,ind3]=ndgrid(1:size(feature,1),1:size(feature,2),1:size(feature,3));
for i=1:length(ind1(:))
     normalfea_raw = normhistMatrix(ind1(i),ind2(i),ind3(i),1:end);
     if featurepara.normalize||featurepara.highdim
         normalfea_norm = (normalfea_raw(:)-norm_mean(:))./norm_scale(:);
     else
         normalfea_norm = normalfea_raw;
         stdf =80;
     end
     normalfea_code = sparse_coding(normalfea_norm(:)',normal_center,stdf,KNN);
     feature(ind1(i),ind2(i),ind3(i),1:end)  = normalfea_code(:);
end

end