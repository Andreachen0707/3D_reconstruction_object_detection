function [feature,ptxMatrix]= mFeature_points(gridIdx,remainder,bincellDim,histSize,KNN_stick,KNN_count,featurepara,raw)
load('./featurecenters/noNomalizedcenter/featurecenter_point.mat'); 
stickfeamun =size(stick_center,2);
[startInd,endInd]=get_start_end_ind(Qnum,stickfeamun);
% create gaussian smoothing filter
cutoff = 1;
sigma = 1;

G = gaussian3d([0 0 0],[sigma sigma sigma]);
% perform blur
[offX,offY,offZ] = ndgrid(-cutoff:cutoff,-cutoff:cutoff,-cutoff:cutoff);
offX = offX(:);
offY = offY(:);
offZ = offZ(:);
offset = [offX offY offZ];
nOffSets = size(offset,1);
nPoints =size(gridIdx,1);
weight = zeros(nOffSets,nPoints);
for i = 1:nOffSets,
     actualOffset = abs(bsxfun(@minus,offset(i,:),remainder));
     weight(i,:) = G(actualOffset(:,1),actualOffset(:,2),actualOffset(:,3));
     isInvalid = sum(actualOffset > cutoff, 2) > 0;
     weight(i,isInvalid) = 0;
end
weight = bsxfun(@rdivide,weight,sum(weight,1));
weight(isnan(weight)) = 0;

% build histogram
ptxMatrix = zeros(histSize);
for i = 1:nOffSets,
    tmpIdx = bsxfun(@plus,gridIdx,offset(i,:));
    isValid = tmpIdx(:,1)>=1 & tmpIdx(:,1)<=histSize(1) & ...
                  tmpIdx(:,2)>=1 & tmpIdx(:,2)<=histSize(2) & ...
                  tmpIdx(:,3)>=1 & tmpIdx(:,3)<=histSize(3);
    ptxMatrix = ptxMatrix + accumarray(tmpIdx(isValid,:),weight(i,isValid),histSize);
end
%change to yxz
ptxMatrix = permute(ptxMatrix,[2 1 3]);
if exist('raw','var')&&raw
    feature = zeros([bincellDim,1036]);
else
    raw =0;
    feature = zeros([bincellDim,size(stick_center,1)+size(count_center,1)]);
end
[ind1,ind2,ind3]=ndgrid(1:size(feature,1),1:size(feature,2),1:size(feature,3));
for i=1:length(ind1(:))
    %cellmatrix =zeros(Qnum,Qnum,Qnum);
    cellmatrix= ptxMatrix((ind1(i)-1)*Qnum+1:min(ind1(i)*Qnum,size(ptxMatrix,1)),...
                          (ind2(i)-1)*Qnum+1:min(ind2(i)*Qnum,size(ptxMatrix,2)),...
                          (ind3(i)-1)*Qnum+1:min(ind3(i)*Qnum,size(ptxMatrix,3)));
    %cellmatrix(1:size(cellmatrix_tmp,1),1:size(cellmatrix_tmp,2),1:size(cellmatrix_tmp,3))= cellmatrix_tmp;
    
    stickfea_raw=cellmatrix(startInd)-cellmatrix(endInd);
    if featurepara.normalize||featurepara.highdim
        stickfea_norm =(stickfea_raw-stick_norm_mean)./stick_norm_scale;
    else
        stickfea_norm = stickfea_raw;
        stdf_stick =28;
    end
   
    
    %count feature 1/2, 1/3 1/6
    countfea1 =shrinkHist(cellmatrix,2);
    countfea2 =shrinkHist(cellmatrix,3);
    countfea3 =sum(cellmatrix(:));
    countfea_raw =[countfea1(:);countfea2(:);countfea3(:)]';
    if featurepara.normalize
        count_norm = (countfea_raw-count_norm_mean)./count_norm_scale;
    else
        count_norm = countfea_raw;
        stdf_count = 90;
    end
    if raw
        feature(ind1(i),ind2(i),ind3(i),1:end)  =[stickfea_norm(:);count_norm(:)];
    else
        stickfea_code = sparse_coding(stickfea_norm,stick_center,stdf_stick,KNN_stick);
        countfea_code =sparse_coding(count_norm,count_center,stdf_count,KNN_count);
        feature(ind1(i),ind2(i),ind3(i),1:end)=[stickfea_code(:);countfea_code(:)]; 
    end
end
end