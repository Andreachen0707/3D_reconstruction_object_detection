function [pointCountIntegral,pointCount] =getIntegralPtCount(bincellDim,cellIdx)
pointCount = zeros(bincellDim);
[ind1,ind2,ind3]=ndgrid(1:size(pointCount,1),1:size(pointCount,2),1:size(pointCount,3));
for i=1:length(ind1(:))
    ptsIncellidx = cellIdx(:,2)==ind1(i)&cellIdx(:,1)== ind2(i)&cellIdx(:,3)==ind3(i);
    pointCount(ind1(i),ind2(i),ind3(i),1:end)= sum(ptsIncellidx(:));
end
pointCountIntegral = cumsum(cumsum(cumsum(pointCount,1),2),3);