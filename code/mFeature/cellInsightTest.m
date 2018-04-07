function [cellInsight,partailInsight] = cellInsightTest(bincellDim,depthMap,Space,Rtilt,crop)
% cell insight test
[X,Y,Z] = ndgrid(1:bincellDim(2)+1,1:bincellDim(1)+1,1:bincellDim(3)+1);

cellCoord = bsxfun(@times,[X(:) Y(:) Z(:)]-1,[Space.s Space.s Space.s]);
cellCoord = bsxfun(@plus,cellCoord,[Space.Rx(1),Space.Ry(1),Space.Rz(1)]);

[cellProj,cellProjDepth] = project3dPtsTo2d(cellCoord,Rtilt,crop);   
cellProj = round(cellProj);
isValid = cellProj(:,1) >= 1 & cellProj(:,1) <= size(depthMap,2) & ...
          cellProj(:,2) >= 1 & cellProj(:,2) <= size(depthMap,1) & ...
          cellProjDepth > 0;
%reshape
Insight =reshape(isValid,size(X));      
cellInsight = Insight(1:end-1,1:end-1,1:end-1)&Insight(2:end,1:end-1,1:end-1)&Insight(1:end-1,2:end,1:end-1)&Insight(1:end-1,1:end-1,2:end);
partailInsight = ~cellInsight&(Insight(1:end-1,1:end-1,1:end-1)|Insight(2:end,1:end-1,1:end-1)|Insight(1:end-1,2:end,1:end-1)|Insight(1:end-1,1:end-1,2:end));

end