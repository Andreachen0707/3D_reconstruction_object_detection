function pointRange = getPointRange(bincellDim,cellIdx,xyzShifted)
        pointRange = nan([bincellDim,6]);
        %pointCount2 = zeros(size(pointCount));
        linearIndPts =sub2ind(bincellDim,cellIdx(:,2),cellIdx(:,1),cellIdx(:,3));
        unique_linearIndPts = unique(linearIndPts);
        for i =1:length(unique_linearIndPts)
            linearIDXPick=linearIndPts == unique_linearIndPts(i);
            cellIdxCurr = cellIdx(find(linearIDXPick,1),:);
            pointRange(cellIdxCurr(2),cellIdxCurr(1),cellIdxCurr(3),:) =[min(xyzShifted(linearIDXPick,:),[],1),max(xyzShifted(linearIDXPick,:),[],1)];
            %pointCount2(cellIdxCurr(2),cellIdxCurr(1),cellIdxCurr(3),:) =sum(linearIndPts==unique_linearIndPts(i));
        end
end
