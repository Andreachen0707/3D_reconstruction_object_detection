function shapeContext = getShapeContext(cell_center,ptsofThisSeg,distance_threshold,disBinNum,directionsBin)
        %linearSep = (log(distance_threshold(2))-log(distance_threshold(1)))/disBinNum;
        linearSep = (distance_threshold(2)-distance_threshold(1))/disBinNum;
        vector  = bsxfun(@minus,ptsofThisSeg,cell_center);
        distance = sqrt(sum(vector.*vector,2));
        direction = bsxfun(@rdivide, vector,distance);
        farind = distance> distance_threshold(2)|distance< distance_threshold(1);
        
        vector(farind,:) =[];
        distance(farind,:) =[];
        direction(farind,:) =[];
        [~, angleIdx] = max(directionsBin * direction',[],1);
        %disIdx = max(min(floor((log(distance) -log(distance_threshold(1)))/linearSep)+1,disBinNum),1);
        disIdx = max(min(floor((distance -distance_threshold(1))/linearSep)+1,disBinNum),1);
        shapeContext =zeros(size(directionsBin,1),disBinNum);
        tmpIdx =[angleIdx(:),disIdx];
        weight =ones(size(tmpIdx,1),1);
        shapeContext =accumarray(tmpIdx,weight,size(shapeContext));
        % normalize by binVolume
        Volume = (([1:disBinNum]*linearSep) + distance_threshold(1)).^2;
        binVolume = Volume - [0,Volume(1:end-1)];
        shapeContext =bsxfun(@rdivide,shapeContext,binVolume);
end