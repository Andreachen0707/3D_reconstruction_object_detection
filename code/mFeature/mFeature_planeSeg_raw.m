function feature = mFeature_planeSeg_raw(xyzShifted,removedpts,cellIdx,bincellDim,imageSize,planeSeg,featurepara)
numbin = 36;
angle = 0:2*pi/numbin:2*pi;
angle =angle(1:numbin);
disThreshold = Inf;
if isfield(featurepara,'removefloor')&&featurepara.removefloor,
   planeSeg(planeSeg==1)=0;
end
feature = zeros([bincellDim,numbin]);
[ind1,ind2,ind3]=ndgrid(1:size(feature,1),1:size(feature,2),1:size(feature,3));
[I,J] = ndgrid(1:imageSize(1),1:imageSize(2));
I =I(~removedpts);
J =J(~removedpts);
planeSeg = planeSeg(~removedpts);
% get the surface normal for each plane 
unique_planeSeg =unique(planeSeg);
for i =1:length(unique_planeSeg)
    if unique_planeSeg(i)>0
        onPlaneIdx = planeSeg==unique_planeSeg(i);
        ptsOnConectedplane = xyzShifted(onPlaneIdx,:);
        if size(ptsOnConectedplane,1)>10,
            [C,V]= PCA(ptsOnConectedplane',3);
            planeNormal = C(:,3);
            if [0,1,0]*planeNormal>0,
               planeNormal =-1*planeNormal;
            end
            planeNormal =planeNormal/norm(planeNormal);
            planeNormal_allplane(:,unique_planeSeg(i)) = planeNormal;
        else
            planeSeg(planeSeg==unique_planeSeg(i)) =0;
        end
    end
   
end
clear onPlaneIdx planeNormal ptsOnConectedplane;
for i=1:length(ind1(:))
    planeShape = zeros(1,numbin);
    ptsIncellidx = find(cellIdx(:,2)==ind1(i)&cellIdx(:,1)== ind2(i)&cellIdx(:,3)==ind3(i));

    % get planId 
    planIdIncell = planeSeg(ptsIncellidx);
    planIdIncell = planIdIncell(planIdIncell>0);
    numOfsegIncell = length(unique(planIdIncell));
    if length(planIdIncell)>10
         planIdIncell_pick = mode(planIdIncell(planIdIncell>0));
         % get points on plane
         onPlaneIdx = planeSeg==planIdIncell_pick;
         % get plane function 
         ptsIncell = xyzShifted(ptsIncellidx,:);
         cell_center = nanmean(ptsIncell(planIdIncell == planIdIncell_pick,:),1); 
         ptsOnConectedplane = xyzShifted(onPlaneIdx,:);
         if ~isempty(ptsOnConectedplane)
            planeNormal =planeNormal_allplane(:,planIdIncell_pick);
            % choose the plane to intersect 
             spaceface =[0,0,1;0 1 0];
             [~,facechoose] = min(abs(spaceface*planeNormal));
             startLineA = cross(planeNormal,spaceface(facechoose,:));
             startLineA = startLineA/norm(startLineA);
             startLineB = cross(planeNormal,startLineA);
             % 
             if facechoose ==1,
                 % z up
                 startLineA = startLineB(3)/abs(startLineB(3))*startLineA;
                 startLineB = startLineB(3)/abs(startLineB(3))*startLineB;
             else
                 % x right   
                 startLineA = startLineA(1)/abs(startLineA(1))*startLineA;
                 startLineB = startLineA(1)/abs(startLineA(1))*startLineB;
             end

             if sum(isnan([startLineA startLineB]))>0,
                 error('wrong normal')
             end

             lineDir = [startLineA(:),startLineB(:),planeNormal(:)]*[cos(angle);sin(angle);zeros(size(angle))];
             projectOnLine = bsxfun(@minus, ptsOnConectedplane, cell_center)*lineDir;
             [~,chooseBin] = max(projectOnLine,[],2);
             projectMatrix = zeros(size(projectOnLine));
             assignInd = sub2ind(size(projectMatrix),1:size(projectMatrix,1),chooseBin');
             projectMatrix(assignInd) = projectOnLine(assignInd);
             planeShape = [max(projectMatrix,[],1)];
             planeShape(planeShape>disThreshold) =disThreshold;
         end
    end
    %planeShape = conv([planeShape(end) planeShape(:)' planeShape(1)],gausswin(3)/sum(gausswin(3)),'valid');
    planefea_raw =planeShape(:)+0.01;
    feature(ind1(i),ind2(i),ind3(i),1:end)  = planefea_raw;
end
end