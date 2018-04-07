function feature = getfeatureTypeDim(featureTest,featurepara,cellInsight,partailInsight,missDepthCell,pointCount,occSource)
if featurepara.codesecondlevel,
     F_dim = 1:300;
elseif strcmp(featurepara.featype,'pointtsdfnormalshapeseg'),
    F_dim = 1:350;
elseif strcmp(featurepara.featype,'pointtsdfnormalshape'),
    F_dim = 1:300;
elseif strcmp(featurepara.featype,'point')
    F_dim = 1:100;
elseif strcmp(featurepara.featype,'tsdf')
    F_dim =101:200;
elseif strcmp(featurepara.featype,'tsdffirstorder')
    F_dim =151:200;
elseif strcmp(featurepara.featype,'shape')
     F_dim =201:250;
elseif strcmp(featurepara.featype,'normal')
     F_dim =251:300;
elseif strcmp(featurepara.featype,'pointtsdfnormalshapePlane')||strcmp(featurepara.featype,'pointtsdfnormalshapePlaneCode')||strcmp(featurepara.featype,'pointtsdfnormalshapePlaneCodeSeg')
     F_dim =1:336;     
elseif strcmp(featurepara.featype,'Plane')||strcmp(featurepara.featype,'PlaneCode')||strcmp(featurepara.featype,'PlaneSeg')||strcmp(featurepara.featype,'PlaneCodeSeg')
    F_dim =301:336;
elseif strcmp(featurepara.featype,'pointtsdfnormalshapecurv')
    F_dim = [1:300,337:386];
elseif strcmp(featurepara.featype,'curv')
     F_dim = [337:386];
elseif strcmp(featurepara.featype,'tsdfllc')
     F_dim = [1:700];
elseif strcmp(featurepara.featype,'other')
    F_dim =1:size(featureTest,4);   
elseif strcmp(featurepara.featype,'custom')
    F_dim = [101:200,251:300];
else
    error('wrong feature type')
end

feature =featureTest(:,:,:,F_dim);
if strcmp(featurepara.featype,'pointtsdfnormalshapePlaneCode')||strcmp(featurepara.featype,'PlaneCode')||strcmp(featurepara.featype,'PlaneCodeSeg')||strcmp(featurepara.featype,'pointtsdfnormalshapePlaneCode')||strcmp(featurepara.featype,'pointtsdfnormalshapePlaneCodeSeg')
    fea_code = linear_coding(feature(:,:,:,end-36+1:end),featurepara);
    feature = cat(4,feature(:,:,:,1:end-36),fea_code);
end
totallyout  =~(cellInsight|partailInsight);

if strcmp(featurepara.outtype,'p0tE')&&exist('partailInsight','var'),
   feaIdx = repmat(partailInsight,[1,1,1,size(feature,4)]);
   feature(feaIdx) = 0;
elseif strcmp(featurepara.outtype,'p0t0')&&exist('partailInsight','var'),
   % set cells that not complete in sight to don't care
   feaIdx = repmat(~cellInsight,[1,1,1,size(feature,4)]);
   feature(feaIdx) = 0;
elseif strcmp(featurepara.outtype,'pft0')&&exist('partailInsight','var'),
   % set cells that complete out of sight to don't care
   feaIdx = repmat(totallyout,[1,1,1,size(feature,4)]);
   feature(feaIdx) = 0;
elseif strcmp(featurepara.outtype,'pEtE')&&exist('partailInsight','var'), 
    % set cells that partial out of sight to empty
    [a,b,c] =   find(~(cellInsight|partailInsight));
    if size(a(:)) >1
        k=min(size(a(:),1),30);
        emptycellfea  =featureTest(a(k),b(k),c(k),:);   
        feaIdx = partailInsight;
        [X,Y,Z,D] =size(feature);
        allF = repmat(emptycellfea,[X Y Z 1]);
        feature(repmat(feaIdx,[1 1 1 D])) = allF(repmat(feaIdx,[1 1 1 D]));  
    end
end

% set the missing depth cell to don't care
if strcmp(featurepara.missing,'0')
    feaIdx = repmat(missDepthCell,[1,1,1,size(feature,4)]);
    feature(feaIdx>0) = 0;    
elseif strcmp(featurepara.missing,'E')
     [a,b,c] =   find(~(cellInsight|partailInsight));
     if size(a(:)) >1
        k=min(size(a(:),1),30);
        emptycellfea  =featureTest(a(k),b(k),c(k),:);   
        feaIdx = missDepthCell;
        [X,Y,Z,D] =size(feature);
        allF = repmat(emptycellfea,[X Y Z 1]);
        feature(repmat(feaIdx,[1 1 1 D])) = allF(repmat(feaIdx,[1 1 1 D]));  
     end
elseif strcmp(featurepara.missing,'f')
    %do nothing
end

% scale feature according to their weight 
emptyCellInd = (pointCount<5)&occSource<0.0001;
if isfield(featurepara,'weightempty')&&featurepara.weightempty~=1;
    D =size(feature,4);
    feature(repmat(emptyCellInd,[1 1 1 D])) = featurepara.weightempty*feature(repmat(emptyCellInd,[1 1 1 D]));  
end
occCellInd =occSource>0.0001;
if isfield(featurepara,'weightocc')&&featurepara.weightocc~=1;
    D =size(feature,4);
    feature(repmat(occCellInd,[1 1 1 D])) = featurepara.weightocc *feature(repmat(occCellInd,[1 1 1 D]));  
end


if featurepara.addempty>0
    feature =cat(4,feature,2*featurepara.addempty*(pointCount<5)-featurepara.addempty);
end

if featurepara.selfocc>0
   feature =cat(4,feature,2*featurepara.selfocc*(occSource>0.0001)-featurepara.selfocc);
end

if featurepara.outadd&&featurepara.missingadd
    if isfield(featurepara,'scalemissing')
        feature =cat(4,feature,featurepara.scalemissing*(~cellInsight|missDepthCell));
    else   
        feature =cat(4,feature,~cellInsight|missDepthCell);
    end
else
    if featurepara.outadd
        feature =cat(4,feature,~cellInsight);
    end
    if  featurepara.missingadd
        feature =cat(4,feature,missDepthCell);
    end
end


