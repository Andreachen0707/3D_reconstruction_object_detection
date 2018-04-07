function [feature,Space,bincellDim,cellInsight,partailInsight,occcell,missDepthcell,pointCountIntegral,pointCount,occsource,pointRange,ImgplaneSeg] = ...
    mFeature_combined(depthMap,XYZworldframe,s,Rtilt,crop,ImgplaneSeg,featurepara,shift,imageNum)
    if ~exist('shift','var')
        shift = [0,0,0];
    end 
    if ~exist('imageNum','var')
        imageNum = 0;
    end 
    
    KNN =50;
    load('.\featurecenters\noNomalizedcenter\featurecenter_tsdf.mat');
    [pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,Space,removedpts]=init_gridInd(XYZworldframe,featurepara.s,Qnum,shift);
    % pts, xyz is removed nan points 
    [pointCountIntegral,pointCount] =getIntegralPtCount(bincellDim,cellIdx);
    pointRange = getPointRange(bincellDim,cellIdx,xyzShifted);
    [cellInsight,partailInsight] = cellInsightTest(bincellDim,depthMap,Space,Rtilt,crop);
    cellInsight = permute(cellInsight,[2 1 3]);
    partailInsight =permute(partailInsight,[2 1 3]);
    
    display(num2str(bincellDim));
    fprintf('bincellDim = %d',bincellDim(1)*bincellDim(2)*bincellDim(3));
    if bincellDim(1)*bincellDim(2)*bincellDim(3)>80000
       display('oversize');
    end

    tic;
    [feature_point,~]= mFeature_points(gridIdx,remainder,bincellDim,histSize,KNN,KNN,featurepara); 
    feature_raw = feature_point; 
    clear feature_point;
    display('feature_point time:',num2str(toc));

    tic;
    [feature_tsdf,occcell,missDepthcell,occsource] =  mFeature_tsdf(depthMap,histSize,bincellDim,Rtilt,crop,Space,KNN,KNN,featurepara);  

    feature_raw =cat(4,feature_raw,feature_tsdf);
    clear feature_tsdf;
    display('feature_tsdf time:',num2str(toc));


    tic;
    feature_shape= mFeature_shape(bincellDim,cellIdx,xyzShifted,Space,Qnum,KNN,featurepara);
    feature_raw = cat(4,feature_raw,feature_shape);
    clear feature_shape;
    display('feature_shape time:',num2str(toc));
    
    tic;
    [feature_normal,normhistMatrix,normals] = mFeature_normal(xyzShifted,cellIdx,bincellDim,KNN,featurepara); 
    feature_raw = cat(4,feature_raw,feature_normal);
    clear feature_normal;
    display('feature_normal time:',num2str(toc));
    feature = feature_raw;    
end