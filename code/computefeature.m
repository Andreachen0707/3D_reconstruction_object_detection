function computefeature()
% Example code to compute features on RMRC dataset (provide 'dataset_detection3d_test.mat').
% The features compute using refined depth map are stored in
% '../NYUdatafeature_release/'
initPath;
[svmpara,featurepara]=setParameters();
Outputpath = ['/n/fs/modelnet/NYUdataSet/NYUdatafeature_s' num2str(featurepara.s) '_' featurepara.featype '/']
mkdir(Outputpath);
NYUnew =load('/n/fs/modelnet/NYUdataSet/dataset_detection3d_test.mat','alignmentMatrices','images','rawDepths');

for imageNum =1:length(NYUnew.alignmentMatrices)
    if ~exist([Outputpath 'feature_' num2str(imageNum) '.mat'],'file')
        Rorg =NYUnew.alignmentMatrices{imageNum};
        imgRgb = NYUnew.images(:,:,:,imageNum);
        imgDepth = NYUnew.rawDepths(:,:,imageNum);
        [rgbTest,XYZworldframeTest]=read_3d_pts_NYU(imgRgb, imgDepth, Rorg);
        Rback = decomposeR (Rorg);
        XYZworldframeTest=[Rback*XYZworldframeTest']';    
        Rtmp = Rorg;
        Rtmp = Rtmp(:,[1 3 2]);
        Rtmp = Rtmp([1 3 2],:);
        Rtilt = Rback * Rtmp;
    
        [featureTest,SpaceTest,bincellDim,cellInsight,partailInsight,occcell,missDepthcell,pointCountIntegral,pointCount,occSource,pointRange]...
            = mFeature_combined(imgDepth,XYZworldframeTest,featurepara.s,Rtilt,[],[],featurepara,[0,0,0],imageNum);
        occSource(~occcell) =0;
        save([Outputpath 'feature_' num2str(imageNum) '.mat'],'featureTest','SpaceTest','Rorg','Rtilt','cellInsight','partailInsight','missDepthcell','pointCountIntegral','occSource','pointCount','pointRange','-v6');        
        fprintf('done\n');
    else
        fprintf('skip%s\n', [Outputpath 'feature_' num2str(imageNum) '.mat'])
    end
end
end


