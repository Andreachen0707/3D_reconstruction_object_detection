 %% training 
if steps.train
   for k=1:length(r(:))
        [Modelfullname] = TrainingMain_test(idxToimageNum,OutputpathTrain, ...
            Rotation(r(k)), Scale(s(k)),Cdepth(d(k)),Cx(x(k)),Tilt(t(k)),Cheight(h(k)),classId,...
            [featurepara.offpath '/' className '/' filename{f(k)} '/' filename{f(k)} '.off'],enlargeBB,enlargewall,...
            enlargefloor(k),featurepara,svmpara,replace,shiftAll(sft(k),:),aspAll(asp(k),:),[]);       
    end
end

fprintf('Training finish\n')
%% testImages using models in testfloder  
if steps.test     
    drawfigure =1;
    Allsvm =dirAll(OutputpathTrain,'.mat');
    for  ii =11,  
        i = idxToimageNum(ii);
        fprintf('testing on Image %d, with %d Models.\n', i,length(Allsvm));
        TestingSingleImage_test(outpathTest,i,Allsvm,featurepara,drawfigure,replace,featurepara.localsearch);
    end  
end

%% gen visulization with AP of this floder
if steps.visul
    modelName = className;
    modelImageFolder=featurepara.modelImageFolder ;
    pointFolder = featurepara.pointpath;
    dataFolder = outpathTest;
    nExamples =60;

    load('groundTruthBbsAll.mat')
    groundTruthBbs = groundTruthBbsAll;
    gtLabels = [groundTruthBbsAll.label]';


    viewScale =0.6;

    calitype ='no';
    testIdRange = idxToimageNum(501:1074);
    trainIdRange = idxToimageNum(1:500);
    
    loadPrevMat =0;
    replacePrevMat=1;
    cali =0;
    prevMatFile= [path2MatRoot ExpirmentName];

    ignoremodelname=1;
    visualizeDetection_load_calibratefunTightBB(prevMatFile,modelImageFolder,OutputpathTrain,pointFolder,dataFolder,testIdRange,trainIdRange,modelName,nExamples,groundTruthBbs,gtLabels,Weboutputpath,viewScale,loadPrevMat,replacePrevMat,cali, calitype,ignoremodelname)
end