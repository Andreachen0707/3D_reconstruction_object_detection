function hardnegativeMining(id)
% /n/fs/vision/ionic/starter.sh hardnegativeMining 7200mb 24:00:00 1 51 1 /n/fs/modelnet/log/
% f =dir('/n/fs//modelnet/SUN3DFeature/Notoilet/feature/');
% for i =3:length(f)
%     if ~exist(['/n/fs//modelnet/SUN3DFeature/Notoilet/images/' sprintf('%015d',str2double(f(i).name(9:end-4))) '.jpg'],'file')
%        delete(['/n/fs//modelnet/SUN3DFeature/Notoilet/feature/' f(i).name])
%     end
% end

ConstsDetection;
classNames = consts.classNames;
className = 'chair';
[~,classId] =ismember(className,classNames);
[idxToimageNum] = split_trainingtestingSetGroup(classId);
initPath;
dbstop if error;
OldModelpath ='/n/fs/modelnet/TrainedModels/chair_Jan22/chair/chair_000000358/';%'/n/fs/modelnet/TrainedModels/toilet_addwall2/';
savesvPath = '/n/fs/modelnet/TrainedModels/chair_Jan22_SV/';
NewOutputpath = [ '/n/fs/modelnet/TrainedModels/chair_Jan22_miningOnGC/'];
postfix = '';
pathToNegs = '/n/fs/modelnet/CGNeg/feature/';%['/n/fs/modelnet/SUN3DFeature/No' className '/']
AllModels = dirAll(OldModelpath,'.mat');



debug.rawdetectionPath ='/n/fs/modelnet/hardnegMining/Training/';
replace =1;
[svmpara,featurepara]=setParameters();

if ~exist('id','var')
    id =1:length(AllModels);
end
for i =id
    debug.showimage =mod(i,5)==1;
    load(AllModels{i});
    Modelfullname =[Model.name '_r' num2str(Model.Rotation) '_s' num2str(Model.Scale) '_d' num2str(Model.Cdepth) '_t' num2str(Model.Tilt) '_x' num2str(Model.Cx) '_h' num2str(Model.Cheight)];
    SVfullpath = [savesvPath classNames{Model.classId} '/' Model.name '/' Modelfullname '_SV.mat'];
    
    try load(SVfullpath),catch, fprintf(['no sv:' AllModels{i} '\n']);continue;end
    
    hardnegativeMiningSingleModel(idxToimageNum,NewOutputpath,pathToNegs,Model,nSV,pSV,svmpara,featurepara,debug,replace); 
    
    if debug.showimage,
        %%gen result web page 
        Modelfullname =[Model.name '_r' num2str(Model.Rotation) '_s' num2str(Model.Scale) '_d' num2str(Model.Cdepth) '_t' num2str(Model.Tilt) '_x' num2str(Model.Cx) '_h' num2str(Model.Cheight)];
        modelName = Modelfullname;
        modelImageFolder='/n/fs/modelnet/Code/Output/WebResult/modelImageFolder/';
        %pointFolder = '/n/fs/modelnet/SUN3DFeature/Notoilet/feature/';
        pointFolder = '/n/fs//modelnet/CGNeg/datadepth/';
        dataFolder = debug.rawdetectionPath;
        otherNegfiles = dirAll(pathToNegs,'.mat');
        nExamples =50;
        viewScale =0.65;
        testIdRange =[1,Inf];
        iFea =0;
        outputpath =['/n/fs/modelnet/Code/Output/WebResult/hardnegMining/' Modelfullname '/'];
        visualizeDetection_noGT(modelImageFolder,pointFolder,dataFolder,testIdRange,modelName,nExamples,outputpath,viewScale,otherNegfiles); 
        %}
    end
end
%}

rootPath = NewOutputpath;
outpathTest =['/n/fs/modelnet/hardnegMining/Testing/' className postfix '/'];
% Get all Model files
folders = regexp(genpath(rootPath), pathsep, 'split');
folders = folders(1:end-1);
cnt = 0;
for f=1:length(folders)
    files = dir([folders{f} '/*.mat']);
    for j=1:length(files)
        if ~files(j).isdir
            cnt = cnt + 1;
            Allsvm{cnt} = fullfile(folders{f},files(j).name);
        end
    end
end
%while 0
    if 1%length(Allsvm) ==length(AllModels),
    drawfigure =0;
    replace =0;
    [idxToimageNum] = split_trainingtestingSetGroup(classId);

    parfor ii =1:1074,  
         i = idxToimageNum(ii);
         fprintf('testing on Image %d, with %d Models.\n', i,length(Allsvm));
         TestingSingleImage_type(outpathTest,i,Allsvm,featurepara,drawfigure,replace,1);
    end
%}
    % show result
    modelName = Modelfullname;
    modelImageFolder='/n/fs/modelnet/Code/Output/WebResult/modelImageFolder/';
    pointFolder = '/n/fs/modelnet/NYUdataSet/NYUdatafeatureNew/';
    dataFolder = outpathTest;
    nExamples =70;
    viewScale =0.65;
    load('/n/fs//modelnet/NYUdataSet/NewGt/groundtruthBBYawNew.mat')
    groundTruthBbs = groundtruthBBYawNew;
    gtLabels = [groundtruthBBYawNew.label]';

    testIdRange =idxToimageNum(501:1074);
    trainIdRange =idxToimageNum(1:500);
    cali =1; calitype ='no';  %testingCali %trainingCaliLocal %no
    outputpath =['/n/fs/modelnet/Code/Output/WebResult/AfterMining/'  modelName postfix '/'];
    OutputpathTrain = NewOutputpath;
    
    visualizeDetection_load_calibratefunTightBB([],modelImageFolder,OutputpathTrain,pointFolder,dataFolder,testIdRange,trainIdRange,modelName,nExamples,groundTruthBbs,gtLabels,outputpath,viewScale,1,0,cali, calitype)
    %break;
    end
%end
end



function Model = hardnegativeMiningSingleModel(idxToimageNum,Outputpath,pathToNegs,Model,nSV,pSV,svmpara,featurepara,debug,replace)
ConstsDetection;
classNames = consts.classNames;
folder =[Outputpath classNames{Model.classId} '/' Model.name '/'];
Modelfullname =[Model.name '_r' num2str(Model.Rotation) '_s' num2str(Model.Scale) '_d' num2str(Model.Cdepth) '_t' num2str(Model.Tilt) '_x' num2str(Model.Cx) '_h' num2str(Model.Cheight)];
Modeltosavefullpath =[folder  Modelfullname '.mat'];
if exist(Modeltosavefullpath,'file')&&~replace, fprintf('skip%s\n',Modeltosavefullpath);
    return;
end

% get all .mat file 
NumOfNegUnflip =500;
otherNegfiles = dirAll(pathToNegs,'.mat');
NumOfNeg =NumOfNegUnflip + length(otherNegfiles);
fprintf('training on %d files\n ',NumOfNeg);


iterInfor = [Model.iter; nan(1,size(Model.iter,2))];
hardnegthr = Model.hardnegthr;
size_pos =Model.size_pos;
hardnegthr_low = -0.95;
size_pos =Model.size_pos;
load('/n/fs/modelnet/Code/detector_change_feature/groundTruthBbsYaw.mat','groundTruthBbsYaw');
allgtnum = sum(ismember([groundTruthBbsYaw.imageNum],idxToimageNum(1:NumOfNegUnflip))&[groundTruthBbsYaw.classId] ==Model.classId);
for iter=1:svmpara.maxiter,
    hardnegstart =tic;
    NegAll =[];
    NexscoreAll =[];
    osall =[];
    oscfall =[];
    scoreall=[];
   if iter ==1, start =NumOfNegUnflip+1;else start =1;end
   for imageNumid=start:NumOfNeg
        if size(NegAll,2)<svmpara.NumOfNexLimit       
            if imageNumid<NumOfNegUnflip+1,
            imageNum = idxToimageNum(imageNumid);
            load([featurepara.featurepath 'feature_' num2str(imageNum) '.mat']);
            if svmpara.exculdepos
                postive_bb =groundTruthBbsYaw([groundTruthBbsYaw.imageNum]==imageNum&[groundTruthBbsYaw.classId] ==Model.classId);
            else
                postive_bb =[];
            end
        else
            %% find hard negtives in this image
            imageNum = imageNumid;
            load(otherNegfiles{imageNumid-NumOfNegUnflip},'SpaceTest','featureTest','pointCountIntegral','cellInsight','missDepthcell','partailInsight','Rtilt')
            if ~exist('pointCount','var'),pointCount =[];end
            postive_bb =[];
        end
        featureNeg = getfeatureType(featureTest,featurepara,cellInsight,partailInsight,missDepthcell);
        SpaceNeg =SpaceTest;
        fprintf('compute hard neg on image %d, iter %d Model: %s \n',imageNumid,iter,Modelfullname);
       
       % not handle occ, no local search, no pading  on training
        pad =0;
        [bbw_3d_d,bbs_3d_f] =detection3D_Range_occ(featureNeg,svmpara.numbox,SpaceNeg,Model,pointCountIntegral,zeros(size(pointCountIntegral)),zeros(Model.size_pos([1,2,3])),0,pad,featurepara.removedontcare);       
        
        if ~isempty(bbw_3d_d)         
                % nms
                indexes  = nmsMe_3d([bbw_3d_d(:,1:3),bbw_3d_d(:,1:3)+bbw_3d_d(:,4:6),bbw_3d_d(:,7)],0.3); 
                bbs_3d_f = bbs_3d_f(indexes,:);
                bbw_3d_d = bbw_3d_d(indexes,:);
                
                if debug.showimage&&iter==1,
                    todsavefloderpath2 =sprintf('%s/%s/%s/%s/',debug.rawdetectionPath,classNames{Model.classId},Model.name,Modelfullname);                    
                    if ~exist(todsavefloderpath2,'dir'), mkdir(todsavefloderpath2);end
                    tosavename2=sprintf('%s/Test%05d-%s.mat',todsavefloderpath2,imageNum,Modelfullname);
                    bbtight3d = transform2tightBB(bbw_3d_d,Model.BBTight);
                    [bb_2d_d,bb2dDraw] = projectStructBbsTo2d(bbtight3d,Rtilt);
                    result=struct('bb3d',bbw_3d_d,'bb2d',bb_2d_d,'bb2dDraw',bb2dDraw,'bbw_3d_f',bbs_3d_f,'bbtight3d',bbtight3d);
                    %result=struct('bb3d',bbw_3d_d,'bb2dDraw',bb2dDraw);
                    save(tosavename2,'result','-v6')
                end
                if ~isempty(postive_bb)
                   neg_bb_w =bbf2w(bbs_3d_f,SpaceNeg);
                   os=bb3dOverlapApprox(neg_bb_w,postive_bb);
                   oscf=bb3dOverlapApprox(neg_bb_w,postive_bb);%bb3dOverlapCloseForm(neg_bb_w,postive_bb);
                   % one gt find one match
                   [maxgtos,matchind] =max(os,[],1);
                   selectos =zeros(size(os));
                   selectos(sub2ind(size(selectos),matchind,1:size(os,2))) = maxgtos;                  
                   maxos =max(selectos,[],2);
                   idxfar = maxos<0.00001; 

                   %repeat for oscf
                   [maxgtoscf,matchindcf] =max(oscf,[],1);
                   selectoscf =zeros(size(oscf));
                   selectoscf(sub2ind(size(selectoscf),matchindcf,1:size(selectoscf,2))) = maxgtoscf;
                   maxoscf =max(selectoscf,[],2);

                   osall =[osall;maxos];
                   oscfall =[oscfall;maxoscf];
                   scoreall =[scoreall;bbs_3d_f(:,7)];
                   bbs_3d_f = bbs_3d_f(idxfar,:);
                   bbw_3d_d = bbw_3d_d(idxfar,:);
                else
                    osall =[osall;zeros(size(bbs_3d_f(:,7)))];
                    oscfall = [oscfall;zeros(size(bbs_3d_f(:,7)))];
                    scoreall =[scoreall;bbs_3d_f(:,7)];                       
                end

                % hard neg 
                if ~isempty(bbs_3d_f)
                    hard_neg_id = bbw_3d_d(:,7) > hardnegthr_low;
                    bbw_3d_d = bbw_3d_d(hard_neg_id,:);
                    bbs_3d_f = bbs_3d_f(hard_neg_id,:);
                end

                % get the neg features 
                Nex = zeros(prod(size_pos),size(bbs_3d_f,1));
                Nexscore = bbs_3d_f(:,7)';
                for i=1:size(bbs_3d_f,1)
                    nex = getTemplate_3d(featureNeg,bbs_3d_f(i,1:6));
                    Nex(:,i) = nex(:);
                end
        else
                Nex =[];
                Nexscore =[];
        end
            NegAll =[NegAll,Nex];
            NexscoreAll =[NexscoreAll,Nexscore];
            fprintf('New hard neigative: %d Totol: %d threshold: %d\n',size(Nex,2),size(NegAll,2),hardnegthr);                        
        else
            fprintf('NegAll exceed max number !\n')
            break;
        end
    end

    thardneg=toc(hardnegstart);
    fprintf('Time to collect hard negs :%f\n',thardneg);
    hardnegthr =update_hardnegthr(osall,scoreall,hardnegthr,svmpara.socrethre,allgtnum);
    if ~isempty(NegAll) 
        NegAll =NegAll(:,NexscoreAll>hardnegthr);
    end
    numNewNex = size(NegAll,2);
    iterInfor =[iterInfor;iter,size(pSV,1),size(nSV,1),hardnegthr,NaN,NaN,numNewNex,imageNumid];

    %% converg or retrain
    if isempty(NegAll)&&imageNumid==NumOfNeg,
        if iter>1,
            fprintf('NegAll are empty !\n');
            break;
        else
            continue;
        end
    else
        % retain SVM
        fprintf('SVM %d retaining start !\n NegAll size: %d \n', iter, size(NegAll,2));
        tic;
        
        if featurepara.removedontcare
            careInd = repmat(Model.careMask,[1,1,1,size_pos(4)]);                      
            NegAll = NegAll(find(careInd(:)),:);
        else           
            careInd =ones(size_pos);
        end
        
        Nex =[NegAll';nSV];
        Pos =[pSV];
        x = [Pos;Nex];
        y = [ones(size(Pos,1),1); -1*ones(size(Nex,1),1)];
        notnanInd=~isnan(sum(x,2));
        x=x(notnanInd,:);
        y=y(notnanInd,:);


        w1 = svmpara.w1;
        clear Pos Nex NegAll featureNeg bbs_3d_f ;
        model = svmtrain(y, x, sprintf('-b 0 -s 0 -t 0 -c %.20f -w1 %.9f',svmpara.c, w1));
        
        w = model.sv_coef' * model.SVs;
        bias = -model.rho;
        template =zeros(size_pos);
        template(careInd) = w;
        
        pSV =model.SVs(find(model.sv_coef>0),:);
        %nSV =[model.SVs(find(model.sv_coef<0),:)];
        nSV =[nSV;model.SVs(find(model.sv_coef<0),:)];
        nSV =unique(nSV,'rows');
        svm =struct('w',template,'bias',bias);

        Model.hardnegthr =hardnegthr;
        Model.svm =svm;
        if numNewNex<5&&imageNumid==NumOfNeg,
            break;
        end
        t=toc;
        fprintf('SVM %d retaining end !\n Time: %f \n', iter,t);
    end
end

Model.iter = iterInfor;
mkdir(folder);
save(Modeltosavefullpath,'Model','-v6');
end