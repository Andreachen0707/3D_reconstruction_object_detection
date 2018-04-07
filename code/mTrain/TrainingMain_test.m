function [Modelfullname,Modeltosavefullpath,folder]= TrainingMain_test(idxToimageNum,Outputpath,Rotation,Scale,Cdepth,Cx,Tilt,Cheight,...
                      classId,offfile,enlargeBB,enlargewall,enlargefloor,featurepara,svmpara,replace,shift,asp, savesvPath)
ConstsDetection;
classNames = consts.classNames;
tic;
ind =find(offfile=='/');
modelname = offfile(ind(end)+1:end);
ind =find(modelname=='.');
modelname =modelname(1:ind(end)-1);
Model.Rotation =Rotation;
Model.Scale =Scale;
Model.classId =classId;
Model.Cdepth =Cdepth;
Model.Cheight =Cheight;
Model.Tilt =Tilt;
Model.enlarge =enlargeBB;
Model.name =modelname;
Model.Cx =Cx;
Model.featurepara = featurepara;

folder =[Outputpath classNames{Model.classId} '/' Model.name '/'];
Modelfullname =[Model.name '_r' num2str(Model.Rotation) '_s' num2str(Model.Scale) '_d' num2str(Model.Cdepth) '_t' num2str(Model.Tilt) ...
                '_x' num2str(Model.Cx) '_h' num2str(Model.Cheight) '_sft' num2str(shift(1)) num2str(shift(2)) num2str(shift(3))...
                '_asp' num2str(asp(1)) num2str(asp(2)) num2str(asp(3)) ];
Model.Modelfullname = Modelfullname;
Modeltosavefullpath =[folder  Modelfullname '.mat'];


if exist(Modeltosavefullpath,'file')&& ~replace
    fprintf('skip training model: %s\n.',Modeltosavefullpath);  
else
%% gen Pos
%  rendering Model depth  map  (offfile,ratio, xzRot,tilt, X(right),Y(depth_in),Z(up),modelsize,addfloor,foorenlarge,enlargewall)
[depth,K,crop,bb3d,vmat,meshObj,BBTight,segMask] = off2im_sys(offfile,featurepara.camera.ratio,Model.Rotation*pi,Model.Tilt,Model.Cx,Model.Cdepth,Cheight,[1;1;1]*Model.Scale,1,enlargefloor,enlargewall);
[rgb,XYZworldframe] = read_3d_pts_general(depth,K,size(depth),[],crop);
if size(XYZworldframe,1) <50,
   return;
end
Rtilt =  getRotationMatrix('tilt',Model.Tilt);
XYZworldframe = (Rtilt*XYZworldframe')';
s =0.1;
% caculate feature for Model
% [feature,Space,bincellDim,cellInsight,partailInsight,occcell,missDepthcell,pointCountIntegral,pointCount]
if featurepara.highdimcombine
    featureparacombine = featurepara;
    featureparacombine.highdim =0;
    featureparacombine.codesecondlevel=0;
    [feature1,SpacePre,bincellDim,cellInsight,partailInsight,occcell,missDepthcell,pointCountIntegral,pointCount,occsource] = mFeature_combined(depth,XYZworldframe,s,Rtilt,[1,1],segMask,featurepara,shift);
    featureparacombine.highdim =1;
    featureparacombine.codesecondlevel=1;
    [feature2,SpacePre,bincellDim,cellInsight,partailInsight,occcell,missDepthcell,pointCountIntegral,pointCount,occsource,pointRange,ImgplaneSeg] = mFeature_combined(depth,XYZworldframe,s,Rtilt,[1,1],segMask,featurepara,shift);
    feature =cat(4,feature1,feature2);
else
    [feature,SpacePre,bincellDim,cellInsight,partailInsight,occcell,missDepthcell,pointCountIntegral,pointCount,occsource,pointRange,ImgplaneSeg] = mFeature_combined(depth,XYZworldframe,s,Rtilt,[1,1],segMask,featurepara,shift);
end
[Space,XYZworldframe, bb3dfea,Model.enlarge]= initModelBB(XYZworldframe,s,enlargeBB,bb3d,SpacePre);

feature = getTemplate_3d(feature,bb3dfea);
Model.occcell = getTemplate_3d(occcell,bb3dfea);
partailInsight = getTemplate_3d(partailInsight,bb3dfea);
cellInsight = getTemplate_3d(cellInsight,bb3dfea);

Model.pointCount = getTemplate_3d(pointCount,bb3dfea);
Model.bb3dfea =bb3dfea;
Model.BBTight =BBTight;
[Model.careMask,~] =getDontcareMaskTight(size(feature),Space,meshObj,Model.pointCount);
Model.Seg = ImgplaneSeg;
Model.feature = feature;
Model.XYZworldframe =XYZworldframe;
feature = getfeatureTypeDim(feature,featurepara,cellInsight,partailInsight,zeros(size(feature,1),size(feature,2),size(feature,3),1),Model.pointCount,Model.occcell);

pos = feature;
size_pos =size(pos);
bb_3d_f =bb3dfea;
Pos_all = pos(:);


postivefloder =[featurepara.modelImageFolder classNames{Model.classId} '/' modelname '/'];
filename=[ postivefloder Modelfullname '.jpg'];
filename2=[postivefloder Modelfullname '.mat'];

if ~exist(postivefloder,'dir'),mkdir(postivefloder);end
randerim = getImagesc(depth);
imwrite(randerim,filename);
save(filename2,'XYZworldframe');
clear randerim feature featureR1 featureR2 featureR3 featureR4 featureR5 featureR6;

t=toc;
fprintf('finish gen postives time : %f\n',t); 
fprintf('Start training model: \n %s\n.',Modeltosavefullpath);  
tic
% random gen Neg and train SVM
load('groundTruthBbsYaw.mat','groundTruthBbsYaw');

NumOfNeg =5;
randNegAll =[];
idxToimageNum(8)=8;
idxToimageNum(9)=9;
idxToimageNum(10) =10;
idxToimageNum(11)=11;
idxToimageNum(12)= 12;
idxToimageNum(13)=13;
idxToimageNum(14)=14;
for ii =1:NumOfNeg
    imageNum = idxToimageNum(ii);
    %filename4=[ featurepara.featurepath 'feature_' num2str(imageNum) '.mat'];
    load([ featurepara.featurepath 'feature_' num2str(imageNum) '.mat']); %ÔØÈëÁËÍ¼Æ¬
    SpaceNeg= SpaceTest;
    featureNeg = getfeatureTypeDim(featureTest,featurepara,cellInsight,partailInsight,missDepthcell,pointCount,occSource);
    if svmpara.exculdepos
        postive_bb = groundTruthBbsYaw([groundTruthBbsYaw.imageNum]==imageNum&[groundTruthBbsYaw.classId] ==classId);
    else
        postive_bb =[];
    end

    if isempty(find((size(featureNeg)-size(pos))<0))
       %random generate negs 
       [randNeg,neg_bb_w] =randgenNeg(size(pos),featureNeg,bb_3d_f,postive_bb,SpaceNeg,round(svmpara.NumOfNexLimit/NumOfNeg));
       randNegAll =[randNegAll,randNeg];
    end
end
t=toc;
NumNex =size(randNegAll,2);
fprintf('Random gen Negs take: %f\n NegSVM first taining start !\n NegAll size: %d \n', t,NumNex);

% remove don't care regoin
if featurepara.removedontcare
    careInd = repmat(Model.careMask,[1,1,1,size_pos(4)]);
    randNegAll = randNegAll(find(careInd(:)),:);
    Pos_all = Pos_all(find(careInd(:)),:);
else
    careInd =ones(size_pos);
end
Pos =Pos_all';
Nex =randNegAll';
x = [Pos;Nex];
y = [ones(size(Pos,1),1); -1*ones(size(Nex,1),1)];
notnanInd=~isnan(sum(x,2));
x=x(notnanInd,:);
y=y(notnanInd,:);
pSV = Pos;
tic;
clear Pos Nex randNegAll pts;
if svmpara.linear 
    model = train(y,sparse(x), sprintf('-s 2 -B %f -c %f -w1 %f -q', 67500,svmpara.c, svmpara.w1));
    w = model.w(1:end-1);
    bias = model.w(end)*67500;
    nSV =x(find(y(2:end).*(x(2:end,:)*w(:)+bias)<1)+1,:);
    nSV =unique(nSV,'rows');
else
    model = svmtrain(y, x, sprintf(['-b 0 -s 0 -t 0 -c %.50f -w1 %.9f'],svmpara.c,svmpara.w1));
    w = model.sv_coef' * model.SVs;
    bias = -model.rho;
    pSV=model.SVs(find(model.sv_coef>0),:);
    nSV=model.SVs(find(model.sv_coef<0),:);
end

template =zeros(size_pos);
template(careInd>0) = w;
svm =struct('w',template,'bias',bias);
iterInfor =[0,size(pSV,1),size(nSV,1),NaN,NaN,NaN,NumNex,ii]; 
Model.svm =svm;
Model.size_pos =size(pos);
clear x y model pos;

t= toc; fprintf('SVM first training ends, takes %f\n',t);
%% add empty box to negative SV
[ei,ej,ek]=find(pointCount<1&(occSource==0),1);
feaempty = featureNeg(ei(1),ej(1),ek(1),:);
emSV = repmat(feaempty,[size_pos([1,2,3]),1]);
emSV = emSV(:);
emSV = emSV(find(careInd(:)),:);
nSV =[nSV;emSV(:)'];
%}
%% Hard negative mining 
NumOfNegUnflip =5;  % change the number of negtive file number
NumOfNeg =NumOfNegUnflip; 
NYUNegfiles =[];otherNegfiles =[];
if isfield(svmpara,'NYUNeg')&&svmpara.NYUNeg,
    NYUNegfiles = dirAll(featurepara.featurepath_NYUneg, '.mat');
    NumOfNeg =NumOfNeg + length(NYUNegfiles);
end   
if isfield(svmpara,'otherNeg')&&svmpara.otherNeg,
    % train on other dataset 
    otherNegfiles = dirAll(featurepara.featurepath_moreneg, '.mat');
    NumOfNeg =NumOfNeg + length(otherNegfiles);
end
fprintf('training on %d files\n ',NumOfNeg);   
hardnegthr_low = -0.8;
hardnegthr =hardnegthr_low;
size_pos =Model.size_pos;
osall =[];
scoreall=[];
load('groundTruthBbsYaw.mat','groundTruthBbsYaw');
allgtnum = sum(ismember([groundTruthBbsYaw.imageNum],idxToimageNum(1:NumOfNegUnflip))&[groundTruthBbsYaw.classId] ==classId);
for iter=1:svmpara.maxiter,
    hardnegstart =tic;
    NegAll =[];
    NexscoreAll =[];
    osall =[];
    oscfall =[];
    scoreall=[];
    for imageNumid=1:NumOfNeg
        if size(NegAll,2)<svmpara.NumOfNexLimit
        tic;        
        %% find hard negtives in this image
        if imageNumid<NumOfNegUnflip+1,
            imageNum = idxToimageNum(imageNumid);
            load([featurepara.featurepath 'feature_' num2str(imageNum) '.mat']);
            if svmpara.exculdepos
                postive_bb =groundTruthBbsYaw([groundTruthBbsYaw.imageNum]==imageNum&[groundTruthBbsYaw.classId] ==classId);
            else
                postive_bb =[];
            end
        elseif imageNumid>NumOfNegUnflip&&imageNumid<length(NYUNegfiles)+NumOfNegUnflip+1
            load(NYUNegfiles{imageNumid-NumOfNegUnflip});
            postive_bb =[];
        else
             load(otherNegfiles{imageNumid-length(NYUNegfiles)-NumOfNegUnflip});
             if ~exist('pointCount','var'),pointCount =[];end
             postive_bb =[];
        end
        featureNeg = getfeatureTypeDim(featureTest,featurepara,cellInsight,partailInsight,missDepthcell,pointCount,occSource);
        SpaceNeg =SpaceTest;
        fprintf('compute hard neg on image %d, iter %d Model: %s \n',imageNumid,iter,Modelfullname);
       
        % not handle occ, no local search, no pading  on training
        pad =0;
        Localsearch =0;
        [bbw_3d_d,bbs_3d_f] =detection3D_Range_selfocc(featureNeg,svmpara.numbox,SpaceNeg,Model,pointCountIntegral,zeros(size(pointCountIntegral)),zeros(Model.size_pos([1,2,3])),Localsearch,pad,featurepara.removedontcare,featurepara);    
        if ~isempty(bbw_3d_d)         
                % nms
                indexes  = nmsMe_3d([bbw_3d_d(:,1:3),bbw_3d_d(:,1:3)+bbw_3d_d(:,4:6),bbw_3d_d(:,7)],0.3); 
                bbs_3d_f = bbs_3d_f(indexes,:);
                bbw_3d_d = bbw_3d_d(indexes,:);
                %not overlap with ground Truth                
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
            t=toc;
            fprintf('New hard neigative: %d Totol: %d threshold: %d \n  Time used : %f\n',size(Nex,2),size(NegAll,2),hardnegthr,t);                        
        else
            fprintf('NegAll exceed max number !\n')
            break;
        end
    end

    thardneg=toc(hardnegstart);
    fprintf('Time to collect hard negs :%f\n',thardneg);

    %% get negatives base on new threshold
    if isfield(svmpara,'update')&&~svmpara.update
        hardnegthr =hardnegthr;
    else
        hardnegthr =update_hardnegthr(osall,scoreall,hardnegthr,svmpara.socrethre,allgtnum);
    end
    if ~isempty(NegAll) 
        NegAll =NegAll(:,NexscoreAll>hardnegthr);
    end
    numNewNex = size(NegAll,2);
    iterInfor =[iterInfor;iter,size(pSV,1),size(nSV,1),hardnegthr,sum(osall>0.25),sum(osall<0.25),numNewNex,imageNumid];

    %% converge or retrain
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
            careInd =ones(size_pos)>0;
        end
        Nex =[NegAll';nSV];
        Pos =[Pos_all';pSV];
        x = [Pos;Nex];
        y = [ones(size(Pos,1),1); -1*ones(size(Nex,1),1)];
        notnanInd=~isnan(sum(x,2));
        x=x(notnanInd,:);
        y=y(notnanInd,:);

        w1 = svmpara.w1;
        clear Pos Nex NegAll featureNeg bbs_3d_f ;
        if svmpara.linear 
            model = train(y,sparse(x), sprintf('-s 2 -B %f -c %f -w1 %f -q', 67500,svmpara.c, svmpara.w1));
            w = model.w(1:end-1);
            bias = model.w(end);
            nSV =[nSV;x(find(y(2:end).*(x(2:end,:)*w(:)+bias)<1)+1,:)];
            nSV =unique(nSV,'rows');
        else
            model = svmtrain(y, x, sprintf('-b 0 -s 0 -t 0 -c %.20f -w1 %.9f',svmpara.c, w1)); 
            w = model.sv_coef' * model.SVs;
            bias = -model.rho;
            pSV =model.SVs(find(model.sv_coef>0),:);
            nSV =[nSV;model.SVs(find(model.sv_coef<0),:)];
            nSV =unique(nSV,'rows');
        end
        template =zeros(size_pos);
        template(careInd>0) = w;
        svm =struct('w',template,'bias',bias);

        Model.hardnegthr =hardnegthr;
        Model.svm =svm;
        if numNewNex<5&&imageNumid==NumOfNeg&&hardnegthr<-0.93,
            break;
        end
        t=toc;
        fprintf('SVM %d retaining end !\n Time: %f \n', iter,t);
    end
     hardnegthr_low = -0.97;
end

%store the os and scores in last iteration for calibration later
Model.osall =osall;
Model.oscfall =oscfall;
Model.socreall =scoreall;
Model.iter = iterInfor;
%Model.modelscore=sum(Pos_all'*Model.svm.w(careInd));
if exist('savesvPath','var')&&~isempty(savesvPath),
   mkdir([savesvPath classNames{Model.classId} '/' Model.name '/'])
   save([savesvPath classNames{Model.classId} '/' Model.name '/' Modelfullname '_SV.mat'],'nSV','pSV');
end
mkdir(folder);
save(Modeltosavefullpath,'Model','-v6');
end
toc
end