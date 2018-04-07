function [bbw,bbs_3d_f] =detection3D_Range_occ(feature,num_dect,Space,Model,pointCountIntegral,Maskfeature,MaskTemp,Localsearch,pad,removedontcare,featurepara)
        svm =Model.svm;
        Cdepth = Model.Cdepth;
        CHeight = Model.Cheight;
        Cx = Model.Cx;

        bb_f_size =[size(svm.w,2)-1,size(svm.w,1)-1,size(svm.w,3)-1];
        
        % only search near by regoin 
        if exist('Cdepth','var')&&~isempty(Cdepth)&&Localsearch
            Cdepth_f =(Cdepth-Space.Ry(1))/Space.s;
            R_d =max(size(svm.w,1),1/Space.s);
            range_f_dmin = max(1,round(Cdepth_f-R_d));
            range_f_dmax = min(size(feature,1),round(Cdepth_f+R_d));
        else
             range_f_dmin =1;             
             range_f_dmax =size(feature,1);
        end
        
        if  exist('CHeight','var')&&~isempty(CHeight)&&Localsearch
            % Cheight is height to ground 
            Cheight_f =(0+0.5-Space.Rz(1))/Space.s;
            range_f_hmin = 1;
            range_f_hmax = min(size(feature,3),max(round(1.5*size(svm.w,3)),round(Cheight_f))); 
        else 
            range_f_hmin =1;
            range_f_hmax =size(feature,3);
        end
        if  exist('Cx','var')&&~isempty(Cx)&&Localsearch
            Cx_f = (Cx -Space.Rx(1))/Space.s;
            R_x =max(size(svm.w,2),1/Space.s);
            range_f_xmin = max(1,round(Cx_f-R_x));          
            range_f_xmax = min(size(feature,2),round(Cx_f+R_x));
        else
            range_f_xmin =1;
            range_f_xmax =size(feature,2);
        end

        featuresearch = feature(range_f_dmin:range_f_dmax,range_f_xmin:range_f_xmax,range_f_hmin:range_f_hmax,:);
        pointCountIntegral = pointCountIntegral(range_f_dmin:range_f_dmax,range_f_xmin:range_f_xmax,range_f_hmin:range_f_hmax,:);
        Maskfeature = Maskfeature(range_f_dmin:range_f_dmax,range_f_xmin:range_f_xmax,range_f_hmin:range_f_hmax,:);

        if pad
           padnumofcell =round(size(svm.w,3));
           padfea = zeros(size(featuresearch,1),size(featuresearch,2),padnumofcell,size(featuresearch,4));
           padfea(:,:,:,end) =1;
           %padpointCountIntegral = repmat(pointCountIntegral(:,:,end),[1,1,padnumofcell]);
           padpointCountIntegral = zeros(size(pointCountIntegral,1),size(pointCountIntegral,2),padnumofcell);
           padMaskfeature = zeros(size(Maskfeature,1),size(Maskfeature,2),padnumofcell);
           Maskfeature =cat(3,padMaskfeature,Maskfeature);
           pointCountIntegral =cat(3,padpointCountIntegral,pointCountIntegral);
           featuresearch =cat(3,padfea,featuresearch);
           range_f_hmin = range_f_hmin -padnumofcell;
        end
        size_featuresearch =size(featuresearch);
        size_w =size(svm.w);
        % check whether the space is too small
        if sum(size_featuresearch([1,2]) -size_w([1,2])<0)||sum(size_featuresearch(3) -size_w(3)<=0)
            sum(size_featuresearch([1,2]) -size_w([1,2])<0)
            sum(size_featuresearch(3) -size_w(3)<=0)
            bbw =[];
            bbs_3d_f =[];
            fprintf('the searching range in feature space is too small !\n');
            return;
        end  
        
        if exist('removedontcare','var')&&removedontcare
           careMask = Model.careMask; 
        else
           careMask = ones(size_w([1,2,3])); 
        end
        SpaceRangeDepth = (range_f_dmin-1)*Space.s+Space.Ry(1);
        
        match = evaluate_score_4d_flag(featuresearch,svm.w,svm.bias, pointCountIntegral,Maskfeature,MaskTemp,careMask,[Space.s,SpaceRangeDepth],featurepara);
        [scores,indexes] = sort(match(:),'descend');
        % throw away detection with low conf
        choosen = scores>-10;
        indexes = indexes(choosen);
        scores  = scores(choosen);
        % get top boxes
        if numel(scores)>num_dect 
            choosen = 1:num_dect;
            indexes = indexes(choosen);
            scores  = scores(choosen);
        else
            choosen=1:numel(scores);
        end
        [indexes1, indexes2, indexes3] = ind2sub(size(match),indexes);
        % change back to the whole world cordinate
        bbs_3d_f = [indexes2+range_f_xmin-1 indexes1+range_f_dmin-1 indexes3+range_f_hmin-1 repmat(bb_f_size,[size(choosen,2),1]) scores(:)];
        bbw =bbf2w(bbs_3d_f,Space);    
end

function match = evaluate_score_4d_flag(featureM,templateM,bias, pointCountIntegral,Maskfeature,MaskTemp,careMask,SpaceRange,featurepara)
sizeTemplate =size(templateM);
sizefeature =size(featureM);
if sizeTemplate(4)~=sizefeature(4),
    size(featureM)
    size(templateM)
    error('svm weight has difference dimention from feature!');
end
EmptyBox = EmptyBoxFlag(pointCountIntegral,sizeTemplate,200);
% sum(EmptyBox(:))/size(EmptyBox(:),1)
% swap dim for C++ to easy read 
featureM = permute(featureM,[4 1 2 3]);
templateM = permute(templateM,[4 1 2 3]);
% // C = fconv(A, B, MA, MB, E, dontCare, SpaceY);
% // A is feature matrix
% // B is filter
% // MA is mask for feature:  0 not occluded, >0 for occluded cell, value depth (y value) of occluder.
% // MB is mask for filter:  1 for occluded cell, 0 otherwise.
% // E is isEmpty:  1 for empty position, 0 otherwise.
% // Space = [s, minY]
match = fconv3d_occ(featureM, templateM, Maskfeature, MaskTemp>0,EmptyBox,~careMask,SpaceRange,featurepara.addempty,featurepara.scalemissing);
match =match+bias;
end