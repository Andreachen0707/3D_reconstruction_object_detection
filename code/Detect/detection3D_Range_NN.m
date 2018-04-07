function [bbw,bbs_3d_f] =detection3D_Range_NN(feature,num_dect,Space,Model,pointCountIntegral,occSource,Modelocc,Localsearch,pad,removedontcare,featurepara)
        bbw =[];
        bbs_3d_f =[];
        Cdepth = Model.Cdepth;
        CHeight = Model.Cheight;
        Cx = Model.Cx;

        bb_f_size =[Model.size_pos(2)-1,Model.size_pos(1)-1,Model.size_pos(3)-1];
        
        % only search near by regoin 
        if exist('Cdepth','var')&&~isempty(Cdepth)&&Localsearch
            Cdepth_f =(Cdepth-Space.Ry(1))/Space.s;
            R_d =max(Model.size_pos(1),1/Space.s);
            range_f_dmin = max(1,round(Cdepth_f-R_d));
            range_f_dmax = min(size(feature,1),round(Cdepth_f+R_d));
        else
             range_f_dmin =1;             
             range_f_dmax =size(feature,1);
        end
        
        if  exist('CHeight','var')&&~isempty(CHeight)&&Localsearch
            % Cheight is height to ground 
            Cheight_f =(0-Space.Rz(1))/Space.s+0.5*Model.size_pos(3);
            range_f_hmin = 1;
            range_f_hmax = min(size(feature,3),max(round(1.5*Model.size_pos(3)),round(Cheight_f))); 
            
        else 
            range_f_hmin =1;
            range_f_hmax =size(feature,3);
        end
        if  exist('Cx','var')&&~isempty(Cx)&&Localsearch
            Cx_f = (Cx -Space.Rx(1))/Space.s;
            R_x =max(Model.size_pos(2),1/Space.s);
            range_f_xmin = max(1,round(Cx_f-R_x));          
            range_f_xmax = min(size(feature,2),round(Cx_f+R_x));
            if isfield(featurepara,'Localview')&&featurepara.Localview
               cameraLine = (0 -Space.Rx(1))/Space.s;
               if Cx<-0.2,
                   range_f_xmax = min(round(cameraLine+0.8*Model.size_pos(2)),range_f_xmax);
               elseif Cx>0.2
                   range_f_xmin = max(round(cameraLine-0.8*Model.size_pos(2)),range_f_xmin);
               end
            end
        else
            range_f_xmin =1;
            range_f_xmax =size(feature,2);
        end

        featuresearch = feature(range_f_dmin:range_f_dmax,range_f_xmin:range_f_xmax,range_f_hmin:range_f_hmax,:);
        pointCountIntegral = pointCountIntegral(range_f_dmin:range_f_dmax,range_f_xmin:range_f_xmax,range_f_hmin:range_f_hmax,:);
        occSource = occSource(range_f_dmin:range_f_dmax,range_f_xmin:range_f_xmax,range_f_hmin:range_f_hmax,:);

        if pad
           padnumofcell =round(Model.size_pos(3));
           padfea = zeros(size(featuresearch,1),size(featuresearch,2),padnumofcell,size(featuresearch,4));
           padfea(:,:,:,end) =1;
           %padpointCountIntegral = repmat(pointCountIntegral(:,:,end),[1,1,padnumofcell]);
           padpointCountIntegral = zeros(size(pointCountIntegral,1),size(pointCountIntegral,2),padnumofcell);
           padMaskfeature = zeros(size(occSource,1),size(occSource,2),padnumofcell);
           occSource =cat(3,padMaskfeature,occSource);
           pointCountIntegral =cat(3,padpointCountIntegral,pointCountIntegral);
           featuresearch =cat(3,padfea,featuresearch);
           range_f_hmin = range_f_hmin -padnumofcell;
        end
        size_featuresearch =size(featuresearch);
        size_w =Model.size_pos;
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
        
        tic;
        if min(size(careMask))==0,
            return;
        end
        match = evaluate_score_4d_ssd(featuresearch,Model,pointCountIntegral,occSource,Modelocc,careMask,[Space.s,SpaceRangeDepth],featurepara);
        toc;
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

function match = evaluate_score_4d_ssd(featureM,Model, pointCountIntegral,occSource,Modelocc,careMask,SpaceRange,featurepara)
    careMask= repmat(careMask,[1,1,1,size(featureM,4)]);
    sizeTemplate =size(careMask);
    sizefeature = size(featureM);

    EmptyBox = EmptyBoxFlag(pointCountIntegral,sizeTemplate,200);
    [dim1,dim2,dim3]=ndgrid(1:sizefeature(1)-sizeTemplate(1)+1,1:sizefeature(2)-sizeTemplate(2)+1,1:sizefeature(3)-sizeTemplate(3)+1);
    match = -100*ones(size(dim2));
    dim1 = dim1(~EmptyBox);
    dim2 = dim2(~EmptyBox);
    dim3 = dim3(~EmptyBox);
    for k = 1:length(dim1(:)),
         boxfea=getTemplate_3d(featureM,[dim2(k),dim1(k),dim3(k),sizeTemplate(2)-1,sizeTemplate(1)-1,sizeTemplate(3)-1]); 
         match(dim1(k),dim2(k),dim3(k))  = (norm(boxfea(:) - Model.pos(:)));
    end
    match = sum(sizeTemplate)*100./(match+1);
end