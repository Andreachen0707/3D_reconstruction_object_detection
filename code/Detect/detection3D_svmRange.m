function [bbw,bbs_3d_f] =detection3D_svmRange(feature,svm,num_dect,Space,Cdepth,CHeight,Cx)
        bb_f_size =[size(svm.w,2)-1,size(svm.w,1)-1,size(svm.w,3)-1];
        
        % only search near by regoin 
        if exist('Cdepth','var')&&~isempty(Cdepth)
            Cdepth_f =(Cdepth-Space.Ry(1))/Space.s;
            R_d =max(size(svm.w,1),1/Space.s);
            range_f_dmin = max(1,round(Cdepth_f-R_d));
            range_f_dmax = min(size(feature,1),round(Cdepth_f+R_d));
        else
             range_f_dmin =1;
             range_f_dmax =size(feature,1);
        end
        
        if  exist('CHeight','var')&&~isempty(CHeight)
            % Cheight is height to ground 
            Cheight_f =(0+0.5-Space.Rz(1))/Space.s;
            range_f_hmin = 1;
            range_f_hmax = min(size(feature,3),min(round(1.5*size(svm.w,3)),round(Cheight_f)));
        else 
            range_f_hmin =1;
            range_f_hmax =size(feature,3);
        end
        if  exist('Cx','var')&&~isempty(Cx)
            Cx_f = (Cx -Space.Rx(1))/Space.s;
            R_x =max(size(svm.w,2),1/Space.s);
            range_f_xmin = max(1,round(Cx_f-R_x));          
            range_f_xmax = min(size(feature,2),round(Cx_f+R_x));
        else
            range_f_xmin =1;
            range_f_xmax =size(feature,2);
        end

                
        featuresearch = feature(range_f_dmin:range_f_dmax,range_f_xmin:range_f_xmax,range_f_hmin:range_f_hmax,:);
        
        % check whether the space is too small
        if ~isempty(find((size(featuresearch)-size(svm.w))<=0))
            bbw =[];
            bbs_3d_f =[];
            fprintf('the searching range in feature space is too small !');
            return;
        end

        match = evaluate_score_4d_svm(featuresearch,svm.w,svm.bias);
        [scores,indexes] = sort(match(:),'descend');

        if numel(scores)>num_dect
            choosen = 1:num_dect;
            indexes = indexes(choosen);
            scores  = scores(choosen);
        else
            choosen=1:numel(scores);
        end
        [indexes1, indexes2, indexes3] = ind2sub(size(match),indexes);
        % change back the whole world cordinate
        bbs_3d_f = [indexes2+range_f_xmin-1 indexes1+range_f_dmin-1 indexes3 repmat(bb_f_size,[size(choosen,2),1]) scores(:)];
        bbw =bbf2w(bbs_3d_f,Space);    
end