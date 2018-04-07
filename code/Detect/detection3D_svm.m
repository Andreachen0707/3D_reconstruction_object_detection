function [bbw,bbs_3d_f] =detection3D_svm(feature,svm,num_dect,Space,useSigmoid,Cdepth)
        if ~exist('useSigmoid','var')
            useSigmoid =0;
        end
        bb_f_size =[size(svm.w,2)-1,size(svm.w,1)-1,size(svm.w,3)-1];
        if ~exist('Cdepth','var')
            % Search the whole space         

            match = evaluate_score_4d_svm(feature,svm.w,svm.bias);
            [scores,indexes] = sort(match(:),'descend');

            if numel(scores)>num_dect
                choosen = 1:num_dect;
                indexes = indexes(choosen);
                scores  = scores(choosen);
            else
                choosen=1:numel(scores);
            end
            [indexes1, indexes2, indexes3] = ind2sub(size(match),indexes);
            bbs_3d_f = [indexes2 indexes1 indexes3 repmat(bb_f_size,[size(choosen,2),1]) scores(:)];
            bbw =bbf2w(bbs_3d_f,Space);    
        else
            % Only search near by regoin (y)
            Cdepth_f =(Cdepth-Space.Ry(1))/Space.s;
            range_f_min = max(1,round(Cdepth_f-1/Space.s));
            range_f_max = min(size(feature,1),round(Cdepth_f+1/Space.s));
            featuresearch = feature(range_f_min:range_f_max,:,:,:);
            if useSigmoid
                match = evaluate_score_4d_svm(featuresearch,svm.w,svm.bias,svm.beta);
            else
                match = evaluate_score_4d_svm(featuresearch,svm.w,svm.bias);
            end
            
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
            bbs_3d_f = [indexes2 indexes1+range_f_min-1 indexes3 repmat(bb_f_size,[size(choosen,2),1]) scores(:)];
            bbw =bbf2w(bbs_3d_f,Space);    
        end
end