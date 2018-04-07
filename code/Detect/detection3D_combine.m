function  [bbw,bbs_3d_f,occCount] =detection3D_combine(featureTest,Allsvm,num_dect,SpaceTest, pointCountIntegral,occSource,cellInsight,partailInsight,pointCount,missDepthcell,Modelocc,localsearch,combine)
for modeltypeId = 1:length(Allsvm)
    load(Allsvm{modeltypeId},'Model');
    if ~isfield(Model,'featurepara')
        [svmpara,featurepara]=setParameters();
        featurepara.featype = Model.featype;  
    else
        featurepara = Model.featurepara;
    end
    featureTestType = getfeatureTypeDim(featureTest,featurepara,cellInsight,partailInsight,missDepthcell,pointCount,occSource);
    if featurepara.occ,
        [~,~,~,matchType{modeltypeId},searchRage,occCountAll] =detection3D_Range_selfocc(featureTestType,1000,SpaceTest,Model,pointCountIntegral,occSource,...
                        Modelocc,localsearch,featurepara.pad,featurepara.removedontcare,featurepara); 
    else
        [~,~,~,matchType{modeltypeId},searchRage,occCountAll] =detection3D_Range_selfocc(featureTestType,1000,SpaceTest,Model,pointCountIntegral,zeros(size(pointCountIntegral)),...
                        zeros(Model.size_pos([1,2,3])),localsearch,featurepara.pad,featurepara.removedontcare,featurepara); 
    end
end
% match final 
if ~isempty(matchType{modeltypeId}),
   match =    combine.w(1)*matchType{1}+ combine.w(2)*matchType{2};
% return bbw,bbs_3d_f
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
occCount = occCountAll(indexes);
[indexes1, indexes2, indexes3] = ind2sub(size(match),indexes);
% change back to the whole world cordinate
bb_f_size =[size(Model.svm.w,2)-1,size(Model.svm.w,1)-1,size(Model.svm.w,3)-1];
bbs_3d_f = [indexes2+searchRage.range_f_xmin-1 indexes1+searchRage.range_f_dmin-1 indexes3+searchRage.range_f_hmin-1 repmat(bb_f_size,[size(choosen,2),1]) scores(:)];
bbw =bbf2w(bbs_3d_f,SpaceTest); 
else
    bbs_3d_f=[];
    bbw=[];
    occCount=[];
end
   
end
     