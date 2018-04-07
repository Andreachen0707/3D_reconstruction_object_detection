function feacode = sparse_coding(fearaw,center,stdf,KNN)
    % calculate distance to center 
    fea1 =pdist2(fearaw,center);  
    
    fea1 = exp(-(fea1.*fea1)/(stdf*stdf));

    %coding 
    if ~isempty(KNN)&&KNN<size(center,1)
        [~,sortIndex] = sort(fea1(:),'descend');
        maxind=sortIndex(1:KNN);
        fea=zeros(size(fea1));
        fea(maxind) = fea1(maxind);
    else
        fea = fea1;
    end
    %normalize 
    if max(fea(:))>0.00000001,
       fea=fea/max(fea(:));
    end
    fea(isnan(fea))=0;
    feacode =fea; 
end

