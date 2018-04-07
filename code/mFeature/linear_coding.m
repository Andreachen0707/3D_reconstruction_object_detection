function feature_coded = linear_coding(feature,featurepara)
if featurepara.floatcode,
    feature_coded = linear_coding_f(feature,featurepara);
else 
    feature_coded = linear_coding_q(feature,featurepara);
end
        
end
function feature_coded = linear_coding_f(feature,featurepara)
        if ~isfield(featurepara,'numofstep')
            featurepara.numofstep =10;
        end
        range = [featurepara.mindistance,featurepara.maxdistance]+0.01;
        step = (featurepara.maxdistance- featurepara.mindistance)/featurepara.numofstep;
        numStep = round((range(2)-range(1))/step)+1;
        feature(feature<range(1)) = range(1);
        feature(feature>range(2)) = range(2);

        feature_q = floor(feature/step)+1;
        feature_r = feature-(feature_q-1).*step;
        codebook = tril(ones(numStep,numStep),-1)';
        codebook_r = eye(numStep);

        feature_coded = zeros(size(feature_q,1),size(feature_q,2),size(feature_q,3),numStep*size(feature_q,4));
        [ind1,ind2,ind3]=ndgrid([1:size(feature_q,1)],[1:size(feature_q,2)],[1:size(feature_q,3)]);
        for k = 1:length(ind1(:))
            D = codebook(:,feature_q(ind1(k),ind2(k),ind3(k),:));
            R = zeros(size(D));
            R(codebook_r(:,feature_q(ind1(k),ind2(k),ind3(k),:))>0) = feature_r(ind1(k),ind2(k),ind3(k),:);
            fea_code = D(:)+R(:);
            feature_coded(ind1(k),ind2(k),ind3(k),:) = fea_code(:);
        end
        
        if ~(isfield(featurepara,'negativefea')&&featurepara.negativefea)
           feature_coded(feature_coded==0) =-1;
        end
end

function feature_coded = linear_coding_q(feature,featurepara)
        if ~isfield(featurepara,'numofstep')
            featurepara.numofstep =10;
        end
        range = [featurepara.mindistance,featurepara.maxdistance]+0.01;
        step = (featurepara.maxdistance- featurepara.mindistance)/featurepara.numofstep;
        numStep = round((range(2)-range(1))/step)+1;
        feature(feature<range(1)) = range(1);
        feature(feature>range(2)) = range(2);

        feature_q = floor(feature/step)+1;
        
        codebook = tril(ones(numStep,numStep),-1)';
        if ~isfield(featurepara,'negativefea')&&featurepara.negativefea
            feature_coded = -1*ones(size(feature_q,1),size(feature_q,2),size(feature_q,3),numStep*size(feature_q,4));
        else
            feature_coded = zeros(size(feature_q,1),size(feature_q,2),size(feature_q,3),numStep*size(feature_q,4));
        end
    %    feature_coded = zeros(size(feature_q,1),size(feature_q,2),size(feature_q,3),numStep*size(feature_q,4));
        [ind1,ind2,ind3]=ndgrid([1:size(feature_q,1)],[1:size(feature_q,2)],[1:size(feature_q,3)]);
        for k = 1:length(ind1(:))
            fea_code = codebook(:,feature_q(ind1(k),ind2(k),ind3(k),:));
            feature_coded(ind1(k),ind2(k),ind3(k),:) = fea_code(:);
        end

%}        
end
