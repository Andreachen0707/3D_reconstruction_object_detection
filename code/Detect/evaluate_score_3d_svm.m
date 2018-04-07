function match = evaluate_score_3d_svm(featureM,templateM,bias,beta)

if size(templateM,4)~=size(featureM,4),
    size(featureM)
    size(templateM)
    error('svm weight has difference dimention from feature!');

end

match =0;
for f =1:size(templateM,4)
    feature =featureM(:,:,:,f);
    template =templateM(:,:,:,f);
    level=size(feature,3)-size(template,3)+1;
    for k=1:level
        match_level = fconvblas(feature(:,:,k:(k+size(template,3)-1)), {template}, 1, 1);
        match_f(:,:,k) = match_level{1};
        %pause;       
    end
    match =match+match_f;
end

match =match+bias;

if exist('beta','var')
   match=1./(1+exp(-beta(1)*(match-beta(2))));
   match =match-0.5;
end