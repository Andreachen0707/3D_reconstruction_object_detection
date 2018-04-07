function match = evaluate_score_4d_svm(featureM,templateM,bias,beta)

if size(templateM,4)~=size(featureM,4),
    size(featureM)
    size(templateM)
    error('svm weight has difference dimention from feature!');
end
match = fconv3d(featureM, {templateM}, 1, 1);
match = match{1};
match =match+bias;
if exist('beta','var')
   match=1./(1+exp(-beta(1)*(match-beta(2))));
   match =match-0.5;
end