function [svm,neg_bb_w] = trainSVM(Pos,size_pos,feature,highscore_bb,svm,bb_3d_f,postive_bb,Space)
if isempty(highscore_bb)
    neg_bb = genNex_bb(bb_3d_f,size(feature),100);
else 
    neg_bb = highscore_bb; 
end
% check not overlap with the ground true bbs
if ~isempty(postive_bb)
   neg_bb_w =bbf2w(neg_bb,Space);
   os=bb3dOverlapApprox(neg_bb_w,postive_bb);
   idxfar = ~sum(os>0.3,2);
   neg_bb_w = neg_bb_w(idxfar,:);
   neg_bb = neg_bb(idxfar,:);
else
   neg_bb_w =[];
end
Nex = zeros(prod(size_pos),size(neg_bb,1));
for i=1:size(neg_bb,1)
   nex = getTemplate_3d(feature,neg_bb(i,:));
   Nex(:,i) = nex(:);
end
if ~isempty(svm)
    Nex =[Nex';svm.nSV(1:min(100,size(svm.nSV,1)),:)];
    Pos =[Pos';svm.pSV(1:min(100,size(svm.pSV,1)),:)];
    w0=svm.w;
    opt =struct('w0',w0(:));
else
    Pos=Pos';
    Nex=Nex'; 
    opt =[];
end

x = [Pos;Nex];
y = [ones(size(Pos,1),1); -1*ones(size(Nex,1),1)];
notnanInd=~isnan(sum(x,2));
x=x(notnanInd,:);
y=y(notnanInd,:);
%[w,b,sv]=primal_svm(1,x,y,1,opt);

%svminputname =[num2str((rand()*1000)) '.mat'];
%save(svminputname,'y','x');
model = svmtrain(y, x, sprintf(['-b 1 -s 0 -t 0 -c %f -w1 %.9f'],0.01, 50));
%delete(svminputname);

w = model.sv_coef' * model.SVs;
bias = -model.rho;
% support vectors: full(model.SVs)
% weight: model.sv_coef (positive and negative)
template=reshape(w,size_pos);
pSV=model.SVs(find(model.sv_coef>0),:);
nSV=model.SVs(find(model.sv_coef<0),:);
%{
pscore=w(:)'*pSV';
nscore=w(:)'*nSV';
pscore=sort(pscore);
nscore=sort(nscore);
pind=max(1,round(0.5*length(pscore)));
nind=max(1,round(0.95*length(nscore)));
c=(0.5*pscore(pind)+0.5*nscore(nind));
if isempty(svm)
    thr=c;
else
    thr=0.7*svm.thr+0.3*c;
end
%}
svm =struct('w',template,'pSV',pSV,'nSV',nSV,'bias',bias);
svm.model =model;
end
function neg_bb = genNex_bb(bb_3d_f,Dim3,n) 
    x= randvalues(round(1:Dim3(2)-bb_3d_f(4)),n);
    y= randvalues(round(1:Dim3(1)-bb_3d_f(5)),n);
    z= randvalues(round(1:Dim3(3)-bb_3d_f(6)),n);
    neg_bb=[x',y',z',repmat(bb_3d_f(4:6),[n,1])];    
end

function out = randvalues(in,k)
% Randomly selects 'k' values from vector 'in' with replacement repeat.
out = [];
N = size(in,2);
if k == 0
  return;
end
 rng(1);
 i1 = ceil(N*rand(1,k));
 out = in(:,i1);
end