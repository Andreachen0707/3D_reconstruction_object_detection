function [randNeg,neg_bb_w,neg_bb] =randgenNeg(size_pos,feature,bb_3d_f,postive_bb,Space,numTogen)
            neg_bb = genNex_bb(bb_3d_f,size(feature),numTogen);
            neg_bb_w =bbf2w(neg_bb,Space);
            % check not overlap with the ground true bbs
            if ~isempty(postive_bb)
                os=bb3dOverlapApprox(neg_bb_w,postive_bb);
                idxfar = ~sum(os>0.3,2);
                neg_bb_w = neg_bb_w(idxfar,:);
                neg_bb = neg_bb(idxfar,:);
            end
            randNeg = zeros(prod(size_pos),size(neg_bb,1));
            for i=1:size(neg_bb,1)
                nex = getTemplate_3d(feature,neg_bb(i,:));
                randNeg(:,i) = nex(:);
            end
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
 %rng(1);
 i1 = ceil(N*rand(1,k));
 out = in(:,i1);
end