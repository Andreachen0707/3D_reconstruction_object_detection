function [startIndLin,endIndLin]=get_start_end_ind(Qnum,sampleNum)

randomLoc=load('randomLoc.mat');
startInd =round((randomLoc.startLoc(:,1:sampleNum))*(Qnum-1)+1);
endInd =round((randomLoc.endLoc(:,1:sampleNum))*(Qnum-1)+1);
startIndLin =sub2ind([Qnum,Qnum,Qnum],...
          startInd(1,:),startInd(2,:),startInd(3,:));
endIndLin =sub2ind([Qnum,Qnum,Qnum],...
          endInd(1,:),endInd(2,:),endInd(3,:));
end