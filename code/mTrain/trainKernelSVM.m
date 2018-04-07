function SVMmodel= trainKernelSVM(Pos,Nex,svmpara)
        x = [Pos;Nex];
        y = [ones(size(Pos,1),1); -1*ones(size(Nex,1),1)];
        clear Pos Nex;
        notnanInd=~isnan(sum(x,2));
        x=x(notnanInd,:);y=y(notnanInd,:);
        SVMmodel = svmtrain(y, x, sprintf(['-b 0 -s 0 -t 2 -c %.50f -w1 %.9f'],svmpara.c,svmpara.w1));
end