function [w,pSV,nSV,bias]=MTrainLibSVM(Pos,Nex,svmpara)
        x = [Pos;Nex];
        y = [ones(size(Pos,1),1); -1*ones(size(Nex,1),1)];
        clear Pos Nex;
        notnanInd=~isnan(sum(x,2));
        x=x(notnanInd,:);y=y(notnanInd,:);
        if isempty(svmpara)
            model = svmtrain(y, x, sprintf(['-b 0 -s 0 -t 0 ']));
        else
            model = svmtrain(y, x, sprintf(['-b 0 -s 0 -t 0 -c %.50f -w1 %.9f'],svmpara.c,svmpara.w1));
        end
        w = model.sv_coef' * model.SVs;
        bias = -model.rho;
        pSV=model.SVs(find(model.sv_coef>0),:);
        nSV=model.SVs(find(model.sv_coef<0),:);
end