function hardnegthr_updated =update_hardnegthr(osall,socreall,hardnegthr,socrethre,allgtnum)
            posthr= 0.3;
            negthr= 0.1;
            posscore =socreall(osall>posthr);
            negscore =socreall(osall<negthr);
            if ~isempty(posscore)&&size(posscore,1)>0.1*allgtnum,
                [s,ind]=sort(posscore,'descend');
                hardnegthr_updated = max(-0.97,posscore(min(size(posscore,1),ind(ceil(socrethre*length(posscore)))))); 
                % threshold can only decrease
                hardnegthr_updated = min(hardnegthr_updated,0.3*hardnegthr-0.7*0.97);
            else
                hardnegthr_updated =0.3*hardnegthr-0.7*0.97;
            end
end