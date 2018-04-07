function bbf =bbw2f(bb_in,Space)
bbf =bb_in;
if ~isempty(bb_in)
bbf(:,1:3) =round((bb_in(:,1:3)-repmat([Space.Rx(1) Space.Ry(1) Space.Rz(1)],[size(bb_in,1),1]))/Space.s)+1;
bbf(:,4:6) =round(bb_in(:,4:6)/Space.s)-1;
bbf(bbf<=0) =1;
bbf =round(bbf);
end
end