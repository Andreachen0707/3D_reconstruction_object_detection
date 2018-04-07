function [pts,xyzShifted,cellIdx,gridIdx,remainder,bincellDim,histSize,Space,removedpts]=init_gridInd(XYZworldframe,s,Qnum,shift)
xmin =nanmin(XYZworldframe(:,1));
ymin =nanmin(XYZworldframe(:,2));
XYZworldframe(XYZworldframe(:,3)<-3,:) = NaN;
zmin =nanmin(XYZworldframe(:,3));
xmin =xmin;
ymin =ymin;
zmin =zmin-0.5*s;



% compute point position in grid
stepSize = s/Qnum;
xyzShifted = bsxfun(@minus,XYZworldframe,[xmin,ymin,zmin]+shift);
removedpts = sum(xyzShifted<0,2)|isnan(XYZworldframe(:,1));
pts =XYZworldframe;
pts = pts(~removedpts,:);
xyzShifted = xyzShifted(~removedpts,:);

gridIdx = floor(bsxfun(@rdivide,xyzShifted,[stepSize,stepSize,stepSize])) + 1; 
remainder = (bsxfun(@rem,xyzShifted,[stepSize,stepSize,stepSize])) / stepSize;

% compute the cell ind
cellIdx = floor(bsxfun(@rdivide,gridIdx-1,[Qnum,Qnum,Qnum]))+1;
bincellDim = max(cellIdx+1,[],1);
histSize = bincellDim*Qnum; % histSize =xyz
bincellDim =[bincellDim(2),bincellDim(1),bincellDim(3)];% bincellDim yxz

% round the Space to be int grid and Cell 
xmax =xmin+bincellDim(2)*s;
ymax =ymin+bincellDim(1)*s;
zmax =zmin+bincellDim(3)*s;
Space=struct('Rx',[xmin, xmax],'Ry',[ymin,ymax],'Rz',[zmin,zmax],'s',s);

end
