function Space= initSpace(XYZworldframe,s,enlarge,bb3dvmat)
    if ~exist('enlarge','var')
        enlarge =0;
    end
    xmin =nanmin(XYZworldframe(:,1))-enlarge;
    xmax=nanmax(XYZworldframe(:,1))+enlarge;
    ymin =nanmin(XYZworldframe(:,2))-enlarge;
    ymax=nanmax(XYZworldframe(:,2))+enlarge;
    zmin =nanmin(XYZworldframe(:,3))-enlarge;
    zmax=nanmax(XYZworldframe(:,3))+enlarge;

    if exist('bb3dvmat','var')
       %zmin =min([nanmax(nanmax(XYZworldframe(:,3)))]) -bb3dvmat(6) -enlarge-0.5*s;
       center =[bb3dvmat(1)+0.5*(bb3dvmat(4)),bb3dvmat(2)+0.5*(bb3dvmat(5)),0.5*(zmax+zmin)];
       xmin = center(1) -0.5*bb3dvmat(4) -enlarge;
       ymin = center(2) -0.5*bb3dvmat(5) -enlarge;
       zmin = center(3) -0.5*bb3dvmat(6) -enlarge; 
       xmax = center(1)+0.5*bb3dvmat(4) + enlarge;
       ymax = center(2)+0.5*bb3dvmat(5) + enlarge;
       zmax = center(3)+0.5*bb3dvmat(6) + enlarge; 
    end
    xmin =xmin-0.5*s;
    ymin =ymin-0.5*s;
    zmin =zmin-0.5*s;
    
    Space=struct('Rx',[xmin, xmax],'Ry',[ymin,ymax],'Rz',[zmin,zmax],'s',s);
end