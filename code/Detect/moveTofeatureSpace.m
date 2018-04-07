function [ptsSpace,XYZworldframeOut]=moveTofeatureSpace(XYZworldframe,Space)
        
        X =XYZworldframe(:,:,1);
        Y= XYZworldframe(:,:,2);
        Z =XYZworldframe(:,:,3);
        X(X<Space.Rx(1)) = NaN;
        X(X>Space.Rx(2)) = NaN;
        Y(Y<Space.Ry(1)) = NaN;
        Y(Y>Space.Ry(2)) = NaN;
        Z(Z<Space.Rz(1)) = NaN;
        Z(Z>Space.Rz(2)) = NaN;
        ISNAN=isnan(X)|isnan(Y)|isnan(Z);
        X(ISNAN)=NaN;
        Y(ISNAN)=NaN;
        Z(ISNAN)=NaN;
        XYZworldframeOut =cat(3,X,Y,Z);
        Xf =round((X+Space.Rx(3))/Space.s+1);
        Yf =round((Y+Space.Ry(3))/Space.s+1);
        Zf =round((Z+Space.Rz(3))/Space.s+1);     
        ptsSpace = cat(3,Xf,Yf,Zf);
end