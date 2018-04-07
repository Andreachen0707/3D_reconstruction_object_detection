% input : 
% example store in : data.mat
% K: camera intrinsic matrix 
% imgDepth: depth map 
% Rtilt: Rotation matrix to align the point could on gravity direction, for
% a room layout the Rtilt can be estimate from [Rtilt,R] = rectify(XYZ);

% Output : 
% featureTest : MxNxKx300 matrix
% SpaceTest : point cloud range
% bincellDim : [M N K]
% cellInsight : indicate whether each cell completely inside view of field
% partailInsight : indicate whether each cell partially inside view of field
% occcell : whether each cell is occluded or not 
% missDepthcell : whether this cell is missing depth
% pointCount : number of points in each cell
% pointCountIntegral : 3D integral image of pointCount
% occSource = for each occluded cell what's the depth of occluder 

%% start
clear all
close all

type = 2;

if type ==1
% fx = 301.433;
% fy = 301.433;
% cx = 159.5;
% cy = 119.281;

fx = 5.1885790117450188e+02;
fy = 5.1946961112127485e+02;
cx = 3.2558244941119034e+02;
cy = 2.5373616633400465e+02;
end

if type ==2
    fx = 616.029;
    fy = 621.495;
    cx = 330.764;
    cy = 229.538;
end

if type ==3
fx = 301.433/2;
fy = 301.433/2;
cx = 159.5/2;
cy = 119.281/2;
end

K = [fx 0 cx; 0 fy cy; 0 0 1];

 imgNum = 3;
 namedepth = ['./pic_new/' num2str(imgNum) '_depth.png'];
 namergb = ['./pic_new/' num2str(imgNum) '_rgb.jpg'];
 imgDepth = imread(namedepth);
% %imgDepth = imread('52.png');
imgDepth = double(imgDepth);
imgDepth = imgDepth/1000;
x = find(imgDepth>3.5);
imgDepth(x) = 0;
% % y = find(imgDepth==0);
% %imgDepth(y) = 0.5;
imgRgb = imread(namergb);

featurepara.featype ='pointtsdfnormalshape'; 
%featurepara.featype ='tsdf';
featurepara.s =0.1;
featurepara.normalize=0;
featurepara.highdim=0;
crop = [1,1];

[rgbTest,XYZworldframeTest,XYZ] = read_3d_pts_general(imgDepth,K,size(imgDepth),imgRgb,crop);
%[rgbTest,XYZworldframe,XYZ] = read_3d_pts_general(imgDepth,K,size(imgDepth),[],crop);
% [Rtilt,R] = rectify(XYZ);
t = [1 0 0 -pi/10];
q(1) = cos(t(4)/2);
q(2:4) = t(1:3)*sin(t(4)/2);
Rtilt =quat2rotm(q);
%Rtilt = [1 0 0; 0 1 0; 0 0 1];

%% resize multiple frames
% XYZworldframeTest = load('XYZ.txt');
% y_least = find(XYZworldframeTest(:,3) <-1.6);
% %z_least = find(XYZworldframeTest(:,2) >5 );
% XYZworldframeTest(:,y_least) = 0;
% %XYZworldframeTest(:,z_least,:) = 0;
% rgbTest = load('rgb.txt');
% rgbTest = uint8(rgbTest);
% rgbTest = im2double(rgbTest);


%% typical
XYZworldframeTest = (Rtilt*XYZworldframeTest')';

figure,vis_point_cloud(XYZworldframeTest,rgbTest,10,5000);

[featureTest,SpaceTest,bincellDim,cellInsight,partailInsight,occcell,missDepthcell,pointCountIntegral,pointCount,occSource]...
            = mFeature_combined(imgDepth,XYZworldframeTest,featurepara.s,Rtilt,crop,[],featurepara,[0,0,0],0);
imgDepth = single(imgDepth);

clear bincellDim crop cx cy fx fy K XYZ featurepara occcell R x y XYZworldframe type t q namedepth namergb
save(['feature_' num2str(imgNum)]);
%save feature_69

