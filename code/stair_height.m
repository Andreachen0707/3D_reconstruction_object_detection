clear all;
load('depth_1.mat');
XYZ = zeros(480,640,3);
XYZ_2 = zeros(480,640,3);
corr_depth = zeros(4,2);
corr_height = zeros(4,2);
corr_depth_test = zeros(128,1);
corr_height_test = zeros(128,1);
idx = zeros(128,4);

fx = 616.029;
fy = 621.495;
cx = 330.764;
cy = 229.538;
tt = 1;
image = [117];

for i = 1:1 
imgNum =image(i);
tt = i;
load(['..\Mydatafeature_large/feature_' num2str(imgNum) '.mat']);
% namedepth = ['./pic_new/' num2str(imgNum) '_depth.png'];
% namergb = ['./pic_new/' num2str(imgNum) '_rgb.jpg'];
% imgDepth = imread(namedepth);
imgDepth_2 = depth;
% imgDepth = rgb2gray(imgDepth);
% imgDepth = double(imgDepth);
%imgDepth_2 = double(imgDepth_2);
% imgDepth = imgDepth/1000;
%imgDepth_2 = imgDepth_2/1000;
% x = find(imgDepth>5);
x_2 = find(imgDepth == 0);
y = find(imgDepth_2>5);
z = find(isnan(imgDepth_2));
% imgDepth(x) = 0.1;
imgDepth(x_2) = 0;
imgDepth_2(y) = 0;
imgDepth_2(z) = 0;


% figure,vis_point_cloud(XYZworldframeTest,rgbTest,10,5000);

v = 640;
w = 480;

for m = 1:v
    for n = 1:w
        XYZ(n,m,3) = imgDepth(n,m);
        XYZ(n,m,1) = (m - cx) * XYZ(n,m,3)/ fx;
        XYZ(n,m,2) = -(n - cy) * XYZ(n,m,3)/ fy;
        
        XYZ_2(n,m,3) = imgDepth_2(n,m);
        XYZ_2(n,m,1) = (m - cx) * XYZ_2(n,m,3)/ fx;
        XYZ_2(n,m,2) = -(n - cy) * XYZ_2(n,m,3)/ fy;
    end
end

% t = [1 0 0 -pi/10];
% q(1) = cos(t(4)/2);
% q(2:4) = t(1:3)*sin(t(4)/2);
% Rtilt =quat2rotm(q);

Rtilt_model = getRotationMatrix('tilt',-0.32);

point(:,:) = zeros(480,3);

for m = 1:v
    for n = 1:w
        point(n,:) = [XYZ(n,m,1) XYZ(n,m,3) XYZ(n,m,2)];
        point(n,:) = (Rtilt*point(n,:)')';
        XYZ(n,m,1) = point(n,1);
        XYZ(n,m,2) = point(n,3);
        XYZ(n,m,3) = point(n,2);
        
         point(n,:) = [XYZ_2(n,m,1) XYZ_2(n,m,3) XYZ_2(n,m,2)];
         point(n,:) = (Rtilt_model*point(n,:)')';
         XYZ_2(n,m,1) = point(n,1);
         XYZ_2(n,m,2) = point(n,3);
         XYZ_2(n,m,3) = point(n,2);
    end
end

depth_model = XYZ_2(201:480,196:445,3);
height_model = XYZ_2(201:480,196:445,2);

for p = 1:2
    depth_model_test = depth_model+0.007*(p-1)*140;
    height_model_test = height_model+0.003*(p-1);
    for t = 1:4
        weight_1 = 201-(p-1)*140;
        weight_2 = 480-(p-1)*140;
        height_1 = 1+(t-1)*130;
        height_2 = 250+(t-1)*130;
        depth_test =  XYZ(weight_1:weight_2,height_1:height_2,3);
        height_test =  XYZ(weight_1:weight_2,height_1:height_2,2);
        corr_depth(t,p) =corr2(depth_test,depth_model_test);
        corr_height(t,p) =corr2(height_test,height_model_test);
    end
    
end


corr_depth_test(i,1) = max(max(corr_depth));
corr_height_test(i,1) = max(max(corr_height));

 row = 2;
 line = 2;
[row,line] = find(corr_depth==corr_depth_test(i,1));
idx(i,1) = row(1,1);
idx(i,2) = line(1,1);
[row,line]=find(corr_height==corr_height_test(i,1));
idx(i,3) = row(end,1);
idx(i,4) = line(end,1);


%% show
figure;
imshow(imgRgb);hold on;
point1 = [201-(idx(tt,2)-1)*140 1+(idx(tt,1)-1)*130];
point2 = [480-(idx(tt,2)-1)*140 1+(idx(tt,1)-1)*130];
point3 = [480-(idx(tt,2)-1)*140 250+(idx(tt,1)-1)*130];
point4 = [201-(idx(tt,2)-1)*140 250+(idx(tt,1)-1)*130];

plot([point1(1,2),point2(1,2)],[point1(1,1),point2(1,1)],'-', 'LineWidth',2,'Color','r');
plot([point1(1,2),point4(1,2)],[point1(1,1),point4(1,1)],'-', 'LineWidth',2,'Color','r');
plot([point2(1,2),point3(1,2)],[point2(1,1),point3(1,1)],'-', 'LineWidth',2,'Color','r');
plot([point4(1,2),point3(1,2)],[point4(1,1),point3(1,1)],'-', 'LineWidth',2,'Color','r');
text(point1(1,2),point1(1,1),sprintf('%.2f',corr_depth_test(tt,1)),'background','w','FontSize',20);
%text(point1(1,2),point1(1,1),sprintf('%.2f',corr_depth(2,2)),'background','w','FontSize',20);

% bar(depth);
%plot(x,depth);
end
save feature_x corr_depth_test;