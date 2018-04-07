% In this demo code will load one examplar model and test it on 3 image
% show the top 3 detections 
% For examples to train and test the full set of models, see full_demo.m
clear all;
close all;

tic;
initPath;
[svmpara,featurepara]=setParameters(); % set paramter
%featurepara.featype = 'custom';

%Allsvm ={'../preTrainedModels\chair\chair\chair_000000358\stair_r0.85_s1_d2_t-0.32_x0_h-0.95_sft000_asp111.mat'};
% rotate = [1];
% 
% for t = rotate
% t = num2str(t);
% name = ['../preTrainedModels\stair\stair\stair_001\stair_r' t '_s0.6_d2.5_t-0.32_x0_h-0.95_sft000_asp111.mat']
% Allsvm ={name};
% 
% 
% drawfigure =1;
% outpathTest ='';
% imageID = [48:120];
% replace =0;
% 
%  for i = imageID
%  TestingSingleImage_test(outpathTest,i,Allsvm,featurepara,drawfigure,replace,featurepara.localsearch);
%  end
% end
% %show the model
 offfile =  '../offData/stair/stair/stair.off'
% load(Allsvm{1});
% depth  = off2im_sys(offfile,1,Model.Rotation*pi,Model.Tilt,Model.Cx,Model.Cdepth,Model.Cheight,[1;1;1]*Model.Scale,1,0.1,0);
 for rotation = 0.6:0.1:1
 depth  = off2im_sys(offfile,1,pi,-0.5,0,2.5,-0.95,[1;1;1]*rotation,1,0.1,0);
figure,
imagesc(depth);
end
title('The depth map of this examplar')