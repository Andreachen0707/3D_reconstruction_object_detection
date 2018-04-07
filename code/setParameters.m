function [svmpara,featurepara]=setParameters()
% paramters for SVM training 
svmpara.socrethre =0.3;
svmpara.numbox =1000;
svmpara.NumOfNexLimit =500;% max number of negatives in each iteration of hard negative mining 
svmpara.otherNeg =0;
svmpara.c =1;
svmpara.maxiter =4; % max number of iteration for hard negative mining 
svmpara.w1 =1000;
svmpara.exculdepos = 1; 
svmpara.update =1; % update hardnegative ming threshold for each iteration
featurepara.addbias =1;
svmpara.linear =0;

% paramters for feature computation 
featurepara.outadd =1;
featurepara.missingadd =1;
featurepara.scalemissing =2;
featurepara.addempty =1;
featurepara.selfocc =0;
%  For cells with missing data: 
%  1. use original feature: f 2. use empty cell feature: E  3.use zeros as feature (don't care):0 
featurepara.missing='0';
%  For cells outside field of view:  missing data: p0tE p0t0 pft0 pEtE pftE 
featurepara.outtype='pft0'; 
% %  path to load precomputed feature (NYU)
% featurepara.featurepath = '../NYUdatafeature_release/'; 
% featurepara.pointpath   = '../NYUdatafeature_release/';

featurepara.featurepath = '../Mydatafeature_large/'; 
featurepara.pointpath   = '../Mydatafeature_large/';

%  path to save postive images
featurepara.modelImageFolder  = './Output/WebResult/modelImageFolder/';
%  path to load precomputed segmetation (NYU)
featurepara.segpath = '../post_segment/';
% path to off files
featurepara.offpath = '../offData/'; 
% point  tsdf  normal  shape, all combined : pointtsdfnormalshape
featurepara.featype ='pointtsdfnormalshape'; 
% whether use the occpation mask to handle cluter
featurepara.removedontcare =1;  
% post pocessing remove empty cell after detection
featurepara.removeEmptybaseOnCell = 0.2; 
% post pocessing remove boxes that are too high 
featurepara.removehigh =0.4; 

%  whether handle occlusion
featurepara.occ =0;

% whether do the local search 
featurepara.localsearch = 0;

%  don't change following 
featurepara.occcellmaxratio =0.0005;
featurepara.occnonlinear =0;
featurepara.weightempty =1;
featurepara.weightfea =1;
featurepara.weightocc =1;
featurepara.pad =0;
featurepara.featurepath_moreneg =''; 
featurepara.normalize =0;
featurepara.codesecondlevel =0;
featurepara.highdim =0;
featurepara.Localview =0;
featurepara.highdimcombine =0;
featurepara.holistic =0;
featurepara.s =0.1;
featurepara.tree =0;
featurepara.kernelsvm = 0;

featurepara.maxdistance =2;
featurepara.mindistance =0;
featurepara.floatcode =1;
featurepara.tril =-1;
featurepara.removefloor =0;
featurepara.camera.ratio =1;