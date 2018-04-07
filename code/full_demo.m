%% set the parameters 
initPath;
className = 'stair';  % className: chair,toilet,bed, table, sofa.
[svmpara,featurepara]=setParameters(); 



ConstsDetection; %change the member in construct,replace box->stair 17.2.26
classNames = consts.classNames;
[~,classId] =ismember(className,classNames);
[idxToimageNum] = split_trainingtestingSetGroup(classId);
replace =1; % whether replace perviouse result
steps.train =1; % do training 
steps.test =0;  % do testing 
steps.visul =0; % generate visulization webpage


%% mode : 
% 1. 'loadmodel' : load the pertrained models and test on the test set.
% 2. 'demo' : train 4 chair example svm and test on the test set. 
% 3. 'useparam': load the parameters to train the full set of models used in paper, 
mode = 'demo_1';


%% set the path 
ExpirmentName = [className '_' mode '_' date]
% output path to trained model
OutputpathTrain = ['../../TrainedModels/' ExpirmentName '/']; 
mkdir(OutputpathTrain);
% output path to testing result for each image
outpathTest =[ '../../Output/MatResult/' ExpirmentName  '/'];
mkdir(outpathTest);
% visualization path 
Weboutputpath =['../../Output/Webpage/'  ExpirmentName  '/'];
mkdir(Weboutputpath);
% save final result boudining box 
path2MatRoot = '../../Output/DetectionMat/';
mkdir(path2MatRoot);

%% 
if strcmp(mode,'loadmodel'),
    OutputpathTrain = ['../perTrainedModels/' className '/'];
    steps.train =0; % skip training 

elseif strcmp(mode,'demo'),
    % example train 4 chairs
    className = 'chair';
    filename = [{'chair_000000358'}];
    Rotation =0.1:0.3:1; % Rotation about gravity direction
    Scale = 1; % height of the model
    Cdepth =[2]; % depth of model
    Cx =[-0.3]; % x position of CG model center
    Cheight =-0.95; % height of CG model center
    Tilt = -0.32; % camera tilt
    shiftAll = [0,0,0]; % slightly shift model within one cell. [sx,sy,sz] <0.1
    aspAll = [1,1,1]; % aspect ratio : always fixed
    [r,s,d,sft,asp,f]=ndgrid([1:length(Rotation)],[1:length(Scale)],[1:length(Cdepth)],[1:size(shiftAll,1)],[1:size(aspAll,1)],[1:length(filename)]);
    x =d;t =x;h =x;
    enlargefloor =0.6*ones(1,length(t(:)));
    enlargeBB =0; 
    enlargewall =0.2;
    
elseif strcmp(mode,'demo_1'),
    %featurepara.featype = 'custom';
    className = 'stair';
    filename = [{'stair'}];
    Rotation = 0.65; % Rotation about gravity direction
    Scale = 0.6; % height of the model
    Cdepth =[2]; % depth of model
    Cx =[0]; % x position of CG model center
    Cheight =-0.95; % height of CG model center
    Tilt = -0.32; % camera tilt
    shiftAll = [0,0,0]; % slightly shift model within one cell. [sx,sy,sz] <0.1
    aspAll = [1,1,1]; % aspect ratio : always fixed
    [r,s,d,sft,asp,f]=ndgrid([1:length(Rotation)],[1:length(Scale)],[1:length(Cdepth)],[1:size(shiftAll,1)],[1:size(aspAll,1)],[1:length(filename)]);
    x =d;t =x;h =x;
    %enlargefloor =0.6*ones(1,length(t(:)));
    enlargefloor = 0;
    enlargeBB =0;
    enlargewall =0;

elseif strcmp(mode,'useparam'),
    load([className '.mat']);
    filename = [{'bed_000000358'}];
    r=1:length(filename);
    s =r; d =r; sft = r; asp =r; f =r; 
    x =d; t =x;h=x;
    enlargeBB =0;
    enlargewall =0.2;
end

test_Para_main;