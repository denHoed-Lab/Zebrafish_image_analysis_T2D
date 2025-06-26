function Config = Config_Demo()
% CONFIG_DEMO - Configuration file for VAST image analysis pipeline
%
% This function returns the configuration structure with all analysis parameters
% organized by analysis type: Vasc → Betacells → Liver → VAST → General
%
% Author: Hanqing Zhang, hanzha@kth.se
% Version: 1.1, 20250619

%% General Analysis Parameters (Base configuration)
Config.ReadfromNAS = 0;
Config.ChannelCode = '2ch'; % mpomeg, ik17 4ch, ins 2ch, 3ch for gene 52,54
Config.VAST2Code = Config.ChannelCode; % Set VAST2Code to match ChannelCode

%% 1. VASCULAR ANALYSIS PARAMETERS
Config.Vasc.ChannelCode = Config.ChannelCode; % Inherit from general config
Config.Vasc.NeutrophilDoG1 = [3 6];
Config.Vasc.NeutrophilDoG2 = [1.5 3];
Config.Vasc.NeutrophilSizeMin = 15;
Config.Vasc.NeutrophilSizeMax = 150;
Config.Vasc.NeutrophilThresh = 0.15;
Config.Vasc.MacrophageDoG1 = [3 8];
Config.Vasc.MacrophageDoG2 = [1.5 4];
Config.Vasc.MacrophageSizeMin = 18;
Config.Vasc.MacrophageSizeMax = 180;
Config.Vasc.MacrophageThresh = 0.18;
Config.Vasc.ROI_OpenkSize = 35;
Config.Vasc.Operation = 'DL_2ch';  % DL 2ch for 
                                   % N4ch(4ch), input (2): maximum projection and slice 25
                                   % BGSub4(4ch),input (2): mp with bgsub gfp, slice 25 with bgsub 25 gfp
                                   % Leakage (4ch),input (3): mp, 25 and gfp
                                   % Leakage2 (4ch),input (2): mp, 25 
                                   % Eug_2ch
                                   % Collect, collect results
Config.Vasc.outputFolderlist = {'Vasc_lipid_leakage_v4'};
Config.Vasc.Resize = [512 512];

%% 2. BETA CELL ANALYSIS PARAMETERS
Config.Betacell.MatlabBetacellsSegFile = '.\DNNs\MATLAB_weights\Betacells2D_resnet50_v100.mat';
Config.Betacell.Deeplabv3.Network.ImageSize = [512 512];
Config.Betacell.Deeplabv3.classNames = ["Betacells", "background"];

% Load beta cell network
LoadNet = load(Config.Betacell.MatlabBetacellsSegFile);
if isfield(LoadNet, 'MyNet')
    disp('Loading beta cell network...')
    net_ = layerGraph(LoadNet.MyNet);
end
Config.Betacell.MATLAB.net = assembleNetwork(net_);

Config.Betacell.pythonVascSegFile = '.\DNNs\python\segment_beta_cells.py';
Config.Betacell.Operation = 'DL_Betacell_vast2';
Config.Betacell.Debug = 0;
Config.Betacell.BetaCellExpr1 = 'betacell';
Config.Betacell.BetaCellExpr2 = 'lng';
Config.Betacell.SaveTrainingdata = 0; % Change to 1 if training data should be saved
Config.Betacell.outputFolderlist = {'LocalResults\Betacell'};
Config.Betacell.ROISize = 600; % VAST2-600, VAST1-300
Config.Betacell.MorphOpenKernelSize = 500; % VAST2-500(22*22 obj) VAST1-5
Config.Betacell.cellRegion_MaxProjBlurKernel = 20; % VAST2-20
Config.Betacell.cellRegion_dilateKernelDiameter = 100; % VAST2-100
Config.Betacell.cellsize_DoGKernel = [7 13]; % VAST2-[7 12]
Config.Betacell.cellsize_minBetaCellVolume = 1000; % VAST2-1000
Config.Betacell.PropStat_Diameter = 15; % VAST2-15
Config.Betacell.PropStat_surfaceArea = 100; % VAST2-100
Config.Betacell.PropStat_MaxIntensity = 2000; % VAST2-2000

% Post-analysis parameters for betaCells_process
Config.Betacell.findRegionFilterSize = 50; % Filter size for region detection
Config.Betacell.regionBinarizeThreshold = 0.5; % Threshold for region binarization
Config.Betacell.postProcessWindowSize = 5; % Window size for post-processing morphological operations
Config.Betacell.smoothKernelDivisor = 1.4; % Divisor for smooth kernel calculation (windowSize/divisor)
Config.Betacell.smoothThreshold1 = 0.5; % First smoothing threshold after Gaussian smoothing
Config.Betacell.watershedDownsizeRatio = 0.6; % Downsize ratio for watershed computational efficiency
Config.Betacell.finalSmoothThreshold = 0.6; % Final smoothing threshold after watershed processing
Config.Betacell.ROIAdjustmentFactor = 1.2; % Factor for ROI size adjustment when image is too small

% Beta cell system parameters
Config.Betacell.VASTsys = 0; % 1 Resize twice for lipid detection in using old vast data

%% 3. LIVER ANALYSIS PARAMETERS
Config.Liver.VAST2Code = Config.ChannelCode; % Inherit from general config
Config.Liver.Operation = 'DL_3ch'; % '2ch' , '4ch', 'Collect' DL_3ch
Config.Liver.VAST1Code = 'Ch2LipidGFP'; % ! Normal() Cross2(PAM,RREB1) Ch2LipidGFP(Christoph_new210701,Eugenia) default : 2ch zstack data in sequence
Config.Liver.ChannelCode = Config.ChannelCode; % mpomeg 4ch, other 2ch
Config.Liver.outputFolderlist = {'LocalResults\Liver'};
Config.Liver.saveResultSuffix = {''}; % {'','_neutrophils','_macrophages'}
Config.Liver.pythonLiverSegFile = '.\DNNs\python\LiverSegmentation.py'; % Liver segmentation
Config.Liver.pythonLiverSegFile2 = '.\DNNs\python\LiverLipid_1ch_segmentation.py'; % Lipid detection

% Inherit neutrophil parameters from Vasc config
Config.Liver.NeutrophilDoG1 = Config.Vasc.NeutrophilDoG1;
Config.Liver.NeutrophilDoG2 = Config.Vasc.NeutrophilDoG2;
Config.Liver.NeutrophilSizeMin = Config.Vasc.NeutrophilSizeMin;
Config.Liver.NeutrophilSizeMax = Config.Vasc.NeutrophilSizeMax;
Config.Liver.NeutrophilThresh = 0.15;

% Inherit macrophage parameters from Vasc config
Config.Liver.MacrophageDoG1 = Config.Vasc.MacrophageDoG1;
Config.Liver.MacrophageDoG2 = Config.Vasc.MacrophageDoG2;
Config.Liver.MacrophageSizeMin = Config.Vasc.MacrophageSizeMin;
Config.Liver.MacrophageSizeMax = Config.Vasc.MacrophageSizeMax;
Config.Liver.MacrophageThresh = 0.18;

Config.Liver.ROI_OpenkSize = 51;
Config.Liver.Resize = [512 512]; % 2
Config.Liver.LipidBG_Factor = 3; % 3
Config.Liver.localmaxima_bias = 10; % 10
Config.Liver.relativeInt_bias = 1.5; % 1.5
Config.Liver.seKernel = 20; % Morphological kernel size for liver segmentation refinement

% File expression patterns for liver analysis
Config.Liver.Expr1 = 'liver'; % First expression for liver file detection
Config.Liver.Expr2 = 'lng'; % Second expression for liver file detection

% Liver-specific file processing parameters
Config.Liver.useThunder = 0; % ! 1 for old data from vast 1

%% 4. VAST ANALYSIS PARAMETERS
Config.VAST.skipSegm = 0;
Config.VAST.outputFolderlist = {'LocalResults\VAST'};
Config.VAST.pythonFile = '.\DNNs\python\FishAnalysis_Image.py';
Config.VAST.useVAST_views = {'_1_1', '_1_4', '_1_7', '_1_10'};
Config.VAST.classNames{1} = ["body", "background"]; % whole body
Config.VAST.classNames{2} = ["sb", "background"]; % swim bladder
Config.VAST.classNames{3} = ["eye", "ear", "back", "heart", "background"]; % eye,ear,backline,pericardium regions
Config.VAST.labelIDs{1} = [255 0];
Config.VAST.labelIDs{2} = [255 0];
Config.VAST.labelIDs{3} = [4 3 2 1 0];
Config.VAST.ImageSize = [512 512];
Config.VAST.LoadTrainedNetFileName{1} = 'VAST_allViews.mat';
Config.VAST.LoadTrainedNetFileName{2} = 'VAST_dorsal.mat';
Config.VAST.LoadTrainedNetFileName{3} = 'VAST_Lateral.mat';
Config.VAST.TrainedNetFilePath = '.\DNNs\MATLAB_weights\';

% File naming and output parameters
Config.VAST.PythonSegmOutputformat = '.tif';
Config.VAST.VASTFoldersNFiles = 'VAST\*';
Config.VAST.VASToutputDetectsufix = '_seg*';
Config.VAST.AnalysisVASTResultsSuffix = '_resultVAST.csv';
Config.VAST.SaveAnalysisResultName = 'VAST_results.xlsx';

% Conversion and measurement parameters
Config.VAST.umConv = 5600/1024; % conversion from pixel to micrometer

% Border detection parameters
Config.VAST.BorderDetectionThreshold = 3; % pixels for border detection (>3 pixels into image edges)

% Tail position analysis parameters
Config.VAST.TailAnalysisFixedRegion = 10; % fixed region size for long fish
Config.VAST.TailAnalysisAdaptiveFactor = 0.5; % factor for adaptive region (half length for short fish)

% Eye morphological operations
Config.VAST.EyeMorphDiskRadius = 3; % radius for disk-shaped structuring element

% Backline analysis parameters
Config.VAST.BacklineMorphDiskRadius = 15; % radius for backline morphological close operation
Config.VAST.BacklineThinIterations = Inf; % iterations for thinning operation

% Pigment analysis parameters
Config.VAST.PigmentThreshold = 0.7; % threshold for enhanced pigment analysis (70%)

% Developmental analysis parameters
Config.VAST.DevelopmentalAngleConversion = 180/pi; % conversion factor for angle calculation

% Load VAST networks
for net_i = 1:length(Config.VAST.classNames)
    LoadNet = load(strcat(Config.VAST.TrainedNetFilePath, Config.VAST.LoadTrainedNetFileName{net_i}));
    if isfield(LoadNet, 'MyNet')
        disp('Loading VAST network...')
        net_ = layerGraph(LoadNet.MyNet);
    end
    Config.VAST.net{net_i} = assembleNetwork(net_);
end

%% 5. PYTHON ENVIRONMENT MANAGEMENT (Global for all Python scripts)
Config.PythonEnv.AutoSetup = 1; % 1: auto-setup Python env, 0: use existing
Config.PythonEnv.EnvName = 'pipeline_python_env'; % Name of conda environment
Config.PythonEnv.PythonVersion = '3.9'; % Python version to use
Config.PythonEnv.RequirementsFile = '.\DNNs\python\requirements.txt'; % Path to requirements file
Config.PythonEnv.ForceReinstall = 0; % 1: force reinstall packages, 0: skip if already installed
Config.PythonEnv.TestMode = 0; % 1: run in test mode with verbose output, 0: normal mode
Config.PythonEnv.FallbackPath = ''; % Fallback Python path if auto-setup fails
Config.PythonEnv.LogFile = 'python_env_setup.log'; % Log file for setup process

%% 6. GENERAL PARAMETERS (Used across multiple analysis types)
% File processing parameters
Config.skipOverwrite = 1; % 0 to overwrite begin from the missing data, 1 to write missing data only
Config.overwriteResults = 0; % 0: skip existing results, 1: overwrite existing results

% System parameters
Config.saveResultSuffix = '_Results';

% Resolution and image processing parameters
Config.resolutionXYZ = [0.419, 0.419, 1.143]; % micrometer, z stack beta cells data
Config.outputResolution = 0.419; % this is fixed to xy plane resolution to get correct number of cells
Config.resolutionTarget = [Config.outputResolution, Config.outputResolution, Config.outputResolution]; % for iso-resolution

end 