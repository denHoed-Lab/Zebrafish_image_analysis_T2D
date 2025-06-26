function [LiverData, halfSample] = Liver_process(wellFolder, I, Config)
% LIVER_PROCESS - Liver tissue and lipid analysis using deep learning segmentation
%
% This function performs comprehensive liver analysis combining:
% 1. Deep learning-based liver region segmentation and refinement
% 2. Lipid spot detection and quantification
% 3. Morphological post-processing for accuracy enhancement
% 4. Statistical analysis and result export
%
% The analysis provides four key quantitative measurements:
% - Liver area (segmented and refined)
% - Lipid area (within liver region)
% - Number of lipid droplets
% - Total lipid intensity
%
% Inputs:
%   wellFolder - Path to the well folder containing liver data
%   I - 3D image stack of liver data (multi-channel)
%   Config - Configuration structure with analysis parameters
%
% Outputs:
%   LiverData - Structure containing analysis results and statistics
%   halfSample - Flag indicating if half-sample processing was used
%
% Author: Hanqing Zhang, hanzha@kth.se
% Version: 1.1, 20250619

%% Initialize variables
halfSample = 0; % Flag for half-sample processing (not used in current version)
LiverData = []; % Initialize results structure

% Extract file information for output naming
filesRAW = wellFolder;
[wellFolder, LiverFile_noExt] = fileparts(filesRAW);
[~, well_name, ~] = fileparts(wellFolder);
wellFolder = [wellFolder, filesep];

% Extract experiment information from path
String_segments = split(wellFolder, filesep);
geneName = String_segments{end-3};
experimentName = String_segments{end-2};
wellname = String_segments{end-1};
Liver_subfolder = wellFolder;

%% Channel separation and preprocessing
% Determine number of channels based on configuration
ChannelNum = 2; % Default to 2 channels

if (strcmp(Config.Liver.VAST2Code, '4ch'))
    % 4-channel data processing
    z_sliceNum = size(I, 3) / 4;
    I_bf = I(:, :, 1:z_sliceNum);
    I_bf_save = ImageNorm(max(I_bf, [], 3));
    I_gfp = I(:, :, z_sliceNum + 1:2 * z_sliceNum);
    I_red = I(:, :, 2 * z_sliceNum + 1:3 * z_sliceNum);
    I_red_save = ImageNorm(max(I_red, [], 3));
    I_lipid = I(:, :, 3 * z_sliceNum + 1:4 * z_sliceNum);
    
    % Create maximum projections for analysis
    Ichmax = zeros(size(I, 1), size(I, 2), 2);
    Ichmax(:, :, 1) = max(I_gfp, [], 3); % GFP channel
    Ichmax(:, :, 2) = max(I_lipid, [], 3); % Lipid channel
else
    % 2-channel data processing
    z_sliceNum = size(I, 3) / 2;
    I = reshape(I, [size(I, 1) size(I, 2) z_sliceNum 2]);
    disp(['Found ' num2str(z_sliceNum) ' slices in each channel']);
    
    % Create maximum projections for analysis
    Ichmax = zeros(size(I, 1), size(I, 2), 2);
    I_gfp = I(:, :, 1:z_sliceNum);
    Ichmax(:, :, 1) = max(I(:, :, :, 1), [], 3); % GFP channel
    Ichmax(:, :, 2) = max(I(:, :, :, 2), [], 3); % Lipid channel
    sliceNum = size(I, 3);
end

%% Process each output folder configuration
for i_output = 1:length(Config.Liver.outputFolderlist)
    %% Prepare channel images for segmentation
    % Select appropriate channel based on output folder configuration
    if (strcmp(Config.Liver.outputFolderlist{i_output}, 'max'))
        Ich1max = ImageNorm(Ichmax(:, :, 1)); % Maximum projection
    elseif (strcmp(Config.Liver.outputFolderlist{i_output}, 'middle'))
        I_midSlice = round(sliceNum / 2);
        Ich1max = ImageNorm(I_gfp(:, :, I_midSlice)); % Middle slice
    else
        Ich1max = ImageNorm(Ichmax(:, :, 1)); % Default to maximum projection
    end
    
    Ich2max = Ichmax(:, :, 2); % Lipid channel
    
    % Store original image dimensions
    numrows = size(Ich1max, 1);
    numcols = size(Ich1max, 2);
    
    % Resize images for deep learning segmentation
    Ich1max_Resize = imresize(Ich1max, Config.Liver.Resize);
    Ich2max_Resize = imresize(Ich2max, Config.Liver.Resize);
    
    %% 1. Liver Region Segmentation and Refinement
    % Setup output folder for segmentation results
    outputFolder = strcat(Liver_subfolder, Config.Liver.outputFolderlist{i_output}, filesep);
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
    % Save input image for Python segmentation
    infileMax = strcat(outputFolder, 'LiverCh1maxProj.tif');
    imwrite(Ich1max_Resize, infileMax);
    
    % Execute Python-based liver segmentation
    [pythonSuccess, pythonOutput] = runPythonSegmentation(Config.Liver.pythonLiverSegFile, infileMax, outputFolder, Config);
    if ~pythonSuccess
        warning('Python liver segmentation failed: %s', pythonOutput);
    end
    
    % Load segmentation results
    [~, fname] = fileparts(infileMax);
    try
        f_segm_filename = [outputFolder, fname, '_seg.tif'];
        Im_liverSeg = imread(f_segm_filename);
    catch
        f_segm_filename = [outputFolder, fname, '_seg.png'];
        Im_liverSeg = imread(f_segm_filename);
    end
    
    % Resize segmentation back to original dimensions
    Im_liverSeg = imresize(Im_liverSeg, [numrows numcols]);
    Im_liverSeg = Im_liverSeg > 0;
    
    % Morphological refinement for accuracy enhancement
    Im_liverSeg = imfill(logical(Im_liverSeg), 'holes'); % Fill holes
    Im_liverSeg = imclose(Im_liverSeg, strel('sphere', 3)); % Morphological closing
    
    % Check if liver region was detected
    if (sum(Im_liverSeg(:)) == 0)
        disp('Liver_process: Cannot find single liver region.');
        halfSample = 0;
        LiverData = [];
        return
    end
    
    % Retain largest connected component and apply final refinement
    Im_liverTmp = keepLargestn(Im_liverSeg);
    liverSegArea = imclose(Im_liverTmp, strel('sphere', round(Config.Liver.seKernel)));
    liverSegArea = imfill(logical(liverSegArea), 'holes');
    
    % Check if liver region touches image edges
    touchEdge = 0;
    Im_liverSeg_edgecheck = bwmorph(logical(liverSegArea), 'dilate');
    [nraw, ncol, ~] = size(liverSegArea);
    checkedge = sum(Im_liverSeg_edgecheck(1:nraw, 1)) + sum(Im_liverSeg_edgecheck(1:nraw, ncol)) + ...
        sum(Im_liverSeg_edgecheck(1, 1:ncol)) + sum(Im_liverSeg_edgecheck(nraw, 1:ncol));
    if (checkedge > 0)
        touchEdge = 1;
    end
    
    %% 2. Lipid Segmentation and Analysis
    % Save lipid channel image for Python processing
    infileName2 = strcat(outputFolder, 'LiverCh2maxProj.tiff');
    imwriteMPTiff(uint8(255 * ImageNorm(Ich2max_Resize)), infileName2, 0);
    
    % Execute Python-based lipid detection
    [pythonSuccess, pythonOutput] = runPythonSegmentation(Config.Liver.pythonLiverSegFile2, infileName2, outputFolder, Config);
    if ~pythonSuccess
        warning('Python lipid detection failed: %s', pythonOutput);
    end
    
    % Load lipid segmentation results
    [~, fname2, ~] = fileparts(infileName2);
    f_segm_filename = [outputFolder, fname2, '_seg.png'];
    try
        Im_liverSeg2 = imread(f_segm_filename);
    catch
        disp('Error loading lipid segmentation results!');
        return
    end
    
    % Resize lipid segmentation to original dimensions
    Im_liverSeg2 = imresize(Im_liverSeg2, [numrows numcols], 'nearest');
    LipidPtyhonResults = Im_liverSeg2 == 2; % Extract lipid class (value 2)
    
    % Combine lipid detection with liver segmentation mask
    lipidSpots = LipidPtyhonResults;
    lipidSpots = and(lipidSpots, liverSegArea); % Ensure lipids are only counted within liver region
    
    % Count lipid droplets by finding center points
    lipidSpotsImg = bwmorph(imfill(lipidSpots, 'holes'), 'shrink', inf);
    
    % Calculate intensity statistics
    AverageLipidsBGIntensity_maxproj2 = mean(Ich2max(and(liverSegArea, not(lipidSpots)))); % Background intensity
    AverageLipidsIntensity_maxproj2 = mean(Ich2max(lipidSpots > 0)); % Lipid intensity
    
    %% 3. Visualization and Result Export
    % Create output filename prefix
    Save_name_main = strcat(geneName, '_', experimentName, '_', well_name, '_');
    
    % Create visualization image with lipid spots and liver boundary
    IRGB = imoverlay(imoverlay(uint8(uint8norm(Ich2max)), bwmorph(lipidSpots, 'remove'), 'blue'), imdilate(edge(liverSegArea), ones(3, 3)), 'red');
    imwrite(uint8(uint8norm(IRGB)), [outputFolder, Save_name_main, 'Liver_Demo.png']);
    
    %% 4. Quantification and Results Compilation
    % Compile all quantitative measurements
    LiverData.NrLipids = numel(find(lipidSpotsImg > 0)); % Number of lipid droplets
    LiverData.LiverArea = numel(find(Im_liverSeg)); % Original liver segmentation area
    LiverData.LiverRefinedArea = numel(find(liverSegArea)); % Refined liver segmentation area
    LiverData.LipidArea = sum(lipidSpots(:)); % Total lipid area
    LiverData.AverageLipidsIntensity = AverageLipidsIntensity_maxproj2; % Average lipid intensity
    LiverData.AverageLipidsBGIntensity = AverageLipidsBGIntensity_maxproj2; % Average background intensity
    LiverData.ROI_edge = touchEdge; % Flag indicating if liver touches image edge
    
    %% Save results to Excel file
    try
        filename = [outputFolder, LiverFile_noExt, Config.saveResultSuffix, '.xlsx'];
        writetable(struct2table(LiverData, 'AsArray', true), filename, 'Sheet', 'Statistics', 'WriteMode', 'overwritesheet');
    catch
        filename = [outputFolder, 'Liver', Config.saveResultSuffix, '.xlsx'];
        writetable(struct2table(LiverData, 'AsArray', true), filename, 'Sheet', 'Statistics', 'WriteMode', 'overwritesheet');
    end
end

end

