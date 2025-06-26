function [BetaCellData, Iregion_globalMaxI, BW_3D_segmR_ROI] = betaCells_process(wellFolder, Ibeta, Config)
% betaCells_process - Beta cell analysis using Python + MATLAB 2D and post analysis 
%
% This function performs comprehensive beta cell analysis combining:
% 1. Python-based 3D segmentation for initial region detection
% 2. MATLAB-based 2D segmentation for detailed cell analysis
% 3. Post-processing and statistical analysis
% 4. Result visualization and data export
%
% Inputs:
%   wellFolder - Path to the well folder containing beta cell data
%   Ibeta - 3D image stack of beta cell data
%   Config - Configuration structure with analysis parameters
%
% Outputs:
%   BetaCellData - Structure containing analysis results and statistics
%   Iregion_globalMaxI - Cropped region of interest for analysis
%   BW_3D_segmR_ROI - 3D binary mask of segmented regions
%
% Author: Hanqing Zhang, hanzha@kth.se
% Version: 1.1, 20250619

%% Initialize and extract configuration parameters
% Extract beta cell analysis parameters from configuration
findGlobalMaxGausFilterSize = Config.Betacell.ROISize; % ROI size for analysis
cellsize_DoGKernel = Config.Betacell.cellsize_DoGKernel; % Difference of Gaussian kernel sizes

% Extract additional processing parameters
findRegionFilterSize = Config.Betacell.findRegionFilterSize; % Filter size for region detection
regionBinarizeThreshold = Config.Betacell.regionBinarizeThreshold; % Threshold for region binarization
postProcessWindowSize = Config.Betacell.postProcessWindowSize; % Window size for post-processing morphological operations
smoothKernelDivisor = Config.Betacell.smoothKernelDivisor; % Divisor for smooth kernel calculation
smoothThreshold1 = Config.Betacell.smoothThreshold1; % First smoothing threshold after Gaussian smoothing
watershedDownsizeRatio = Config.Betacell.watershedDownsizeRatio; % Downsize ratio for watershed computational efficiency
finalSmoothThreshold = Config.Betacell.finalSmoothThreshold; % Final smoothing threshold after watershed processing
ROIAdjustmentFactor = Config.Betacell.ROIAdjustmentFactor; % Factor for ROI size adjustment when image is too small

% Extract size property thresholds
% PropStat_Diameter = Config.Betacell.PropStat_Diameter; % Minimum diameter for cell statistics

%% Image preprocessing and format conversion
% Convert image to single precision for processing
Ibeta = single(Ibeta); % Ensure 16-bit input data is properly handled

% Extract file information for output naming
[wellFolder, betaCellFile] = fileparts(wellFolder);
wellFolder = [wellFolder, filesep];
[~, betaCellFile_noExt] = fileparts(betaCellFile);

%% Region of Interest (ROI) Detection
% Create maximum intensity projection for ROI detection
Ibeta_max = max(Ibeta, [], 3);

% Create output directory if it doesn't exist
if ~exist([wellFolder, Config.Betacell.outputFolderlist{1}, filesep], 'dir')
    mkdir([wellFolder, Config.Betacell.outputFolderlist{1}, filesep]);
end

% Adjust ROI size if image is too small
s_fix = min(size(Ibeta, 1), size(Ibeta, 2)); % Get minimum dimension
if (s_fix < findGlobalMaxGausFilterSize)
    findGlobalMaxGausFilterSize = round(s_fix / ROIAdjustmentFactor) + 1;
end

%% Find regions of interest using intensity-based detection
% Apply Gaussian filtering to smooth the maximum projection
Ibeta_max_filt = ImageNorm(imfilter(Ibeta_max, ones(findRegionFilterSize, findRegionFilterSize)));

% Binarize image to find candidate regions
Region_cand = imbinarize(Ibeta_max_filt, regionBinarizeThreshold);

% Apply morphological closing to connect nearby regions
Region_cand = imclose(Region_cand, ones(findRegionFilterSize, findRegionFilterSize));

% Extract region properties (centroid, area, mean intensity)
Region_cand_props = regionprops(Region_cand, Ibeta_max_filt, 'WeightedCentroid', 'Area', 'MeanIntensity');

%% Determine the best region for analysis
MultiRegions = 1; % Default to single region

if (length(Region_cand_props) > 1)
    % Multiple regions found - select the one with highest intensity-area product
    MultiRegions = length(Region_cand_props);
    max_v = 0;
    
    for i = 1:length(Region_cand_props)
        max_cand = Region_cand_props(i).Area * Region_cand_props(i).MeanIntensity;
        if (max_cand > max_v)
            max_v = max_cand;
            r_col = round(Region_cand_props(i).WeightedCentroid(1));
            r_row = round(Region_cand_props(i).WeightedCentroid(2));
        end
    end
    
elseif (length(Region_cand_props) == 1)
    % Single region found
    r_col = round(Region_cand_props(1).WeightedCentroid(1));
    r_row = round(Region_cand_props(1).WeightedCentroid(2));
    
else
    % No regions found - save diagnostic image and return empty results
    imSave = uint8(uint8norm(Ibeta_max));
    imwrite(imSave, [wellFolder, Config.Betacell.outputFolderlist{1}, filesep, betaCellFile_noExt, '_maxproj.jpg'])
    disp('Cannot find Betacells region');
    BetaCellData = [];
    Iregion_globalMaxI = [];
    BW_3D_segmR_ROI = [];
    return
end

%% Extract Region of Interest (ROI) with boundary checking
% Calculate ROI boundaries
rows = [r_row - round(findGlobalMaxGausFilterSize / 2):r_row + round(findGlobalMaxGausFilterSize / 2) - 1];
cols = [r_col - round(findGlobalMaxGausFilterSize / 2):r_col + round(findGlobalMaxGausFilterSize / 2) - 1];

% Handle boundary conditions for rows
roi_edge = 0; % Flag for edge cases
if rows(1) < 1
    shiftx = 1 - rows(1);
    rows = rows + shiftx;
    roi_edge = 1;
end
if rows(end) > size(Ibeta_max, 1)
    shiftx = rows(end) - size(Ibeta_max, 1);
    rows = rows - shiftx;
    roi_edge = 1;
end

% Handle boundary conditions for columns
if cols(1) < 1
    shifty = 1 - cols(1);
    cols = cols + shifty;
    roi_edge = 1;
end
if cols(end) > size(Ibeta_max, 2)
    shifty = cols(end) - size(Ibeta_max, 2);
    cols = cols - shifty;
    roi_edge = 1;
end

%% Save ROI visualization
% Create visualization image with ROI rectangle
imSave = uint8(uint8norm(Ibeta_max));
imRect = zeros(size(imSave));
imRect(rows, cols) = 1;
imRect = imdilate(edge(imRect), ones(2, 2)); % Create rectangle outline
imSave(imRect == 1) = 255; % Draw white rectangle
imwrite(imSave, [wellFolder, Config.Betacell.outputFolderlist{1}, filesep, betaCellFile_noExt, '_maxproj.jpg']);

%% Extract ROI for analysis
Iregion_globalMaxI = Ibeta(rows, cols, :);
disp('Image Segmentation...');

%% Python-based 3D Segmentation
% Setup output folder for Python segmentation
outputFolder = strcat(wellFolder, Config.Betacell.outputFolderlist{1}, filesep);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Save ROI image for Python processing
infileName = strcat(wellFolder, Config.Betacell.outputFolderlist{1}, filesep, betaCellFile_noExt, '.tiff');
imwriteMPTiff(uint8(uint8norm(Iregion_globalMaxI)), infileName, 0);

% Execute Python segmentation using wrapper function
[pythonSuccess, pythonOutput] = runPythonSegmentation(Config.Betacell.pythonVascSegFile, infileName, outputFolder, Config);
if ~pythonSuccess
    error('Python segmentation failed: %s', pythonOutput);
end

% Load Python segmentation results
f_segm_filename = [outputFolder, betaCellFile_noExt, '_bw.tif'];
Im_liverSeg = imreadMPTiff(f_segm_filename, 1);

%% MATLAB-based 2D Segmentation
% Initialize 2D segmentation array
ImSeg2D = zeros(size(Iregion_globalMaxI));

% Process each slice with MATLAB network
for i_stack = 1:size(Iregion_globalMaxI, 3)
    % Resize slice to network input size
    Itmp = imresize(Iregion_globalMaxI(:, :, i_stack), Config.Betacell.Deeplabv3.Network.ImageSize);
    
    % Perform semantic segmentation using MATLAB network
    pxdsResults = semanticseg(Itmp, Config.Betacell.MATLAB.net);
    
    % Resize results back to original size and create binary mask
    ImSeg2D(:, :, i_stack) = imresize((pxdsResults ~= Config.Betacell.Deeplabv3.classNames(2)), [size(Iregion_globalMaxI, 1) size(Iregion_globalMaxI, 2)]);
end

% Create binary masks for further processing
BW_3D_segmR_ROI = logical(Im_liverSeg); % Python 3D segmentation result
bw_2D_SegmR = logical(ImSeg2D); % MATLAB 2D segmentation result

%% Post-processing and refinement
disp('Post-processing...');
% Smooth ROI using morphological 
windowSize = postProcessWindowSize;
se2 = strel('sphere', windowSize);
BW_3D_segmR_ROI = imclose(BW_3D_segmR_ROI, se2);
bw_2D_SegmR = imclose(bw_2D_SegmR, se2);

% Combine Python and MATLAB segmentation results
bw_2D_SegmR = and(BW_3D_segmR_ROI, bw_2D_SegmR);

%% Isotropic voxel analysis for accurate measurements
% Resize to isotropic voxels based on resolution
bw_2D_SegmR = logical(imresize3(single(bw_2D_SegmR), round([size(bw_2D_SegmR) .* Config.resolutionXYZ / (Config.resolutionXYZ(1))]), 'method', 'nearest'));
Iregion_globalMaxI = imresize3(Iregion_globalMaxI, round([size(Iregion_globalMaxI) .* Config.resolutionXYZ / (Config.resolutionXYZ(1))]));

% Label connected components
bw_labeled = bwlabeln(bw_2D_SegmR > 0);
bw_2D_SegmR = bw_labeled > 0;

%% Smooth cell contours for better shape analysis
windowSize = postProcessWindowSize;
se1 = strel('sphere', windowSize);
kernel = windowSize / smoothKernelDivisor;

% Apply 3D Gaussian smoothing
bw_2D_SegmR_smoothGray = smooth3(single(bw_2D_SegmR), 'guassian', [windowSize, windowSize, windowSize], kernel);
bw_tmp = bw_2D_SegmR_smoothGray > smoothThreshold1;
bw_2D_SegmR = imopen(bw_tmp, se1);

%% Watershed segmentation for cell separation
% Apply Difference of Gaussian (DoG) filter for intensity-based separation
Ibeta_crop_doG = imgaussian(Iregion_globalMaxI, cellsize_DoGKernel(1)) - imgaussian(Iregion_globalMaxI, cellsize_DoGKernel(2) + 1);
Ibeta_crop_doG(Ibeta_crop_doG < 0) = 0; % Remove negative values

% Resize for computational efficiency
downsize_ratio = watershedDownsizeRatio;
bw_2D_SegmR_resize = imresize3(bw_2D_SegmR, downsize_ratio, 'nearest');
DistanceTransCellShape_resize = bwdist(~bw_2D_SegmR_resize); % Distance transform
Ibeta_watershed_matrix = imresize3(Ibeta_crop_doG, downsize_ratio);
Ibeta_watershed_matrix =-(DistanceTransCellShape_resize .* Ibeta_watershed_matrix);
% Apply watershed using both shape and intensity information
watershed_D_resize = watershed(Ibeta_watershed_matrix);
watershed_D_resize(~bw_2D_SegmR_resize) = 0;
bw_final = watershed_D_resize > 0;

% Remove small objects and resize back
bw_final = imopen(bw_final, se2);
bw_labeled = bwlabeln(bw_final);
bw_labeled = imresize3(bw_labeled, [size(bw_2D_SegmR)], 'nearest');

% Final smoothing after resizing
bw_2D_SegmR_smoothGray = smooth3(single(bw_labeled > 0), 'guassian', [3, 3, 3], 0.6);
bw_2D_SegmR = bw_2D_SegmR_smoothGray > finalSmoothThreshold;
bw_labeled = bwlabeln(bw_2D_SegmR > 0);

%% Save final post-processed binary mask (replacing Python output)
% Create binary mask with cell values = 255
final_binary_mask = uint8(bw_2D_SegmR) * 255;

% Resize back to original TIFF dimensions
original_size = size(Im_liverSeg); % Get original size from Python output
final_binary_mask_resized = imresize3(final_binary_mask, original_size, 'nearest');

% Save the final post-processed binary mask
final_mask_filename = [outputFolder, betaCellFile_noExt, '_bw.tif'];
imwriteMPTiff(final_binary_mask_resized, final_mask_filename, 0);

%% Statistical analysis of detected cells
% Extract 3D region properties for each cell
cellProp = regionprops3(bw_labeled, Iregion_globalMaxI, 'Volume', 'PrincipalAxisLength', 'SurfaceArea', 'EquivDiameter', 'MeanIntensity', 'MaxIntensity', 'MinIntensity', 'Centroid', 'WeightedCentroid', 'ConvexVolume');

%% Islet-level analysis using alpha shape
if (~isempty(cellProp))
    % Create alpha shape from cell centroids
    shp = alphaShape(cellProp.Centroid);
    
    if (shp.Alpha ~= Inf)
        % Generate islet binary mask from alpha shape
        islet_BW = alphaShape2im_fast(shp, size(bw_labeled) / 2, size(bw_labeled));
        
        % Calculate islet properties
        isletProp = regionprops3(or(islet_BW, bw_2D_SegmR), 'ConvexImage', 'ConvexVolume', 'SurfaceArea', 'BoundingBox', 'EquivDiameter', 'PrincipalAxisLength');
        
        % Remove single-voxel artifacts
        islet_delete = isletProp.ConvexVolume == 1;
        isletProp(islet_delete, :) = [];
        
        if (~isempty(isletProp))
            % Extract islet measurements
            isletVolume = isletProp.ConvexVolume(1);
            isletSurfaceArea = isletProp.SurfaceArea(1);
            isletEquivDiameter = isletProp.EquivDiameter(1);
            isletPrincipleAxisLength = isletProp.PrincipalAxisLength(1, :);
        end
    else
        % Invalid alpha shape - set default values
        isletVolume = 0;
        isletSurfaceArea = 0;
        isletEquivDiameter = 0;
        isletPrincipleAxisLength = 0;
    end
else
    % No cells detected - set default values
    isletVolume = 0;
    isletSurfaceArea = 0;
    isletEquivDiameter = 0;
    isletPrincipleAxisLength = 0;
end

%% Compile analysis results
cellNumBefore = size(cellProp, 1);

% Initialize result structure with default values
BetaCellData.Volume = numel(find(bw_2D_SegmR > 0));
BetaCellData.CellDiameter_Average = 0;
BetaCellData.CellDiameter_Stdev = 0;
BetaCellData.CellVolume_Average = 0;
BetaCellData.CellVolume_Stdev = 0;
BetaCellData.CellSurfaceArea_Average = 0;
BetaCellData.CellSurfaceArea_Stdev = 0;
BetaCellData.CellIntensity_Average = 0;
BetaCellData.CellIntensity_Stdev = 0;

%% Calculate cell-level statistics
if (~isempty(cellProp) && size(cellProp, 1) >= 1)
    cellNum = size(cellProp, 1);
    
    % Initialize arrays for cell measurements
    CellAverageDiameter = zeros(cellNum, 1);
    CellAverageVolume = zeros(cellNum, 1);
    CellAverageSurfaceArea = zeros(cellNum, 1);
    CellAverageIntensity = zeros(cellNum, 1);
    
    % Extract measurements for each cell
    for i = 1:cellNum
        CellAverageDiameter(i) = cellProp.EquivDiameter(i);
        CellAverageVolume(i) = cellProp.Volume(i);
        CellAverageSurfaceArea(i) = cellProp.SurfaceArea(i);
        CellAverageIntensity(i) = cellProp.MeanIntensity(i);
    end
    
    % Calculate mean values
    BetaCellData.CellDiameter_Average = round(mean(CellAverageDiameter(:)));
    BetaCellData.CellVolume_Average = round(mean(CellAverageVolume(:)));
    BetaCellData.CellSurfaceArea_Average = round(mean(CellAverageSurfaceArea(:)));
    BetaCellData.CellIntensity_Average = round(mean(CellAverageIntensity(:)));
    
    % Calculate standard deviations
    BetaCellData.CellDiameter_Stdev = round(std(CellAverageDiameter(:)));
    BetaCellData.CellVolume_Stdev = round(std(CellAverageVolume(:)));
    BetaCellData.CellSurfaceArea_Stdev = round(std(CellAverageSurfaceArea(:)));
    BetaCellData.CellIntensity_Stdev = round(std(CellAverageIntensity(:)));
else
    cellNum = 0;
end

%% Compile final results
BetaCellData.NrCells = cellNum;
BetaCellData.Total_Cell_Intensity = sum(Iregion_globalMaxI(bw_2D_SegmR));
BetaCellData.AverageIntensity_Mean = round(mean(Iregion_globalMaxI(:)));
BetaCellData.AverageIntensity_Stdev = round(std(Iregion_globalMaxI(:)));
BetaCellData.AverageIntensityMasked_Mean = round(mean(Iregion_globalMaxI(bw_2D_SegmR)));
BetaCellData.AverageIntensityMasked_Stdev = round(std(Iregion_globalMaxI(bw_2D_SegmR)));

% Islet-level measurements
BetaCellData.isletNr = MultiRegions;
BetaCellData.isletVolume = round(isletVolume);
BetaCellData.isletSurfaceArea = round(isletSurfaceArea);
BetaCellData.isletEquivDiameter = round(isletEquivDiameter);
BetaCellData.isletPrincipleAxisLength = round(isletPrincipleAxisLength);
BetaCellData.ROI_edge = roi_edge;

%% Save results to Excel file
filename = [wellFolder, Config.Betacell.outputFolderlist{1}, filesep, betaCellFile_noExt, Config.saveResultSuffix, '.xlsx'];

% Convert results to table format
cellPropAll = struct2table(BetaCellData);

% Prepare cell property table with expanded coordinate columns
cellVarNames = cellProp.Properties.VariableNames;
cellVarNames = {cellVarNames{1} strcat(cellVarNames{2}, '_x') strcat(cellVarNames{2}, '_y') strcat(cellVarNames{2}, '_z') ...
    cellVarNames{3} strcat(cellVarNames{4}, '_x') strcat(cellVarNames{4}, '_y') strcat(cellVarNames{4}, '_z') ...
    cellVarNames{5:6} strcat(cellVarNames{7}, '_x') strcat(cellVarNames{7}, '_y') strcat(cellVarNames{7}, '_z') ...
    cellVarNames{8:end}};

% Convert cell properties to array format
CcellProp = table2cell(cellProp);
cellTableArray = [];

try
    % Standard conversion for all data types
    for i = 1:size(CcellProp, 1)
        cellTableArray(i, :) = round([CcellProp{i, 1:end}]);
    end
catch
    % Handle special cases with cell arrays
    for i = 1:size(CcellProp, 1)
        cellTableArray(i, :) = round([CcellProp{i, 1:end-2} cellProp{i, end-1}{1} cellProp{i, end}{1}]);
    end
end

% Create cell property table
if (~isempty(cellTableArray))
    cellPropR = array2table(cellTableArray, 'VariableNames', cellVarNames);
else
    cellPropR = table();
end

% Delete existing file and write new results
if isfile(filename)
    delete(filename);
end
% Write results to Excel file with multiple sheets
writetable(cellPropAll, filename, 'Sheet', 'Statistics', 'WriteMode', 'overwritesheet');
writetable(cellPropR, filename, 'Sheet', 'Cell Properties', 'WriteMode', 'overwritesheet');
disp('Done!');

%% Create visualization results
% Create labeled image with cell boundaries
bw_labeledwrite = zeros(size(bw_labeled));
for i = 1:size(Iregion_globalMaxI, 3)
    bw_labeledwrite(:, :, i) = uint16(bw_labeled(:, :, i)) .* uint16(imdilate(edge(bw_2D_SegmR(:, :, i)), ones(5, 5)));
end

% Create RGB visualization
if (max(bw_labeledwrite(:)) > 0)
    % Convert labels to RGB colors
    RGB = label2rgb3d(bw_labeledwrite, 'prism', [0 0 0], 'shuffle');
    bw_notedge = bw_labeledwrite <= 1;
    RGB = permute(RGB, [1 2 4 3]);
    RGB_original = zeros(size(RGB));
    
    % Combine with original intensity image
    Iregion_globalMaxI_norm = ImageNorm(Iregion_globalMaxI);
    for r_col = 1:size(RGB, 4)
        RGB(:, :, :, r_col) = RGB(:, :, :, r_col) + repmat(Iregion_globalMaxI_norm(:, :, r_col), [1 1 3]) .* double(repmat(bw_notedge(:, :, r_col), [1 1 3]));
        RGB_original(:, :, :, r_col) = repmat(Iregion_globalMaxI_norm(:, :, r_col), [1 1 3]);
    end
    RGB = permute(RGB, [1 2 4 3]);
    RGB_original = permute(RGB_original, [1 2 4 3]);
else
    % No cells detected - use original image
    RGB = repmat(ImageNorm(Iregion_globalMaxI), [1 1 1 3]);
    RGB_original = repmat(ImageNorm(Iregion_globalMaxI), [1 1 1 3]);
end
%% Save visualization results as video
RGB_all = [RGB_original RGB];
v = VideoWriter([wellFolder, Config.Betacell.outputFolderlist{1}, filesep, betaCellFile_noExt, '_res'], 'Motion JPEG AVI');
open(v);
for k = 1:size(RGB, 3)
    writeVideo(v, squeeze(uint8(255 * RGB_all(:, :, k, :))));
end
close(v);

end