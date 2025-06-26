function VAST_process(inputFolder, Config)
% VAST_PROCESS - VAST brightfield image analysis using deep learning segmentation
%
% This function performs comprehensive fish analysis using 4 views out of 12 VAST views:
% - 2 dorsal views (#1, #7) and 2 lateral views (#4, #10)
% - Uses DeepLab with ResNet backbone trained on local VAST datasets
% - Provides morphological measurements, border detection, and developmental analysis
%
% Analysis Components:
% 1. Body Segmentation for all views (MaxFeret diameter, area, border detection, tail position)
% 2. Dorsal View Analysis (swim bladder detection, pigment region analysis)
% 3. Lateral View Analysis (eye, ear, pericardium, backline with curvature analysis)
%
% Inputs:
%   inputFolder - Path to the experiment folder containing well subdirectories
%   Config - Configuration structure with all VAST analysis parameters
%
% Outputs:
%   Saves individual well results as CSV files and aggregated results as Excel
%
% Author: Hanqing Zhang, hanzha@kth.se
% Version: 1.1, 20250619

%% Initialize analysis parameters
% Use config-defined parameters instead of hardcoded values
umConv = Config.VAST.umConv; % conversion from pixel to micrometer

% Output file naming conventions
PythonSegmOutputformat = Config.VAST.PythonSegmOutputformat;
VASTFoldersNFiles = Config.VAST.VASTFoldersNFiles;
VASToutputDetectsufix = Config.VAST.VASToutputDetectsufix;
AnalysisresultsfolderPathName = Config.VAST.outputFolderlist{1}; % Use config-defined output folder
AnalysisVASTResultsSuffix = Config.VAST.AnalysisVASTResultsSuffix;
SaveAnalysisResultName = Config.VAST.SaveAnalysisResultName;

% Setup folder paths
saveRootFolder = inputFolder;
inputExperimentFolder = inputFolder;

% Get experiment information
sourceWellFolders = dir(inputExperimentFolder);
[~, experimentName, ~] = fileparts(inputExperimentFolder);

% Filter for well directories only
sourceWellFolders = sourceWellFolders(~ismember({sourceWellFolders(:).name}, {'.', '..'}));
sourceWellFolders = sourceWellFolders([sourceWellFolders(:).isdir]);

%% Process each well if segmentation is enabled
if (Config.VAST.skipSegm == 0)
    resultData = struct();
    
    for i = 1:numel(sourceWellFolders)
        % Initialize result data structure for current well
        resultData.ExperimentName = strcat(experimentName, '_', sourceWellFolders(i).name);
        resultData.Wellfolder = sourceWellFolders(i).name;
        
        % Determine source and save folders based on NAS configuration
        if (Config.ReadfromNAS == 1)
            String_segments = split(sourceWellFolders(i).folder, filesep);
            geneName = String_segments{end-1};
            experimentName = String_segments{end};
            wellname = sourceWellFolders(i).name;
            sourceFolder = [sourceWellFolders(i).folder, filesep, sourceWellFolders(i).name, filesep];
            saveFolder = strcat(Config.localNASrootPath, filesep, geneName, filesep, experimentName, filesep, wellname, filesep);
        else
            String_segments = split(sourceWellFolders(i).folder, filesep);
            geneName = String_segments{end-1};
            experimentName = String_segments{end};
            wellname = sourceWellFolders(i).name;
            sourceFolder = [sourceWellFolders(i).folder, filesep, sourceWellFolders(i).name, filesep];
            saveFolder = sourceFolder;
        end
        
        welldate = 'NaT'; % Default datetime value
        Write2FileVAST = 1; % Flag to control file writing
        
        %% Process each VAST view (4 views: 2 dorsal + 2 lateral)
        for VASTview = 1:numel(Config.VAST.useVAST_views)
            % Search for VAST images for current view
            VASTfiles = dir([sourceFolder, VASTFoldersNFiles, Config.VAST.useVAST_views{VASTview}, '.*']);
            
            if (~isempty(VASTfiles))
                % Extract file date for metadata
                try
                    welldate = VASTfiles(1).date;
                catch
                    welldate = 'NaT';
                end
                
                % Create results folder
                resultFolder = [saveFolder, AnalysisresultsfolderPathName, filesep];
                if (~isfolder(resultFolder))
                    mkdir(resultFolder);
                end
                
                % Handle file naming issues (spaces in filenames)
                VASTFile = [sourceFolder, 'VAST', filesep, VASTfiles(1).name];
                [path, name, n_ext] = fileparts(VASTFile);
                new_name = strrep(name, ' ', '_'); % Fix badly named files with spaces
                VASTFile_new = strcat(path, filesep, new_name, n_ext);
                
                if (~strcmp(VASTFile_new, VASTFile))
                    movefile(VASTFile, VASTFile_new, 'f'); % Replace old file
                    VASTfiles(1).name = strcat(new_name, n_ext); % Update filename
                else
                    VASTFile_new = VASTFile;
                end
                
                % Check if segmentation already exists
                [~, vast_Img_name, ~] = fileparts(VASTFile_new);
                resultFile_precheck = dir([saveFolder, filesep, AnalysisresultsfolderPathName, filesep, vast_Img_name, VASToutputDetectsufix]);
                resultFile_precheck2 = isfile([saveRootFolder, filesep, sourceWellFolders(i).name, filesep, AnalysisresultsfolderPathName, filesep, sourceWellFolders(i).name, '_resultVAST.csv']);
                
                if ((resultFile_precheck2 || ~isempty(resultFile_precheck)) && Config.overwriteResults == 0)
                    disp('Segmentation already exists, skipping...');
                    continue;
                else
                    Write2FileVAST = 0;
                    
                    % Execute Python-based segmentation
                    [pythonSuccess, pythonOutput] = runPythonSegmentation(Config.VAST.pythonFile, VASTFile_new, resultFolder, Config);
                    if ~pythonSuccess
                        warning('Python segmentation failed for VAST: %s', pythonOutput);
                    end
                end
                
                % Load Python segmentation results
                resultFile = dir(strcat(saveFolder, AnalysisresultsfolderPathName, filesep, vast_Img_name, VASToutputDetectsufix));
                if (isempty(resultFile))
                    bw_python = [];
                else
                    resultFile = resultFile(1);
                    bw_python = imread(fullfile(resultFile.folder, resultFile.name));
                    bw_python = bw_python > 0;
                end
                
                % Load and analyze the VAST image
                I = imread(VASTFile_new);
                [~, filename, fext] = fileparts(VASTFile_new);
                
                % Determine view type (dorsal vs lateral)
                match1 = contains(strcat(filename, fext), '_1_1.tiff', 'IgnoreCase', true);
                match2 = contains(strcat(filename, fext), '_1_7.tiff', 'IgnoreCase', true);
                
                if (match1 == 1 || match2 == 1)
                    Operation_ = 0; % Dorsal view
                else
                    Operation_ = 1; % Lateral view
                end
                
                %% 1. BODY SEGMENTATION FOR ALL VIEWS
                % Resize image for deep learning segmentation
                I_re = imresize(I, Config.VAST.ImageSize);
                
                % Perform semantic segmentation using appropriate network
                if (Operation_ == 0) % Dorsal views
                    pxdsResults{1} = semanticseg(I_re, Config.VAST.net{1}); % Body segmentation
                    pxdsResults{2} = semanticseg(I_re, Config.VAST.net{2}); % Swim bladder segmentation
                else % Lateral views
                    pxdsResults{1} = semanticseg(I_re, Config.VAST.net{1}); % Body segmentation
                    pxdsResults{2} = semanticseg(I_re, Config.VAST.net{3}); % Eye, ear, backline, pericardium
                end
                
                % Extract body segmentation mask
                if (isempty(bw_python))
                    bwbodySegm = pxdsResults{1} == Config.VAST.classNames{1}(1); % Use MATLAB segmentation
                else
                    bwbodySegm = bw_python; % Use Python segmentation
                end
                
                % Resize body mask to original dimensions
                [nrow, ncol, ~] = size(I);
                bw_body = logical(imresize(bwbodySegm, [nrow, ncol], 'nearest'));
                
                %% Body Analysis (Common for all views)
                if (sum(bw_body(:)) ~= 0)
                    % Keep largest connected component
                    bw_body = logical(keepLargestn(bw_body));
                    
                    % Calculate MaxFeret diameter (longest dimension, less sensitive to shape variation)
                    MaxFeret = regionprops(bw_body > 0, 'MaxFeretProperties','Centroid');
                    
                    % Border detection: check if fish touches image border (>3 pixels)
                    borderThreshold = Config.VAST.BorderDetectionThreshold;
                    isTouchingBorder = numel(find(bw_body(:, [1:borderThreshold end-borderThreshold+1:end]) > 0)) > 0;
                    isTouchingBorderSize = sum(max(bw_body(:, [1:borderThreshold end-borderThreshold+1:end]), [], 2) > 0);
                    
                    % Store basic body measurements
                    resultData.(['Length_um', Config.VAST.useVAST_views{VASTview}]) = MaxFeret.MaxFeretDiameter * umConv;
                    resultData.(['Area_um', Config.VAST.useVAST_views{VASTview}]) = sum(sum(bw_body > 0)) * umConv^2;
                    resultData.(['IncompleteFish', Config.VAST.useVAST_views{VASTview}]) = num2str(isTouchingBorder);
                    resultData.(['IncompleteFish_EdgeSize_pix', Config.VAST.useVAST_views{VASTview}]) = num2str(isTouchingBorderSize);
                    
                    %% Tail Position Analysis (Common for all views)
                    % Calculate column-wise summation of body pixels (thickness profile)
                    fish_row_check_body_weight = sum(bw_body, 1);
                    fish_row_check_body = find(fish_row_check_body_weight > 0);
                    
                    % Adaptive region size: 10 pixels for long fish, half length for short fish
                    if (length(fish_row_check_body) > Config.VAST.TailAnalysisFixedRegion)
                        tale_check_region = Config.VAST.TailAnalysisFixedRegion;
                    else
                        tale_check_region = round(length(fish_row_check_body) * Config.VAST.TailAnalysisAdaptiveFactor);
                    end
                    
                    % Compare pixel density at both ends to determine tail position
                    if sum(fish_row_check_body_weight(fish_row_check_body(1:tale_check_region))) > sum(fish_row_check_body_weight(fish_row_check_body(end-tale_check_region+1:end)))
                        tale_position = fish_row_check_body(end); % Tail at end
                    else
                        tale_position = fish_row_check_body(1); % Tail at beginning
                    end
                    
                    %% 2. DORSAL VIEW ANALYSIS (Operation_ == 0)
                    if (Operation_ == 0)
                        % Initialize contour image for visualization
                        ContourImg = zeros(nrow, ncol);
                        ContourImg(bw_body) = 255;
                        ContourImg = edge(ContourImg);
                        
                        %% Swim Bladder Detection and Analysis
                        bw_swimbladder = logical(imresize(pxdsResults{2} == Config.VAST.classNames{2}(1), [nrow, ncol], 'nearest'));
                        ContourImg(bw_swimbladder) = 2;
                        
                        if (sum(bw_swimbladder(:)) ~= 0)
                            bw_swimbladder = logical(keepLargestn(bw_swimbladder));
                            fish_row_check_sb = find(sum(bw_swimbladder, 1) > 0);
                            
                            %% Pigment Region Analysis
                            % Define pigment region between tail and swim bladder
                            % Determine which end of swim bladder forms correct anatomical relationship
                            if (abs(fish_row_check_sb(1) - tale_position) > abs(fish_row_check_sb(end) - tale_position))
                                if (fish_row_check_sb(end) > tale_position)
                                    pigment_region = [tale_position fish_row_check_sb(end)];
                                else
                                    pigment_region = [fish_row_check_sb(end) tale_position];
                                end
                            else
                                if (fish_row_check_sb(1) > tale_position)
                                    pigment_region = [tale_position fish_row_check_sb(1)];
                                else
                                    pigment_region = [fish_row_check_sb(1) tale_position];
                                end
                            end
                            
                            % Create pigment region mask
                            bw_pigm = bw_body;
                            bw_pigm(:, 1:pigment_region(1)) = 0;
                            bw_pigm(:, pigment_region(end):end) = 0;
                            
                            % Convert to grayscale for intensity analysis
                            if (size(I, 3) == 3)
                                I_g = rgb2gray(I);
                                I_g_thresh = rgb2gray(I);
                            else
                                I_g = I;
                                I_g_thresh = I;
                            end
                            
                            % Calculate intensity statistics
                            Inetisy_min = min(I_g_thresh(:));
                            Inetisy_max = max(I_g_thresh(:));
                            
                            % Apply pigment region mask
                            I_g(~bw_pigm) = 0;
                            I_pigm_list = I_g(bw_pigm);
                            
                            % Thresholding for enhanced analysis
                            if (strcmp(class(I_g), 'uint8'))
                                I_pigm = 255 - I_g_thresh;
                            elseif (strcmp(class(I_g), 'uint16'))
                                I_pigm = 65535 - I_g_thresh;
                            else
                                I_pigm = 1 - I_g_thresh;
                            end
                            
                            I_pigm_inv_list = I_pigm(bw_pigm);
                            min_value_I_pigm = min(I_pigm_inv_list);
                            
                            try
                                I_pigm(~bw_pigm) = min_value_I_pigm;
                            catch
                                I_pigm(~bw_pigm) = 0;
                            end
                            
                            % Apply 70% threshold for enhanced pigment analysis
                            I_pigm_thresh = imbinarize(ImageNorm(double(I_pigm)), Config.VAST.PigmentThreshold);
                            I_g_thresh(~I_pigm_thresh) = 0;
                            
                            % Store pigment region measurements
                            resultData.(['Img_MinIntensity', Config.VAST.useVAST_views{VASTview}]) = Inetisy_min;
                            resultData.(['Img_MaxIntensity', Config.VAST.useVAST_views{VASTview}]) = Inetisy_max;
                            resultData.(['Pigm_TotalIntensity', Config.VAST.useVAST_views{VASTview}]) = sum(sum(double(I_g)));
                            resultData.(['Pigm_Area_pixel', Config.VAST.useVAST_views{VASTview}]) = sum(sum(double(bw_pigm)));
                            resultData.(['Pigm_Area_um', Config.VAST.useVAST_views{VASTview}]) = sum(sum(double(bw_pigm))) * umConv^2;
                            resultData.(['Pigm_MeanIntensityPerPixel', Config.VAST.useVAST_views{VASTview}]) = mean(I_pigm_list(:));
                            resultData.(['Pigm_Thresh70_TotalIntensity', Config.VAST.useVAST_views{VASTview}]) = sum(sum(double(I_g_thresh)));
                            resultData.(['Pigm_Thresh70_Area_um', Config.VAST.useVAST_views{VASTview}]) = sum(sum(double(I_pigm_thresh))) * umConv^2;
                            resultData.(['Pigm_Thresh70_MeanIntensity', Config.VAST.useVAST_views{VASTview}]) = mean(I_g_thresh(:));
                        else
                            % No swim bladder detected
                            I_pigm_thresh = zeros(size(bw_swimbladder));
                            resultData.(['Img_MinIntensity', Config.VAST.useVAST_views{VASTview}]) = 0;
                            resultData.(['Img_MaxIntensity', Config.VAST.useVAST_views{VASTview}]) = 0;
                            resultData.(['Pigm_TotalIntensity', Config.VAST.useVAST_views{VASTview}]) = 0;
                            resultData.(['Pigm_Area_pixel', Config.VAST.useVAST_views{VASTview}]) = 0;
                            resultData.(['Pigm_Area_um', Config.VAST.useVAST_views{VASTview}]) = 0;
                            resultData.(['Pigm_MeanIntensityPerPixel', Config.VAST.useVAST_views{VASTview}]) = 0;
                            resultData.(['Pigm_Thresh70_TotalIntensity', Config.VAST.useVAST_views{VASTview}]) = 0;
                            resultData.(['Pigm_Thresh70_Area_um', Config.VAST.useVAST_views{VASTview}]) = 0;
                            resultData.(['Pigm_Thresh70_MeanIntensity', Config.VAST.useVAST_views{VASTview}]) = 0;
                        end
                        
                        % Create visualization image
                        if (size(I, 3) == 1)
                            Original_Img = repmat(I, [1 1 3]);
                        else
                            Original_Img = I;
                        end
                        IRGB = labeloverlay(labeloverlay(Original_Img, ContourImg), edge(I_pigm_thresh), 'Colormap', 'hot');
                        
                    %% 3. LATERAL VIEW ANALYSIS (Operation_ == 1)
                    elseif (Operation_ == 1)
                        % Initialize contour image for visualization
                        ContourImg = zeros(nrow, ncol);
                        ContourImg(bw_body) = 255;
                        ContourImg = edge(ContourImg);
                        
                        %% Eye Detection and Analysis
                        eyeMask = logical(imresize(pxdsResults{2} == Config.VAST.classNames{3}(1), [nrow ncol], 'nearest'));
                        if (sum(eyeMask(:)) ~= 0)
                            eyeMask = logical(keepLargestn(eyeMask));
                        end
                        % Morphological refinement with disk-shaped structuring element (radius 3)
                        eyeMask = imclose(eyeMask, strel('disk', Config.VAST.EyeMorphDiskRadius));
                        eyeMask = logical(RefineShape(eyeMask, I));
                        ContourImg(eyeMask) = Config.VAST.labelIDs{3}(1);
                        
                        % Calculate eye measurements
                        Eyeprop = regionprops(eyeMask, 'Centroid', 'MaxFeretProperties', 'MinFeretProperties', 'EquivDiameter');
                        
                        if (~isempty(Eyeprop))
                            resultData.(['Eye_maxDiameter_MaxFeret', Config.VAST.useVAST_views{VASTview}]) = Eyeprop(1).MaxFeretDiameter * umConv;
                            resultData.(['Eye_meanDiameter_EquivCircle', Config.VAST.useVAST_views{VASTview}]) = Eyeprop(1).EquivDiameter * umConv;
                            resultData.(['Eye_minDiameter_MaxFeret', Config.VAST.useVAST_views{VASTview}]) = Eyeprop(1).MinFeretDiameter * umConv;
                            resultData.(['Eye_Area', Config.VAST.useVAST_views{VASTview}]) = sum(sum(double(eyeMask))) * umConv^2;
                        else
                            resultData.(['Eye_maxDiameter_MaxFeret', Config.VAST.useVAST_views{VASTview}]) = 0;
                            resultData.(['Eye_meanDiameter_EquivCircle', Config.VAST.useVAST_views{VASTview}]) = 0;
                            resultData.(['Eye_minDiameter_MaxFeret', Config.VAST.useVAST_views{VASTview}]) = 0;
                            resultData.(['Eye_Area', Config.VAST.useVAST_views{VASTview}]) = 0;
                        end
                        
                        %% Ear Detection and Analysis
                        EarMask = logical(imresize(pxdsResults{2} == Config.VAST.classNames{3}(2), [nrow ncol], 'nearest'));
                        if (sum(EarMask(:)) ~= 0)
                            EarMask = logical(keepLargestn(EarMask));
                        end
                        ContourImg(EarMask) = Config.VAST.labelIDs{3}(2);
                        
                        % Record ear centroid for developmental age assessment
                        Earprop = regionprops(EarMask, 'Centroid');
                        if (~isempty(Earprop))
                            XY_Ear = Earprop(1).Centroid;
                        end
                        
                        %% Backline Analysis
                        BacklineMask = logical(imresize(pxdsResults{2} == Config.VAST.classNames{3}(3), [nrow ncol], 'nearest'));
                        if (sum(BacklineMask(:)) > 1)
                            BacklineMask = logical(keepLargestn(BacklineMask));
                        end
                        
                        if (sum(BacklineMask(:)) > 1)
                            % Morphological close operation with disk-shaped structuring element (radius 15)
                            BacklineMask = imclose(BacklineMask, strel('disk', Config.VAST.BacklineMorphDiskRadius));
                            % Thin backline to single-pixel width
                            BacklineMask_thin = bwmorph(BacklineMask, 'thin', Config.VAST.BacklineThinIterations);
                            BacklineMask_Endp = bwmorph(BacklineMask_thin, 'endpoints');
                            ContourImg(BacklineMask) = Config.VAST.labelIDs{3}(3);
                            
                            % Calculate Feret diameters
                            BacklineMaskprop = regionprops(BacklineMask, 'Centroid', 'MaxFeretProperties', 'MinFeretProperties');
                            BacklineMaskThinprop = regionprops(BacklineMask_thin, 'PixelList');
                            BacklineMaskEndpprop = regionprops(BacklineMask_Endp, 'PixelList');
                            
                            XY_backline = BacklineMaskThinprop(1).PixelList;
                            XY_lineEnds(1, :) = BacklineMaskEndpprop(1).PixelList(1, :);
                            try
                                XY_lineEnds(2, :) = BacklineMaskEndpprop(1).PixelList(2, :);
                            catch
                                XY_lineEnds(2, :) = BacklineMaskEndpprop(2).PixelList(1, :);
                            end
                            
                            % Add tail position points to backline for curvature analysis
                            xcorr_tail = find(bw_body(:, tale_position) == 1);
                            if (~isempty(xcorr_tail))
                                XY_append = [repmat(tale_position, [length(xcorr_tail) 1]) xcorr_tail];
                                XY_backline = [XY_append; XY_backline];
                            end
                            
                            % Curvature analysis using circle fitting algorithm (Pratt method)
                            Par = CircleFitByPratt(XY_backline);
                            
                            resultData.(['b_MaxFeret', Config.VAST.useVAST_views{VASTview}]) = BacklineMaskprop(1).MaxFeretDiameter * umConv;
                            resultData.(['a_MinFeret', Config.VAST.useVAST_views{VASTview}]) = BacklineMaskprop(1).MinFeretDiameter * umConv;
                            resultData.(['Curvature_CircleFit', Config.VAST.useVAST_views{VASTview}]) = Par(3) * umConv;
                        else
                            XY_backline = [];
                            resultData.(['b_MaxFeret', Config.VAST.useVAST_views{VASTview}]) = 0;
                            resultData.(['a_MinFeret', Config.VAST.useVAST_views{VASTview}]) = 0;
                            resultData.(['Curvature_CircleFit', Config.VAST.useVAST_views{VASTview}]) = 0;
                        end
                        
                        %% Pericardium Detection and Analysis
                        pericardMask = logical(imresize(pxdsResults{2} == Config.VAST.classNames{3}(4), [nrow ncol], 'nearest'));
                        if (sum(pericardMask(:)) ~= 0)
                            pericardMask = logical(keepLargestn(pericardMask));
                        end
                        ContourImg(pericardMask) = Config.VAST.labelIDs{3}(4);
                        
                        resultData.(['Pericard_Area', Config.VAST.useVAST_views{VASTview}]) = sum(sum(double(pericardMask))) * umConv^2;
                        
                        %% Developmental Age Assessment
                        if (~isempty(Eyeprop) && ~isempty(Earprop) && ~isempty(XY_backline))
                            % Calculate developmental angle between reference lines
                            d1 = norm(XY_lineEnds(1, :) - XY_Ear);
                            d2 = norm(XY_lineEnds(2, :) - XY_Ear);
                            
                            % Determine which endpoint is closer to ear
                            if (d1 > d2)
                                XY1 = XY_lineEnds(2, :);
                                XY2 = XY_lineEnds(1, :);
                            else
                                XY1 = XY_lineEnds(1, :);
                                XY2 = XY_lineEnds(2, :);
                            end
                            
                            XY3 = XY_Ear;
                            XY4 = Eyeprop(1).Centroid;
                            
                            % Calculate angle between backline and ear-eye line
                            theta1 = atan((XY2(1) - XY1(1)) / (XY2(2) - XY1(2)));
                            theta2 = atan((XY4(1) - XY3(1)) / (XY4(2) - XY3(2)));
                            theta = abs((theta1 - theta2) * Config.VAST.DevelopmentalAngleConversion);
                            
                            % Create visualization with reference lines
                            RGB = insertShape(I, 'Line', [XY1 XY2; XY3 XY4], 'LineWidth', 3, 'Color', 'green');
                            RGB = insertShape(RGB, 'Circle', [Par], 'LineWidth', 3, 'Color', 'red');
                            
                            % Determine bending direction using cross product analysis
                            d_d1 = (MaxFeret(1).Centroid(2) - XY1(2)) * (XY2(1) - XY1(1)) - (MaxFeret(1).Centroid(1) - XY1(1)) * (XY2(2) - XY1(2));
                            d_d2 = (Par(2) - XY1(2)) * (XY2(1) - XY1(1)) - (Par(1) - XY1(1)) * (XY2(2) - XY1(2));
                            
                            if (d_d1 * d_d2 > 0)
                                Bent_in = 1; % Bent inward
                            else
                                Bent_in = 0; % Bent outward
                            end
                            
                            resultData.(['Angle_BtwLines', Config.VAST.useVAST_views{VASTview}]) = theta;
                            resultData.(['Bent_direction', Config.VAST.useVAST_views{VASTview}]) = Bent_in;
                        else
                            % Missing required structures for developmental assessment
                            if (size(I, 3) == 1)
                                RGB = repmat(I, [1 1 3]);
                            else
                                RGB = I;
                            end
                            resultData.(['Angle_BtwLines', Config.VAST.useVAST_views{VASTview}]) = 0;
                            resultData.(['Bent_direction', Config.VAST.useVAST_views{VASTview}]) = 0;
                        end
                        
                        IRGB = labeloverlay(RGB, ContourImg);
                    end
                    
                else
                    % No body detected - set all measurements to zero
                    resultData.(['Length_um', Config.VAST.useVAST_views{VASTview}]) = 0;
                    resultData.(['Area_um', Config.VAST.useVAST_views{VASTview}]) = 0;
                    resultData.(['IncompleteFish', Config.VAST.useVAST_views{VASTview}]) = 0;
                    resultData.(['IncompleteFish_EdgeSize_pix', Config.VAST.useVAST_views{VASTview}]) = 0;
                    
                    if (Operation_ == 0) % Dorsal view
                        resultData.(['Img_MinIntensity', Config.VAST.useVAST_views{VASTview}]) = 0;
                        resultData.(['Img_MaxIntensity', Config.VAST.useVAST_views{VASTview}]) = 0;
                        resultData.(['Pigm_TotalIntensity', Config.VAST.useVAST_views{VASTview}]) = 0;
                        resultData.(['Pigm_Area_pixel', Config.VAST.useVAST_views{VASTview}]) = 0;
                        resultData.(['Pigm_Area_um', Config.VAST.useVAST_views{VASTview}]) = 0;
                        resultData.(['Pigm_MeanIntensityPerPixel', Config.VAST.useVAST_views{VASTview}]) = 0;
                        resultData.(['Pigm_Thresh70_TotalIntensity', Config.VAST.useVAST_views{VASTview}]) = 0;
                        resultData.(['Pigm_Thresh70_Area_um', Config.VAST.useVAST_views{VASTview}]) = 0;
                        resultData.(['Pigm_Thresh70_MeanIntensity', Config.VAST.useVAST_views{VASTview}]) = 0;
                    else % Lateral view
                        resultData.(['Eye_maxDiameter_MaxFeret', Config.VAST.useVAST_views{VASTview}]) = 0;
                        resultData.(['Eye_meanDiameter_EquivCircle', Config.VAST.useVAST_views{VASTview}]) = 0;
                        resultData.(['Eye_minDiameter_MaxFeret', Config.VAST.useVAST_views{VASTview}]) = 0;
                        resultData.(['Eye_Area', Config.VAST.useVAST_views{VASTview}]) = 0;
                        resultData.(['b_MaxFeret', Config.VAST.useVAST_views{VASTview}]) = 0;
                        resultData.(['a_MinFeret', Config.VAST.useVAST_views{VASTview}]) = 0;
                        resultData.(['Curvature_CircleFit', Config.VAST.useVAST_views{VASTview}]) = 0;
                        resultData.(['Pericard_Area', Config.VAST.useVAST_views{VASTview}]) = 0;
                        resultData.(['Angle_BtwLines', Config.VAST.useVAST_views{VASTview}]) = 0;
                        resultData.(['Bent_direction', Config.VAST.useVAST_views{VASTview}]) = 0;
                    end
                    
                    IRGB = I;
                end
                
                % Save visualization image
                ImageFileFolder = [saveRootFolder, filesep, sourceWellFolders(i).name, filesep, AnalysisresultsfolderPathName, filesep];
                ImageFilename = [ImageFileFolder, sourceWellFolders(i).name, Config.VAST.useVAST_views{VASTview}, 'res.png'];
                
                if ~exist((ImageFileFolder), 'dir')
                    mkdir(ImageFileFolder);
                end
                imwrite(uint8(IRGB), ImageFilename);
            end
        end
        
        % Save individual well results
        if (Write2FileVAST == 0)
            resultData.Date = welldate;
            resultFilename = [saveRootFolder, filesep, sourceWellFolders(i).name, filesep, AnalysisresultsfolderPathName, filesep, sourceWellFolders(i).name, '_resultVAST.csv'];
            
            % Ensure the output folder exists
            resultFolder = [saveRootFolder, filesep, sourceWellFolders(i).name, filesep, AnalysisresultsfolderPathName, filesep];
            if (~isfolder(resultFolder))
                mkdir(resultFolder);
            end
            
            writetable(struct2table(resultData, 'AsArray', true), resultFilename, 'Delimiter', ';')
        end
    end
end

%% Collect and aggregate all results for the experiment
resultData = struct();
inc = 1;
for i = 1:numel(sourceWellFolders)
    resultData.ExperimentName = strcat(experimentName, '_', sourceWellFolders(i).name);
    
    if (Config.ReadfromNAS == 1)
        String_segments = split(sourceWellFolders(i).folder, filesep);
        geneName = String_segments{end-1};
        experimentName = String_segments{end};
        wellname = sourceWellFolders(i).name;
        resultData.Wellfolder = wellname;
        resultFileName = strcat(Config.localNASrootPath, filesep, geneName, filesep, experimentName, filesep, wellname, filesep, AnalysisresultsfolderPathName, filesep, sourceWellFolders(i).name, AnalysisVASTResultsSuffix);
    else
        resultData.Wellfolder = sourceWellFolders(i).name;
        resultFileName = [sourceWellFolders(i).folder, filesep, sourceWellFolders(i).name, filesep, AnalysisresultsfolderPathName, filesep, sourceWellFolders(i).name, AnalysisVASTResultsSuffix];
    end
    
    if (isfile(resultFileName))
        importresult{inc} = readtable(resultFileName);
        inc = inc + 1;
    end
end

% Combine all results with error handling for NaN values
try
    if (numel(sourceWellFolders) ~= 0)
        allTables = vertcat(importresult{:});
    else
        return;
    end
catch
    disp('ERROR: check NAN values in the results!')
    for i = 1:length(importresult)
        check_table = importresult{i};
        if (~isdatetime(check_table.Date)) % logical index
            check_table.Date
            check_table.Wellfolder
        end
    end
    for i = 1:length(importresult)
        if (~isdatetime(importresult{i}.Date))
            importresult{i}.Date
            importresult{i}.Wellfolder
            importresult{i}.Date = NaT;
        end
    end
    allTables = vertcat(importresult{:});
end

% Save aggregated results
filename = strcat(saveRootFolder, filesep, SaveAnalysisResultName);
writetable(allTables, filename, 'Sheet', 'Vast', 'WriteMode', 'overwritesheet');
disp("VAST image analysis completed");

end