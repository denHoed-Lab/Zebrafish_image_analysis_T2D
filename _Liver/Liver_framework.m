function logFile = Liver_framework(rootfolder, Config)
% LIVER_FRAMEWORK - Main framework for liver tissue and lipid analysis
%
% This function orchestrates the complete liver analysis pipeline:
% 1. Scans for liver data files in the analysis folder
% 2. Processes each file using deep learning segmentation
% 3. Collects and aggregates results across all wells
% 4. Saves comprehensive analysis results
%
% The analysis provides four key quantitative measurements:
% - Liver area (segmented and refined)
% - Lipid area (within liver region)
% - Number of lipid droplets
% - Total lipid intensity
%
% Inputs:
%   rootfolder - Path to the folder containing experiment data
%   Config - Configuration structure with all analysis parameters
%
% Outputs:
%   logFile - Cell array containing analysis log messages
%
% Author: Hanqing Zhang, hanzha@kth.se
% Version: 1.1, 20250619

%% Initialize variables
valid_digit = 1; % Precision for rounding results

% Get all subdirectories (wells) in the analysis folder
folders = dir(rootfolder);
folders = folders([folders(:).isdir]); % Keep only directories
folders = folders(~ismember({folders(:).name}, {'.', '..'})); % Remove . and ..

analyzed_liverdata = 0; % Counter for processed files

%% Setup logging system
LogfileName = 'LiverLog.txt';
Logfullpath = strcat(rootfolder, filesep, LogfileName);

% Create log file if it doesn't exist
if (~isfile(Logfullpath))
    mkdir(rootfolder);
    writematrix([], Logfullpath, 'WriteMode', 'overwrite');
end

logFile = []; % Initialize log array
loadprivFile = 0; % Flag for re-analyzing previous files

%% Process each well (subdirectory)
for i = 1:numel(folders)
    folderpath = [folders(i).folder, filesep, folders(i).name, filesep];
    
    % Get all files in the well folder
    d = dir(folderpath);
    
    % Find files matching liver naming patterns
    match = contains({d.name}, Config.Liver.Expr1, 'IgnoreCase', true); % First expression (e.g., 'liver')
    files = d(match); % Files matching the pattern
    
    % Log current well being processed
    stringTmp = ['Analyzing Well ' folders(i).name];
    logFile = addtolog(logFile, stringTmp, Logfullpath);
    
    % Check if well has liver data
    if (numel(files) == 0)
        [~, foldername, ~] = fileparts(folderpath(1:end-1));
        stringTmp = [foldername ' has no Liver data.'];
        logFile = addtolog(logFile, stringTmp, Logfullpath);
        continue
    end
    
    %% Process each liver file in the well
    for fId = 1:numel(files)
        % Check if file matches second expression (e.g., 'lng')
        isThunder = contains(files(fId).name, Config.Liver.Expr2, 'IgnoreCase', true);
        
        % Determine which files to process based on configuration
        if Config.Liver.useThunder
            if (isThunder)
                file2read = [folderpath, files(fId).name];
            else
                continue;
            end
        else
            if (isThunder)
                continue
            else
                file2read = [folderpath, files(fId).name];
            end
        end
        
        % Extract file information for output naming
        [~, name_nonext, ~] = fileparts(files(fId).name);
        [~, name_folder, ~] = fileparts(files(fId).folder);
        
        %% Determine result file path
        checkResultFile = strcat(folderpath, Config.Liver.outputFolderlist{1}, filesep, name_nonext, Config.saveResultSuffix, '.xlsx');
        
        %% Check if file has already been analyzed (skip if overwrite is disabled)
        if (Config.overwriteResults ~= 1)
            if isfile(checkResultFile)
                analyzed_liverdata = analyzed_liverdata + 1;
                previousFileName = file2read;
                
                [~, logfilename, ~] = fileparts(checkResultFile);
                stringTmp = [logfilename ' already exists'];
                logFile = addtolog(logFile, stringTmp, Logfullpath);
                
                % Allow rewrite of previous data if needed
                loadprivFile = 1;
                continue
            end
        end
        
        %% Analyze liver files using specified operation mode
        [logpath, logfilename, ~] = fileparts(file2read);
        [~, wellfilename, ~] = fileparts(logpath);
        stringTmp = ['Analyzing ' file2read];
        logFile = addtolog(logFile, stringTmp, Logfullpath);
        
        % Execute analysis based on operation mode
        switch Config.Liver.Operation
            case 'DL_3ch'
                Iliver = imreadMPTiff(file2read, 1);
                fileResultsName = strcat(rootfolder, filesep, wellfilename, filesep, logfilename);
                [~, ~] = Liver_process(fileResultsName, Iliver, Config);
            otherwise
                error('No valid operation')
        end
        
        analyzed_liverdata = analyzed_liverdata + 1;
        
        %% Re-analyze previous file if needed (for incomplete analysis)
        if (loadprivFile == 1)
            if (Config.skipOverwrite ~= 1)
                Config.overwriteResults = 1;
            else
                continue
            end
            
            [~, logfilename, ~] = fileparts(previousFileName);
            stringTmp = ['Analyzing ' logfilename];
            logFile = addtolog(logFile, stringTmp, Logfullpath);
            
            % Re-execute analysis on previous file
            switch Config.Liver.Operation
                case 'DL_3ch'
                    Iliver = imreadMPTiff(file2read, 1);
                    fileResultsName = strcat(rootfolder, filesep, wellfilename, filesep, logfilename);
                    [~, ~] = Liver_process(fileResultsName, Iliver, Config);
                otherwise
                    error('No valid operation')
            end
            
            loadprivFile = 0;
        end
    end
end

%% Log analysis completion
stringTmp = ['Analyzed ' num2str(analyzed_liverdata) ' files.'];
logFile = addtolog(logFile, stringTmp, Logfullpath);

stringTmp = 'Liver DL data collection completed.';
logFile = addtolog(logFile, stringTmp, Logfullpath);

%% Collect and aggregate results from all wells
usedfileCount = 0;
usedfileNames = {};

% Initialize structure to store liver information
LiverInfo = struct('WellName', [], 'FileName', [], 'LocalName', [], 'ModTime', [], 'Statistics', []);
LiverCount = 1;

% Process each well for result collection
for i = 1:numel(folders)
    folderpath = [folders(i).folder, filesep, folders(i).name, filesep];
    
    % Create local name for the well
    [~, folder_name, ~] = fileparts(folders(i).folder);
    local_name = strcat(folder_name, '_', folders(i).name);
    
    % Get VAST image timestamp for metadata
    vastfolder = strcat(folderpath, 'VAST', filesep);
    try
        vast_imds = imageDatastore(vastfolder, 'IncludeSubfolders', 0, 'FileExtensions', {'.jpg', '.tif', '.bmp', '.png', '.tiff'});
        FileNames = vast_imds.Files;
        fileImginfo = imfinfo(FileNames{1});
        ImageDateTime = fileImginfo.FileModDate;
    catch
        disp(['LiverAnalysis: Cannot find VAST folder, ' vastfolder]);
        ImageDateTime = 'NaT';
    end
    
    % Find liver result files
    filesliver = dir([folderpath, '*Liver*']);
    
    for fId = 1:numel(filesliver)
        if (filesliver(fId).isdir == 1)
            continue; % Skip directories
        end
        
        [~, name_nonext, ~] = fileparts(filesliver(fId).name);
        [~, well_name, ~] = fileparts(filesliver(fId).folder);
        
        usedfileCount = usedfileCount + 1;
        usedfileNames{end + 1} = name_nonext;
        
        % Determine result file path
       checkResultFile = strcat(folderpath, Config.Liver.outputFolderlist{1}, filesep, name_nonext, Config.saveResultSuffix, '.xlsx');
        
        % Load result data if file exists
        if isfile(checkResultFile)
            LiverInfo(LiverCount).WellName = folders(i).name;
            LiverInfo(LiverCount).FileName = name_nonext;
            LiverInfo(LiverCount).LocalName = local_name;
            LiverInfo(LiverCount).ModTime = ImageDateTime;
            LiverInfo(LiverCount).Statistics = readtable(checkResultFile);
            LiverCount = LiverCount + 1;
        end
    end
end

%% Aggregate and save comprehensive results
if (~isempty(LiverInfo(1).WellName))
    % Initialize structure for saving
    structSave = struct();
    
    % Get field names for data organization
    MainFieldNames = fieldnames(LiverInfo);
    init_i = numel(MainFieldNames) - 1; % Number of main fields (excluding Statistics)
    subFieldNames = LiverInfo(1).Statistics.Properties.VariableNames; % Statistical field names
    Allfields = {MainFieldNames{1:end-1} subFieldNames{:}}; % Combined field names
    
    % Reorganize data structure for table format
    for i = 1:numel(LiverInfo)
        % Copy main fields (WellName, FileName, etc.)
        for j = 1:init_i
            structSave(i).(Allfields{j}) = LiverInfo(i).(Allfields{j});
        end
        
        % Copy statistical fields
        for j = 1:numel(subFieldNames)
            structSave(i).(Allfields{init_i + j}) = LiverInfo(i).Statistics{1, j};
        end
    end
    
    % Convert to table format
    TableSave = struct2table(structSave);
    
    %% Save results to Excel file
    filename = strcat(rootfolder, filesep, 'Liver_results.xlsx');
    writetable(TableSave, filename, 'Sheet', 'Liver', 'WriteMode', 'overwritesheet');
    disp("Liver image analysis completed");
end

end

%% Helper function for logging
function logFile = addtolog(logFile, stringTmp, fullpath)
% ADDTOLOG - Add message to log file and display it
%
% Inputs:
%   logFile - Current log cell array
%   stringTmp - Message to add
%   fullpath - Path to log file
%
% Outputs:
%   logFile - Updated log cell array

logFile{end + 1} = stringTmp;
savestring = ['Time: ' datestr(datetime) ' , Op: ' logFile{end}];

try
    writematrix(convertCharsToStrings(savestring), fullpath, 'WriteMode', 'append');
catch
    error('logFile: cannot write to log file')
end

disp(savestring);
end
