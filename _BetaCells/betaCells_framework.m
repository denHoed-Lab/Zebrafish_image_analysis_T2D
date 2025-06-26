function logFile = betaCells_framework(analysisfolder, Config)
% betaCells_framework - Main framework for beta cell analysis using deep learning
% The framework is to adapt networks for different experiemntal conventions
%
% This function orchestrates the complete beta cell analysis pipeline:
% 1. Scans for beta cell data files in the analysis folder
% 2. Processes each file using the specified operation mode
% 3. Collects and aggregates results across all wells
% 4. Saves comprehensive analysis results
%
% Inputs:
%   analysisfolder - Path to the folder containing experiment data
%   Config - Configuration structure with all analysis parameters
%
% Outputs:
%   logFile - Cell array containing analysis log messages
%
% Author: Hanqing Zhang, hanzha@kth.se
% Version: 1.1, 20250619

%% Initialize variables and setup
valid_digit = 1; % Precision for rounding results

% Get all subdirectories (wells) in the analysis folder
folders = dir(analysisfolder);
folders = folders([folders(:).isdir]); % Keep only directories
folders = folders(~ismember({folders(:).name}, {'.', '..'})); % Remove . and ..

analyzed_betacells = 0; % Counter for processed files

%% Setup logging system
LogfileName = 'BetacellLog.txt';
Logfullpath = strcat(analysisfolder, filesep, LogfileName);

% Create log file if it doesn't exist
if (~isfile(Logfullpath))
    mkdir(analysisfolder);
    writematrix([], Logfullpath, 'WriteMode', 'overwrite');
end

logFile = []; % Initialize log array
loadprivFile = 0; % Flag for re-analyzing previous files

%% Process each well (subdirectory)
for i = 1:numel(folders)
    % Construct full path to current well folder
    folderpath = [folders(i).folder, filesep, folders(i).name, filesep];
    
    % Get all files in the well folder
    d = dir(folderpath);
    
    % Find files matching beta cell naming patterns
    match1 = contains({d.name}, Config.Betacell.BetaCellExpr1, 'IgnoreCase', true); % First expression (e.g., 'betacell')
    match2 = contains({d.name}, Config.Betacell.BetaCellExpr2, 'IgnoreCase', true); % Second expression (e.g., 'lng')
    files = d(and(match1, match2)); % Files matching both patterns
    
    % Log current well being processed
    stringTmp = ['Analyzing Well ' folders(i).name];
    logFile = addtolog(logFile, stringTmp, Logfullpath);
    
    % Check if well has beta cell data
    if (numel(files) == 0)
        [~, foldername, ~] = fileparts(folderpath(1:end-1));
        stringTmp = [foldername ' has no betacell data.'];
        logFile = addtolog(logFile, stringTmp, Logfullpath);
    end
    
    %% Process each beta cell file in the well
    for fId = 1:numel(files)
        % Construct full path to current file
        file2read = [folderpath, files(fId).name];
        [~, name_nonext, ~] = fileparts(files(fId).name);
        
        %% Determine result file path based on data source
        % Locally saved extraction - use local folder path
        checkResultFile = strcat(folderpath, Config.Betacell.outputFolderlist{1}, filesep, name_nonext, Config.saveResultSuffix, '.xlsx');

        %% Check if file has already been analyzed (skip if overwrite is disabled)
        if (Config.overwriteResults ~= 1)
            if isfile(checkResultFile)
                analyzed_betacells = analyzed_betacells + 1;
                previousFileName = file2read;
                
                [~, logfilename, ~] = fileparts(checkResultFile);
                stringTmp = [logfilename ' already exists'];
                logFile = addtolog(logFile, stringTmp, Logfullpath);
                
                % Allow rewrite of previous data if needed
                loadprivFile = 1;
                continue
            end
        end
        
        %% Analyze beta cells using specified operation mode
        [~, logfilename, ~] = fileparts(file2read);
        stringTmp = ['Analyzing ' logfilename];
        logFile = addtolog(logFile, stringTmp, Logfullpath);
        
        % Load image data
        Ibeta = imreadMPTiff(file2read, 1);
        
        % Execute analysis based on operation mode
        switch Config.Betacell.Operation
            case 'DL_Betacell_vast2'
                [~, ~] = betaCells_process(file2read, Ibeta, Config);
            otherwise
                error('No valid operation')
        end
        
        analyzed_betacells = analyzed_betacells + 1;
        
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
            switch Config.Betacell.Operation
                case 'DL_Betacell_vast2'
                    [~, ~] = betaCells_process(file2read, Ibeta, Config);
                otherwise
                    error('No valid operation')
            end

            analyzed_betacells = analyzed_betacells + 1;
        end
    end
end

%% Log analysis completion
stringTmp = ['Analyzed ' num2str(analyzed_betacells) ' files.'];
logFile = addtolog(logFile, stringTmp, Logfullpath);

stringTmp = 'Betacell analysis completed.';
logFile = addtolog(logFile, stringTmp, Logfullpath);

%% Collect and aggregate results from all wells
usedfileCount = 0;
usedfileNames = {};

% Initialize structure to store beta cell information
BetacellInfo = struct('WellName', [], 'FileName', [], 'LocalName', [], 'ModTime', [], 'Statistics', []);
BetacellCount = 1;

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
        disp(['betaCells_v4: Cannot find VAST folder, ' vastfolder]);
        ImageDateTime = 'NaT';
    end
    
    % Find beta cell result files
    filesBeta = dir([folderpath, '*Betacell*']);
    
    for fId = 1:numel(filesBeta)
        if (filesBeta(fId).isdir == 1)
            continue; % Skip directories
        end
        
        [~, name_nonext, ~] = fileparts(filesBeta(fId).name);
        
        usedfileCount = usedfileCount + 1;
        usedfileNames{end + 1} = name_nonext;
        
        % Determine result file path based on data source
        checkResultFile = strcat(folderpath, Config.Betacell.outputFolderlist{1}, filesep, name_nonext, Config.saveResultSuffix, '.xlsx');

        % Load result data if file exists
        if isfile(checkResultFile)
            BetacellInfo(BetacellCount).WellName = folders(i).name;
            BetacellInfo(BetacellCount).FileName = name_nonext;
            BetacellInfo(BetacellCount).LocalName = local_name;
            BetacellInfo(BetacellCount).ModTime = ImageDateTime;
            BetacellInfo(BetacellCount).Statistics = readtable(checkResultFile);
            BetacellCount = BetacellCount + 1;
        end
    end
end

%% Aggregate and save comprehensive results
if (~isempty(BetacellInfo(1).WellName))
    % Initialize structure for saving
    structSave = struct();
    
    % Get field names for data organization
    MainFieldNames = fieldnames(BetacellInfo);
    init_i = numel(MainFieldNames) - 1; % Number of main fields (excluding Statistics)
    subFieldNames = BetacellInfo(1).Statistics.Properties.VariableNames; % Statistical field names
    Allfields = {MainFieldNames{1:end-1} subFieldNames{:}}; % Combined field names
    
    % Reorganize data structure for table format
    for i = 1:numel(BetacellInfo)
        % Copy main fields (WellName, FileName, etc.)
        for j = 1:init_i
            structSave(i).(Allfields{j}) = BetacellInfo(i).(Allfields{j});
        end
        
        % Copy statistical fields with error handling
        for j = 1:numel(subFieldNames)
            try
                structSave(i).(Allfields{init_i + j}) = BetacellInfo(i).Statistics{1, j};
            catch
                structSave(i).(Allfields{init_i + j}) = 0; % Default value if data missing
            end
        end
    end
    % Convert to table format
    TableSave = struct2table(structSave);
    
    %% Convert pixel units to micrometers for meaningful measurements
    TableSaveMicro = TableSave;
    
    % Volume measurements (pixel^3 to μm^3)
    TableSaveMicro{:, init_i + 1} = round(valid_digit * TableSaveMicro{:, init_i + 1} .* Config.outputResolution^3) / valid_digit; % Total volume
    TableSaveMicro{:, init_i + 4} = round(valid_digit * TableSaveMicro{:, init_i + 4} .* Config.outputResolution^3) / valid_digit; % Cell volume 1
    TableSaveMicro{:, init_i + 5} = round(valid_digit * TableSaveMicro{:, init_i + 5} .* Config.outputResolution^3) / valid_digit; % Cell volume 2
    TableSaveMicro{:, init_i + 16} = round(valid_digit * TableSaveMicro{:, init_i + 16} .* Config.outputResolution^3) / valid_digit; % Volume all
    
    % Diameter measurements (pixels to μm)
    TableSaveMicro{:, init_i + 2} = round(valid_digit * TableSaveMicro{:, init_i + 2} .* Config.outputResolution) / valid_digit; % Equivalent diameter 1
    TableSaveMicro{:, init_i + 3} = round(valid_digit * TableSaveMicro{:, init_i + 3} .* Config.outputResolution) / valid_digit; % Equivalent diameter 2
    TableSaveMicro{:, init_i + 18} = round(valid_digit * TableSaveMicro{:, init_i + 18} .* Config.outputResolution) / valid_digit; % Equivalent diameter 3
    TableSaveMicro{:, init_i + 19} = round(valid_digit * TableSaveMicro{:, init_i + 19} .* Config.outputResolution) / valid_digit; % Axis length 1
    TableSaveMicro{:, init_i + 20} = round(valid_digit * TableSaveMicro{:, init_i + 20} .* Config.outputResolution) / valid_digit; % Axis length 2
    TableSaveMicro{:, init_i + 21} = round(valid_digit * TableSaveMicro{:, init_i + 21} .* Config.outputResolution) / valid_digit; % Axis length 3
    
    % Surface area measurements (pixel^2 to μm^2)
    TableSaveMicro{:, init_i + 6} = round(valid_digit * TableSaveMicro{:, init_i + 6} .* Config.outputResolution^2) / valid_digit; % Surface area 1
    TableSaveMicro{:, init_i + 7} = round(valid_digit * TableSaveMicro{:, init_i + 7} .* Config.outputResolution^2) / valid_digit; % Surface area 2
    TableSaveMicro{:, init_i + 17} = round(valid_digit * TableSaveMicro{:, init_i + 17} .* Config.outputResolution^2) / valid_digit; % Surface area 3
    
    %% Save results to Excel file
    if (Config.ReadfromNAS == 1)
        % Save to NAS location
        String_segments = split(analysisfolder, filesep);
        geneName = String_segments{end-1};
        experimentName = String_segments{end};
        rootfolderTmp = strcat(Config.localNASrootPath, filesep, geneName, filesep, experimentName);
        filename = strcat(rootfolderTmp, filesep, 'Betacell_results.xlsx');
    else
        % Save to local folder
        filename = strcat(analysisfolder, filesep, 'Betacell_results.xlsx');
    end
    
    % Write both pixel and micrometer data to separate sheets
    writetable(TableSave, filename, 'Sheet', 'Betacells', 'WriteMode', 'overwritesheet');
    writetable(TableSaveMicro, filename, 'Sheet', 'BetacellsMicrometer', 'WriteMode', 'overwritesheet');
    disp("Beta cells image analysis completed");
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
