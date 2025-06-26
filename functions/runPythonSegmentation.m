function [success, output] = runPythonSegmentation(pythonScript, inputFile, outputFolder, Config)
% RUNPYTHONSEGMENTATION - Execute Python segmentation script with environment management
%
% Inputs:
%   pythonScript - Path to Python script to execute
%   inputFile - Input file path
%   outputFolder - Output folder path
%   Config - Configuration structure
%
% Outputs:
%   success - Boolean indicating if execution was successful
%   output - Command output (for debugging)

success = false;
output = '';

try
    % Setup Python environment if needed
    if isfield(Config, 'PythonEnv')
        [pythonPath, envSuccess] = setupPythonEnvironment(Config);
        if ~envSuccess
            error('Failed to setup Python environment');
        end
    else
        % Fallback to existing Python path from Config (for backward compatibility)
        if isfield(Config, 'pythonEnvPath')
            pythonPath = Config.pythonEnvPath;
            % Remove quotes if present
            pythonPath = strrep(pythonPath, '"', '');
            % fprintf('Using fallback Python path: %s\n', pythonPath);  % Silenced debug output
        else
            error('No Python path configured. Please set up Config.PythonEnv or Config.pythonEnvPath');
        end
    end
    
    % Validate inputs
    if ~exist(pythonScript, 'file')
        error('Python script not found: %s', pythonScript);
    end
    
    if ~exist(inputFile, 'file')
        error('Input file not found: %s', inputFile);
    end
    
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
    % Construct command - clean all paths thoroughly
    pythonScript = strrep(pythonScript, '"', '');
    inputFile = strrep(inputFile, '"', '');
    outputFolder = strrep(outputFolder, '"', '');
    
    % Normalize paths for cross-platform compatibility
    % Convert backslashes to forward slashes and remove trailing slashes
    inputFile = strrep(inputFile, '\', '/');
    outputFolder = strrep(outputFolder, '\', '/');
    
    % Remove trailing slashes to avoid path issues
    while endsWith(inputFile, '/')
        inputFile = inputFile(1:end-1);
    end
    while endsWith(outputFolder, '/')
        outputFolder = outputFolder(1:end-1);
    end
    
    % Ensure output folder exists
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
    % Build command string
    cmd = sprintf('"%s" "%s" "%s" "%s"', pythonPath, pythonScript, inputFile, outputFolder);
    
    % Execute command
    % fprintf('Executing Python segmentation...\n');  % Silenced debug output
    % fprintf('Command: %s\n', cmd);  % Silenced debug output
    
    [status, output] = system(cmd);
    
    if status == 0
        success = true;
        % fprintf('Python segmentation completed successfully.\n');  % Silenced debug output
    else
        error('Python segmentation failed with status %d:\n%s', status, output);
    end
    
    % Check if output files were created - handle different output patterns
    [~, inputName, ~] = fileparts(inputFile);
    
    % Define possible output file patterns based on the Python script
    possibleOutputs = {
        fullfile(outputFolder, [inputName '_bw.tif']),      % Original expected pattern
        fullfile(outputFolder, [inputName '_seg.tif']),     % Common segmentation pattern
        fullfile(outputFolder, [inputName '_seg.png']),     % PNG segmentation pattern
        fullfile(outputFolder, [inputName '_bw.png']),      % PNG binary pattern
        fullfile(outputFolder, [inputName '_result.tif']),  % Result pattern
        fullfile(outputFolder, [inputName '_result.png'])   % PNG result pattern
    };
    
    % Check if any of the expected output files exist
    outputFileFound = false;
    for i = 1:length(possibleOutputs)
        if exist(possibleOutputs{i}, 'file')
            outputFileFound = true;
            break;
        end
    end
    
    if ~outputFileFound
        warning('Expected output file not found for input: %s', inputName);
        warning('Checked for: %s', strjoin(possibleOutputs, ', '));
        % List files in output directory for debugging
        outputFiles = dir(outputFolder);
        % fprintf('Files in output directory:\n');  % Commented out to silence terminal output
        % for i = 1:length(outputFiles)
        %     if ~outputFiles(i).isdir
        %         fprintf('  %s\n', outputFiles(i).name);
        %     end
        % end
    end
    
catch ME
    % fprintf('Error in runPythonSegmentation: %s\n', ME.message);  % Silenced debug output
    if isfield(Config, 'PythonEnv') && Config.PythonEnv.TestMode
        % fprintf('Full error details:\n');  % Silenced debug output
        % disp(getReport(ME, 'extended'));  % Silenced debug output
    end
    success = false;
    output = ME.message;
end

end 