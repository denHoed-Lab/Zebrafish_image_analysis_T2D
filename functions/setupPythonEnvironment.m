function [pythonPath, success] = setupPythonEnvironment(Config)
% SETUPPYTHONENVIRONMENT - Setup and manage Python environment for VAST Pipeline
%
% Inputs:
%   Config - Configuration structure with Python environment settings
%
% Outputs:
%   pythonPath - Path to Python executable
%   success - Boolean indicating if setup was successful
%
% This function:
% 1. Checks if Python environment exists
% 2. Creates environment if needed
% 3. Installs required packages
% 4. Tests the environment
% 5. Returns Python executable path

success = false;
pythonPath = '';

try
    % Get configuration
    if ~isfield(Config, 'PythonEnv')
        error('Python environment configuration not found in Config.PythonEnv');
    end

    pyConfig = Config.PythonEnv;

    % Check if auto-setup is enabled
    if ~pyConfig.AutoSetup
        % Use existing Python path
        if ~isempty(pyConfig.FallbackPath) && exist(pyConfig.FallbackPath, 'file')
            pythonPath = pyConfig.FallbackPath;
            success = true;
            % fprintf('Using existing Python path: %s\n', pythonPath);  % Silenced debug output
            return;
        else
            error('Auto-setup disabled but no valid Python path provided');
        end
    end

    % Check if conda is available
    [condaStatus, ~] = system('conda --version');
    if condaStatus ~= 0
        % Try common install locations
        basePaths = {
            fullfile(getenv('USERPROFILE'), 'miniconda3');
            fullfile(getenv('USERPROFILE'), 'anaconda3');
            fullfile(getenv('USERPROFILE'), 'Anaconda3');
            fullfile(strcat('C:',filesep,'ProgramData'),'anaconda3');
            };

        found = false;
        for i = 1:numel(basePaths)
            condaDir = basePaths{i};
            condaExe = fullfile(condaDir, 'Scripts', 'anaconda.exe');
            if isfile(condaExe)
                % Temporarily prepend conda paths to PATH
                oldPath = getenv('PATH');
                newPath = [ ...
                    fullfile(condaDir, 'Library', 'bin') ';' ...
                    fullfile(condaDir, 'Scripts') ';' ...
                    condaDir ';' ...
                    oldPath ...
                    ];
                setenv('PATH', newPath);

                % Check if conda works now
                [status, result] = system('conda --version');
                if status == 0
                    disp(['Conda found: ' result]);
                    found = true;
                    break;
                end
            end
        end

        if ~found
            error('Conda not found. Please install Anaconda/Miniconda or ensure it is accessible.');
        end
    end

    % Get current directory for relative paths
    currentDir = pwd;

    % Check if environment already exists - more reliable method
    [envStatus, envOutput] = system(['conda env list']);
    envExists = contains(envOutput, pyConfig.EnvName);

    if envExists && ~pyConfig.ForceReinstall
        % fprintf('Python environment "%s" already exists. Skipping setup.\n', pyConfig.EnvName);  % Silenced debug output

        % Verify the environment is working by checking if we can get Python path
        [pathStatus, pathOutput] = system(['conda run -n ' pyConfig.EnvName ' python -c "import sys; print(sys.executable)"']);
        if pathStatus == 0
            pythonPath = strtrim(pathOutput);

            % Quick check if key packages are installed
            [pkgStatus, pkgOutput] = system(['conda run -n ' pyConfig.EnvName ' python -c "import torch, cv2, numpy; print(''Packages OK'')"']);
            if pkgStatus == 0
                success = true;
                % fprintf('Using existing Python environment: %s\n', pythonPath);  % Silenced debug output

                % Test the environment if in test mode
                if pyConfig.TestMode
                    % fprintf('Testing existing Python environment...\n');  % Silenced debug output
                    testResult = testPythonEnvironment(pythonPath, currentDir);
                    if ~testResult
                        warning('Existing environment test failed. Consider setting ForceReinstall=1');
                    else
                        % fprintf('Existing environment test passed!\n');  % Silenced debug output
                    end
                end
                return;
            else
                % fprintf('Environment exists but missing required packages. Reinstalling packages...\n');  % Silenced debug output
                % Continue to package installation section
            end
        else
            % fprintf('Environment exists but appears corrupted. Recreating...\n');  % Silenced debug output
            envExists = false; % Force recreation
        end
    end

    % Create or recreate environment
    if envExists
        % fprintf('Removing existing environment "%s"...\n', pyConfig.EnvName);  % Silenced debug output
        system(['conda env remove -n ' pyConfig.EnvName ' -y']);
    end

    % Create new environment if needed
    if ~envExists
        % fprintf('Creating Python environment "%s" with Python %s...\n', pyConfig.EnvName, pyConfig.PythonVersion);  % Silenced debug output
        createCmd = sprintf('conda create -n %s python=%s -y', pyConfig.EnvName, pyConfig.PythonVersion);
        [createStatus, createOutput] = system(createCmd);

        if createStatus ~= 0
            error('Failed to create Python environment: %s', createOutput);
        end
    end

    % Install packages from requirements file
    if exist(pyConfig.RequirementsFile, 'file')
        % fprintf('Installing Python packages from %s...\n', pyConfig.RequirementsFile);  % Silenced debug output

        % Get absolute path to requirements file
        if isfile(pyConfig.RequirementsFile)
            reqPath = pyConfig.RequirementsFile;
        else
            reqPath = fullfile(currentDir, pyConfig.RequirementsFile);
        end

        installCmd = sprintf('conda run -n %s pip install -r "%s"', pyConfig.EnvName, reqPath);
        [installStatus, installOutput] = system(installCmd);

        if installStatus ~= 0
            error('Failed to install Python packages: %s', installOutput);
        end
    else
        warning('Requirements file not found: %s', pyConfig.RequirementsFile);
    end

    % Get Python executable path
    [pathStatus, pathOutput] = system(['conda run -n ' pyConfig.EnvName ' python -c "import sys; print(sys.executable)"']);
    if pathStatus ~= 0
        error('Failed to get Python executable path: %s', pathOutput);
    end

    pythonPath = strtrim(pathOutput);

    % Test the environment
    if pyConfig.TestMode
        % fprintf('Testing Python environment...\n');  % Silenced debug output
        testResult = testPythonEnvironment(pythonPath, currentDir);
        if ~testResult
            error('Python environment test failed');
        end
        % fprintf('Python environment test passed!\n');  % Silenced debug output
    end

    success = true;
    % fprintf('Python environment setup completed successfully.\n');  % Silenced debug output
    % fprintf('Python executable: %s\n', pythonPath);  % Silenced debug output

catch ME
    % fprintf('Error in setupPythonEnvironment: %s\n', ME.message);  % Silenced debug output
    if isfield(Config, 'PythonEnv') && Config.PythonEnv.TestMode
        % fprintf('Full error details:\n');  % Silenced debug output
        % disp(getReport(ME, 'extended'));  % Silenced debug output
    end
    success = false;
end

end

function success = testPythonEnvironment(pythonPath, currentDir)
% Test the Python environment by running a simple test script

success = false;

try
    % Create a simple test script
    testScript = fullfile(currentDir, 'test_python_env.py');
    fid = fopen(testScript, 'w');
    if fid == -1
        error('Cannot create test script');
    end

    % Write test script content
    fprintf(fid, 'import sys\n');
    fprintf(fid, 'import torch\n');
    fprintf(fid, 'import cv2\n');
    fprintf(fid, 'import numpy as np\n');
    fprintf(fid, 'import tifffile\n');
    fprintf(fid, 'import scipy\n');
    fprintf(fid, 'import skimage\n');
    fprintf(fid, 'print("All required packages imported successfully")\n');
    fprintf(fid, 'print(f"Python version: {sys.version}")\n');
    fprintf(fid, 'print(f"PyTorch version: {torch.__version__}")\n');
    fprintf(fid, 'print(f"OpenCV version: {cv2.__version__}")\n');
    fclose(fid);

    % Run test script
    testCmd = sprintf('"%s" "%s"', pythonPath, testScript);
    [testStatus, testOutput] = system(testCmd);

    % Clean up test script
    if exist(testScript, 'file')
        delete(testScript);
    end

    if testStatus == 0
        % fprintf('Test output:\n%s\n', testOutput);  % Silenced debug output
        success = true;
    else
        % fprintf('Test failed with output:\n%s\n', testOutput);  % Silenced debug output
        success = false;
    end

catch ME
    % fprintf('Error in testPythonEnvironment: %s\n', ME.message);  % Silenced debug output
    success = false;
end

end