%% Test Python Environment Setup for VAST Pipeline
% This script tests and sets up the Python environment required for the VAST Pipeline
% Run this script before running the main analysis to ensure Python is properly configured

clear all; close all;
cPath=mfilename('fullpath');
[GUIpath,~]=fileparts(cPath);
lPath_f={'functions'};
for lPath_i=lPath_f
    lPath_i=lPath_i{1};
    lPath_i=strcat(GUIpath,'\',lPath_i);
    addpath(genpath(lPath_i))
end
fprintf('=== VAST Pipeline Python Environment Test ===\n\n');
% Create a minimal configuration for testing
Config.PythonEnv.AutoSetup = 1;
Config.PythonEnv.EnvName = 'pipeline_python_env';
Config.PythonEnv.PythonVersion = '3.9';
Config.PythonEnv.RequirementsFile = '.\DNNs\python\requirements.txt';
Config.PythonEnv.ForceReinstall = 0;
Config.PythonEnv.TestMode = 1; % Enable test mode for verbose output
Config.PythonEnv.FallbackPath = '';
Config.PythonEnv.LogFile = 'python_env_setup.log';

% Test Python environment setup
fprintf('Testing Python environment setup...\n');
[pythonPath, success] = setupPythonEnvironment(Config);

if success
    fprintf('\n✅ Python environment setup successful!\n');
    fprintf('Python executable: %s\n', pythonPath);
    
    % Test importing required modules
    fprintf('\nTesting Python module imports...\n');
    testModules(pythonPath);
    
    % Test running a simple segmentation script
    fprintf('\nTesting segmentation script execution...\n');
    testSegmentationScript(pythonPath);
    
    fprintf('\n🎉 All tests passed! Your Python environment is ready for VAST Pipeline.\n');
else
    fprintf('\n❌ Python environment setup failed!\n');
    fprintf('Please check the error messages above and ensure:\n');
    fprintf('1. Anaconda or Miniconda is installed\n');
    fprintf('2. You have internet connectivity for package installation\n');
    fprintf('3. You have sufficient disk space for the Python environment\n');
end

function testModules(pythonPath)
    % Test importing required Python modules
    modules = {'torch', 'cv2', 'numpy', 'tifffile', 'scipy', 'skimage'};
    
    for i = 1:length(modules)
        module = modules{i};
        testCmd = sprintf('"%s" -c "import %s; print(''%s imported successfully'')"', pythonPath, module, module);
        [status, output] = system(testCmd);
        
        if status == 0
            fprintf('✅ %s\n', output);
        else
            fprintf('❌ Failed to import %s\n', module);
        end
    end
end

function testSegmentationScript(pythonPath)
    % Test running a simple segmentation script
    try
        % Create a simple test image
        testImage = rand(100, 100, 10, 'single') * 255;
        testImage = uint8(testImage);
        
        % Save test image
        testInputFile = 'test_input.tiff';
        testOutputFolder = 'test_output';
        
        if ~exist(testOutputFolder, 'dir')
            mkdir(testOutputFolder);
        end
        
        % Save as TIFF
        tiffwrite(testInputFile, testImage);
        
        % Test running the beta cell segmentation script
        scriptPath = '.\DNNs\python\segment_beta_cells.py';
        if exist(scriptPath, 'file')
            cmd = sprintf('"%s" "%s" "%s" "%s"', pythonPath, scriptPath, testInputFile, testOutputFolder);
            [status, output] = system(cmd);
            
            if status == 0
                fprintf('✅ Segmentation script executed successfully\n');
            else
                fprintf('❌ Segmentation script failed:\n%s\n', output);
            end
        else
            fprintf('⚠️  Segmentation script not found: %s\n', scriptPath);
        end
        
        % Clean up test files
        if exist(testInputFile, 'file')
            delete(testInputFile);
        end
        if exist(testOutputFolder, 'dir')
            rmdir(testOutputFolder, 's');
        end
        
    catch ME
        fprintf('❌ Error testing segmentation script: %s\n', ME.message);
    end
end

function tiffwrite(filename, data)
    % Simple TIFF writer for testing
    t = Tiff(filename, 'w');
    tagstruct.ImageLength = size(data, 1);
    tagstruct.ImageWidth = size(data, 2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 8;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';
    
    for i = 1:size(data, 3)
        t.setTag(tagstruct);
        t.write(data(:,:,i));
        if i < size(data, 3)
            t.writeDirectory();
        end
    end
    t.close();
end 