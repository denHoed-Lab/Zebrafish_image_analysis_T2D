# Detailed Usage Guide

This guide explains how to use the VAST Image Analysis Pipeline for analyzing zebrafish embryo images in detail

## Quick Start

### 1. Basic Usage

The simplest way to run the pipeline is:

```matlab
% Add pipeline to path (if not already done)
addpath(genpath('.'));

% Run analysis on test data
Analysis_main;
```

This will:
- Process all samples in the `Test data/Test_experiment_demo/` folder
- Perform beta cell, liver, and VAST analysis
- Generate results in `LocalResults` subfolders
- Create summary Excel files

### 2. Custom Data Analysis

To analyze your own data:

```matlab
% Modify the source path in Analysis_main.m
source_path{1} = '.\Your_Data_Folder\Your_Experiment';

% Or modify directly in the script
Config = Config_Demo();
inputFolder = '.\Your_Data_Folder\Your_Experiment';

% Run individual analyses
BetalogFile = betaCells_framework(inputFolder, Config);
LiverlogFile = Liver_framework(inputFolder, Config);
VAST_process(inputFolder, Config);
```

## Data Preparation

### Required Data Structure

Your data must follow this structure:

```
Your_Experiment_Folder/
├── Sample_01/
│   ├── betacell_data.tif          # Beta cell fluorescence data
│   ├── liver_data.tif             # Liver fluorescence data
│   └── VAST/
│       ├── VAST_1_1.tiff          # Dorsal view 1
│       ├── VAST_1_4.tiff          # Lateral view 1
│       ├── VAST_1_7.tiff          # Dorsal view 2
│       └── VAST_1_10.tiff         # Lateral view 2
├── Sample_02/
│   └── ...
└── Sample_03/
    └── ...
```

### File Naming Conventions

#### Beta Cell Data
- **Required**: Files containing "betacell" or "lng" in the filename
- **Format**: Multi-channel TIFF files
- **Channels**: GFP channel for beta cell detection

#### Liver Data
- **Required**: Files containing "liver" or "lng" in the filename
- **Format**: Multi-channel TIFF files
- **Channels**: 
  - Channel 1: GFP (liver tissue)
  - Channel 2: Lipid fluorescence

#### VAST Data
- **Required**: Files named exactly `VAST_1_1.tiff`, `VAST_1_4.tiff`, `VAST_1_7.tiff`, `VAST_1_10.tiff`
- **Format**: Brightfield TIFF images
- **Views**:
  - `VAST_1_1.tiff`: Dorsal view 1
  - `VAST_1_4.tiff`: Lateral view 1
  - `VAST_1_7.tiff`: Dorsal view 2
  - `VAST_1_10.tiff`: Lateral view 2

## Configuration

### Modifying Analysis Parameters

Edit `Config/Config_Demo.m` to customize analysis parameters:

```matlab
% Beta cell analysis parameters
Config.Betacell.ROISize = 600;                    % ROI size for analysis
Config.Betacell.cellsize_minBetaCellVolume = 1000; % Minimum cell volume
Config.Betacell.PropStat_Diameter = 15;           % Minimum cell diameter

% Liver analysis parameters
Config.Liver.LipidBG_Factor = 3;                  % Lipid background factor
Config.Liver.relativeInt_bias = 1.5;              % Relative intensity bias

% VAST analysis parameters
Config.VAST.umConv = 5600/1024;                   % Pixel to micrometer conversion
Config.VAST.PigmentThreshold = 0.7;               % Pigment detection threshold
```

### Output Configuration

```matlab
% Output folder structure
Config.Betacell.outputFolderlist = {'LocalResults\Betacell'};
Config.Liver.outputFolderlist = {'LocalResults\Liver'};
Config.VAST.outputFolderlist = {'LocalResults\VAST'};

% File overwrite settings
Config.overwriteResults = 0;  % 0: skip existing, 1: overwrite
Config.skipOverwrite = 1;     % 0: overwrite missing, 1: skip missing
```

## Running Individual Analyses

### Beta Cell Analysis Only

```matlab
Config = Config_Demo();
inputFolder = '.\Your_Data_Folder\Your_Experiment';

% Run beta cell analysis
BetalogFile = betaCells_framework(inputFolder, Config);
```

### Liver Analysis Only

```matlab
Config = Config_Demo();
inputFolder = '.\Your_Data_Folder\Your_Experiment';

% Run liver analysis
LiverlogFile = Liver_framework(inputFolder, Config);
```

### VAST Analysis Only

```matlab
Config = Config_Demo();
inputFolder = '.\Your_Data_Folder\Your_Experiment';

% Run VAST analysis
VAST_process(inputFolder, Config);
```

## Understanding Results

### Beta Cell Analysis Results

**Location**: `LocalResults\Betacell\`

**Files Generated**:
- `Statistics.xlsx`: Cell counts, volumes, intensities
- `Cell Properties.xlsx`: Individual cell measurements
- `*_Demo.png`: Visualization images
- `*_seg.tif`: Segmentation masks

**Key Metrics**:
- `TotalCells`: Number of detected beta cells
- `TotalVolume`: Combined volume of all cells
- `AverageIntensity`: Mean fluorescence intensity
- `IsletCount`: Number of islet clusters

### Liver Analysis Results

**Location**: `LocalResults\Liver\`

**Files Generated**:
- `Statistics.xlsx`: Liver and lipid measurements
- `*_Demo.png`: Liver boundary and lipid spot overlays
- `*_seg.tif`: Liver segmentation masks

**Key Metrics**:
- `LiverArea`: Total liver tissue area
- `LipidArea`: Total lipid droplet area
- `NrLipids`: Number of lipid droplets
- `AverageLipidsIntensity`: Mean lipid fluorescence

### VAST Analysis Results

**Location**: `LocalResults\VAST\`

**Files Generated**:
- `VAST_results.xlsx`: Morphological measurements
- `*_seg.tif`: Body and organ segmentation masks
- `*_res.tif`: Analysis result overlays

**Key Metrics**:
- `BodyLength`: Fish body length in micrometers
- `EyeArea`: Total eye area
- `EarArea`: Total ear area
- `DevelopmentalAngle`: Developmental stage angle

## Advanced Usage

### Batch Processing Multiple Experiments

```matlab
% Process multiple experiments
experiments = {'Experiment_1', 'Experiment_2', 'Experiment_3'};

for i = 1:length(experiments)
    inputFolder = fullfile('.\Data', experiments{i});
    fprintf('Processing %s...\n', experiments{i});
    
    % Run all analyses
    BetalogFile = betaCells_framework(inputFolder, Config);
    LiverlogFile = Liver_framework(inputFolder, Config);
    VAST_process(inputFolder, Config);
end
```

### Custom Configuration for Different Data Types

```matlab
% Create custom configuration
Config = Config_Demo();

% Adjust for different image resolutions
Config.resolutionXYZ = [0.5, 0.5, 1.0];  % Your image resolution
Config.VAST.umConv = 5000/1024;          % Adjust conversion factor

% Adjust for different cell types
Config.Betacell.cellsize_minBetaCellVolume = 500;  % Smaller cells
Config.Betacell.PropStat_Diameter = 10;            % Smaller diameter

% Run analysis with custom settings
Analysis_main;
```

### Error Handling and Logging

```matlab
% Enable detailed logging
Config.Betacell.Debug = 1;
Config.Liver.Debug = 1;

% Run with error handling
try
    BetalogFile = betaCells_framework(inputFolder, Config);
    fprintf('Beta cell analysis completed successfully\n');
catch ME
    fprintf('Beta cell analysis failed: %s\n', ME.message);
    % Continue with other analyses
end
```

## Best Practices

### 1. Data Organization
- Keep consistent naming conventions
- Organize data by experiment and sample
- Use descriptive folder names

### 2. Configuration Management
- Create backup copies of working configurations
- Document parameter changes
- Test configurations on small datasets first

### 3. Quality Control
- Review segmentation results visually
- Check for edge cases (touching boundaries, etc.)
- Validate measurements against known standards

### 4. Performance Optimization
- Use appropriate image sizes for your analysis
- Close other applications during processing
- Consider using GPU acceleration if available

### 5. Result Validation
- Compare results across multiple runs
- Check for consistency in measurements
- Validate against manual measurements when possible

## Troubleshooting

### Common Issues

1. **No Results Generated**
   - Check file naming conventions
   - Verify data folder structure
   - Review log files for errors

2. **Poor Segmentation Quality**
   - Adjust threshold parameters
   - Check image quality and contrast
   - Verify channel assignments

3. **Memory Issues**
   - Reduce image size in configuration
   - Process fewer images at once
   - Close other applications

4. **Python Errors**
   - Verify Python environment setup
   - Check package installations
   - Review Python script paths

### Getting Help

- Check log files in analysis folders
- Review MATLAB error messages
- Verify configuration parameters
- Contact: hanzha@kth.se

## Examples

### Example 1: Standard Analysis
```matlab
% Standard analysis workflow
addpath(genpath('.'));
Config = Config_Demo();
Analysis_main;
```

### Example 2: Custom Data Path
```matlab
% Analyze custom data location
Config = Config_Demo();
source_path{1} = '.\My_Experiment_Data';
inputFolder = source_path{1};

BetalogFile = betaCells_framework(inputFolder, Config);
LiverlogFile = Liver_framework(inputFolder, Config);
VAST_process(inputFolder, Config);
```

### Example 3: Parameter Optimization
```matlab
% Optimize parameters for specific data
Config = Config_Demo();

% Adjust for high-resolution data
Config.Betacell.ROISize = 800;
Config.Betacell.cellsize_minBetaCellVolume = 2000;

% Adjust for different lipid staining
Config.Liver.LipidBG_Factor = 2.5;
Config.Liver.relativeInt_bias = 1.2;

% Run analysis
Analysis_main;
``` 
