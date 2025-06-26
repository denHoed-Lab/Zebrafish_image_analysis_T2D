# Python Environment Setup for VAST Pipeline

This directory contains Python scripts and neural networks for image segmentation in the VAST Pipeline.

## Overview

The VAST Pipeline uses Python-based deep learning models for:
- Beta cell segmentation (`segment_beta_cells.py`)
- Liver segmentation (`LiverSegmentation.py`)
- Liver lipid segmentation (`LiverLipid_1ch_segmentation.py`)
- Fish analysis (`FishAnalysis_Image.py`)

## Automatic Setup (Recommended)

The pipeline now includes automatic Python environment management. Simply run:

```matlab
% Test and setup Python environment
testPythonSetup

% Or run the main analysis (will auto-setup if needed)
Analysis_main
```

## Manual Setup

If you prefer manual setup or need to troubleshoot:

### Prerequisites

1. **Anaconda or Miniconda** - Download from [Anaconda](https://www.anaconda.com/products/distribution) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

2. **Internet Connection** - Required for downloading Python packages

### Setup Steps

1. **Create Python Environment**:
   ```bash
   conda create -n pipeline_python_env python=3.9 -y
   conda activate pipeline_python_env
   ```

2. **Install Dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

3. **Test Installation**:
   ```bash
   python -c "import torch, cv2, numpy, tifffile, scipy, skimage; print('All packages imported successfully')"
   ```

## Configuration

The Python environment can be configured in `Analysis_main.m`:

```matlab
% Python Environment Management (Global for all Python scripts)
Config.PythonEnv.AutoSetup = 1;           % 1: auto-setup, 0: use existing
Config.PythonEnv.EnvName = 'pipeline_python_env'; % Environment name
Config.PythonEnv.PythonVersion = '3.9';   % Python version
Config.PythonEnv.RequirementsFile = '.\DNNs\python\requirements.txt';
Config.PythonEnv.ForceReinstall = 0;      % 1: force reinstall packages
Config.PythonEnv.TestMode = 0;            % 1: verbose output for debugging
```

**Note**: This configuration is global and applies to all Python scripts in the pipeline (Beta cells, Liver, VAST, etc.).

## Required Python Packages

- **torch** (>=1.9.0) - PyTorch deep learning framework
- **torchvision** (>=0.10.0) - Computer vision utilities
- **opencv-python** (>=4.5.0) - Image processing
- **numpy** (>=1.21.0) - Numerical computing
- **tifffile** (>=2021.7.0) - TIFF file handling
- **scipy** (>=1.7.0) - Scientific computing
- **scikit-image** (>=0.18.0) - Image processing algorithms
- **gryds** (>=0.0.3) - Image transformations

## Troubleshooting

### Common Issues

1. **Conda not found**:
   - Install Anaconda/Miniconda
   - Add conda to system PATH
   - Restart MATLAB

2. **Package installation fails**:
   - Check internet connection
   - Try updating pip: `pip install --upgrade pip`
   - Check Python version compatibility

3. **CUDA/GPU issues**:
   - The scripts are configured to use CPU by default
   - For GPU support, install CUDA-enabled PyTorch

4. **Memory issues**:
   - Reduce batch size in segmentation scripts
   - Process smaller image regions

### Debug Mode

Enable debug mode for detailed output:

```matlab
Config.PythonEnv.TestMode = 1;
```

### Manual Python Path

If auto-setup fails, specify a manual Python path:

```matlab
Config.PythonEnv.AutoSetup = 0;
Config.PythonEnv.FallbackPath = 'C:/path/to/your/python.exe';
```

## File Structure

```
DNNs/python/
├── requirements.txt          # Python dependencies
├── segment_beta_cells.py     # Beta cell segmentation
├── LiverSegmentation.py      # Liver segmentation
├── LiverLipid_1ch_segmentation.py  # Liver lipid segmentation
├── FishAnalysis_Image.py     # Fish analysis
├── Augment.py               # Image augmentation utilities
├── modeling/                # Neural network architectures
│   ├── deeplab.py
│   ├── Unet2D.py
│   ├── Unet3D.py
│   └── ...
├── utils/                   # Utility functions
│   ├── Augment.py
│   └── implot.py
└── weights/                 # Pre-trained model weights
    ├── Beta_cells_weights_2022-02-17.pth
    ├── Liver_weights_zstack_v1.10_20210302.pth
    └── ...
```

## Testing

Run the test script to verify your setup:

```matlab
testPythonSetup
```

This will:
1. Create/verify Python environment
2. Install required packages
3. Test module imports
4. Test segmentation script execution

## Support

For issues with Python environment setup:
1. Check the error messages in MATLAB console
2. Enable debug mode: `Config.PythonEnv.TestMode = 1`
3. Check the log file: `python_env_setup.log`
4. Ensure all prerequisites are installed 