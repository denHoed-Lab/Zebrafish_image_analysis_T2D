# Zebrafish Image Analysis Pipeline for Type 2 Diabetes traits

A comprehensive MATLAB-based image analysis pipeline for zebrafish larvae analysis, combining deep learning segmentation with traditional image processing techniques.

**This software was developed for the research study:**
> "_In vivo_, large-scale perturbation screen of T2D candidate genes highlights potential therapeutic targets"  
> *Mujica and Zhang et al., 2025*

The pipeline enables automated analysis of zebrafish larvae images from VAST (Versatile Automated Screening Technology) systems, providing quantitative measurements for functional genomics studies of type 2 diabetes candidate genes.

## Overview

This pipeline performs automated analysis of zebrafish larvae images from VAST (Versatile Automated Screening Technology) systems, providing:

- **Beta Cell Detection**: Automated detection and quantification of beta cells in zebrafish larvae
- **Liver Analysis**: Liver tissue segmentation and lipid droplet quantification
- **VAST Shape Analysis**: Morphological analysis of brightfield images including body segmentation, eye/ear detection, and developmental assessment

## Features

- **Multi-modal Analysis**: Combines Python deep learning models with MATLAB image processing
- **Automated Environment Setup**: Self-contained Python environment management
- **Comprehensive Results**: Generates quantitative statistics, visualizations, and detailed reports
- **Cross-platform Compatibility**: Tested on Windows systems with MATLAB and Python integration

## System Requirements

### Software Requirements
- **Windows Operating System**
  - Tested on Windows 11
  - Not yet suitable to run on macOS (working on it)
- **MATLAB R2020a or later** with the following toolboxes:
  - Image Processing Toolbox
  - Deep Learning Toolbox
  - Computer Vision Toolbox
  - Statistics and Machine Learning Toolbox
  - Tested on version 2022a, 2022b and 2023b
- **Python 3.9** (automatically managed by the pipeline)
- **Anaconda or Miniconda** for Python environment management

### Hardware Requirements
- **RAM**: Minimum 8GB, recommended 16GB+
- **Storage**: 10GB+ free space for Python environment, data and results
- **GPU**: Optional but recommended for faster deep learning inference

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/denHoed-Lab/Zebrafish_image_analysis_T2D.git
```
### 2. Download test data and weights
1. Open the ./Test data/ and use the path provided in the readme file to download the data from Google Drive
2. Repeat the same for the following paths:
  - ./DNNs/MATLAB weights/
  - ./DNNs/python/weights/

### 3. MATLAB Setup
1. Open MATLAB and navigate to the cloned directory
2. Add the pipeline to MATLAB path:
   ```matlab
   addpath(genpath('.'))
   ```

### 4. Python Environment Setup
1. Run the Python environment setup script:
   ```matlab
   RunPythonSetupfirst
   ```
   This script will:
   - Create a conda environment named `pipeline_python_env`
   - Install required Python packages (PyTorch, OpenCV, NumPy, etc.)
   - Verify the installation

### 4. Verify Installation
Run the test script to ensure everything is working:
```matlab
Analysis_main
```
### 5. Run time
Excluding the data downloads, setting up the python enviroment takes 23.73s

## Usage

### Quick Start
1. **Prepare your data** in the following structure:
   ```
   Test data/
   └── Test_experiment_demo/
       ├── C01/
       │   ├── betacell_data.tif
       │   ├── liver_data.tif
       │   └── VAST/
       │       ├── VAST_1_1.tiff
       │       ├── VAST_1_4.tiff
       │       ├── VAST_1_7.tiff
       │       └── VAST_1_10.tiff
       └── C02/
           └── ...
   ```

2. **Run the analysis**:
   ```matlab
   Analysis_main
   ```

3. **Check results** in the `LocalResults` folder:
   ```
   Test data/Test_experiment_demo/C01/
   ├── LocalResults/
   │   ├── Betacell/
   │   │   ├── Statistics.xlsx
   │   │   ├── Cell Properties.xlsx
   │   │   └── Visualization files
   │   ├── Liver/
   │   │   ├── Statistics.xlsx
   │   │   └── Segmentation files
   │   └── VAST/
   │       ├── VAST_results.xlsx
   │       └── Analysis files
   └── Individual well results
   ```
4. **Run time**:
   On a computer with a Intel(R) Core(TM) i9-9900K CPU @ 3.60GHz processor, 32GB RAM and an 11GB graphics card:
- It takes ~15min to analyze the three samples in the test data.
   
6. **Average time to run the analysis per sample and per model**:
-  Image Type:  Time per Image
-  Beta cells:  ~1 min 37 sec
-  Liver:       ~2 min 13 sec
-  VAST_        ~1 min



### Configuration

The pipeline uses a centralized configuration file (`Config/Config_Demo.m`) where you can modify:

- **Analysis parameters**: Thresholds, kernel sizes, detection parameters
- **File paths**: Input/output directories, network weights locations
- **Python environment**: Environment name, Python version, package requirements

### Data Format Requirements

#### Beta Cell Analysis
- **Format**: Multi-channel TIFF files
- **Naming**: Files containing "betacell" or "lng" in filename
- **Channels**: GFP channel for beta cell detection

#### Liver Analysis
- **Format**: Multi-channel TIFF files
- **Naming**: Files containing "liver" or "lng" in filename
- **Channels**: GFP (liver) and lipid channels

#### VAST Analysis
- **Format**: Brightfield TIFF images
- **Naming**: VAST_1_1.tiff, VAST_1_4.tiff, VAST_1_7.tiff, VAST_1_10.tiff
- **Views**: Dorsal (1,7) and lateral (4,10) views

## Output Files

### Beta Cell Analysis
- **XLSX Files**: Cell counts, volumes, intensities, islet measurements
- **Visualization**: 3D cell segmentation videos, ROI detection images
- **Segmentation**: Binary masks of detected cells

### Liver Analysis
- **XLSX Files**: Liver area, lipid counts, intensity measurements
- **Visualization**: Liver boundary and lipid spot overlays
- **Segmentation**: Liver region and lipid droplet masks

### VAST Analysis
- **XLSX Files**: Body measurements, eye/ear properties, developmental angles
- **Visualization**: Analysis overlays on original images
- **Segmentation**: Body, swim bladder, and organ masks

## Troubleshooting

### Common Issues

1. **Python Environment Not Found**
   - Solution: Run `RunPythonSetupfirst` again
   - Check if conda is properly installed and in PATH

2. **MATLAB Toolbox Missing**
   - Solution: Install required MATLAB toolboxes
   - Check MATLAB version compatibility

3. **Memory Issues**
   - Solution: Reduce image size or batch processing
   - Close other applications to free memory

4. **File Path Errors**
   - Solution: Ensure data follows the expected folder structure
   - Check file naming conventions

### Getting Help

- Check the log files generated in each analysis folder
- Review the configuration parameters in `Config/Config_Demo.m`
- Ensure all dependencies are properly installed

## Citation

If you use this pipeline in your research, please cite our work:

### Primary Citation (Research Paper)
```bibtex
@article{mujica_t2d_2025,
  title={In vivo, large-scale perturbation screen of T2D candidate genes highlights potential therapeutic targets},
  author={Mujica, Endrina and Zhang, Hanqing and Emmanouilidou, Anastasia and Cook, Naomi and Metzendorf, Christoph and Mazzaferro, Eugenia and Bandaru, Manoj and Costa, Joao Campos and Alavioon, Ghazal and Klingström, Tiffany and Zalamitsou, Chrysoula and van Zuydam, Natalie and Dronkers, Tessa and Hariri, Mehran and Mathews, Bobby and Pakula, Adrianna and Rottner, Antje K and Sun, Han and Knowles, Joshua W and Ingelsson, Erik and McCarthy, Mark I and van de Bunt, Martijn and Frederiksen, Klaus Stensgaard and Gloyn, Anna L and Brooke, Hannah L and Larsson, Anders and Vienberg, Sara Gry and Flannick, Jason and Allalou, Amin and den Hoed, Marcel},
  journal={[Journal Name]},
  year={2025},
  volume={[Volume]},
  number={[Number]},
  pages={[Pages]},
  doi={[DOI]},
  publisher={[Publisher]}
}
```

## Related Publications

If you use this pipeline in your research, please also cite the underlying technologies and methods:

### MATLAB
```
MATLAB. (2020). Version R2020b. Natick, Massachusetts: The MathWorks Inc.
```

### Python Libraries
```
Harris, C.R., Millman, K.J., van der Walt, S.J. et al. (2020). Array programming with NumPy. Nature 585, 357–362.
```

```
Bradski, G. (2000). The OpenCV Library. Dr. Dobb's Journal of Software Tools.
```

```
Paszke, A., Gross, S., Massa, F. et al. (2019). PyTorch: An Imperative Style, High-Performance Deep Learning Library. Advances in Neural Information Processing Systems 32.
```

## Acknowledgments

When using this pipeline, please acknowledge:

- The VAST (Versatile Automated Screening Technology) system developers
- The MATLAB and Python communities for the underlying libraries
- Research collaborators and contributors to this project
- The Beijer Laboratory and Department of Immunology, Genetics and Pathology at Uppsala University
- SciLifeLab for computational resources and support

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Contact

- **Author**: Endrina Mujica
- **Email**: endrina.mujica@igp.uu.se
- **Institution**: Uppsala University, Sweden

- **Author**: Hanqing Zhang
- **Email**: hanzha@kth.se
- **Institution**: KTH Royal Institute of Technology

- **Author**: Amin Allalou
- **Email**: amin.allalou@it.uu.se
- **Institution**: Uppsala University, Sweden

- **Author**: Marcel den Hoed
- **Email**: marcel.den_hoed@igp.uu.se
- **Institution**: Uppsala University, Sweden 
