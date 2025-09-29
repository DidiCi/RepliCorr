## RepliCorr
**RepliCorr** is a MATLAB-based analysis tool designed to investigate DNA replication dynamics by examining the spatial correlation of replication patterns along DNA fibers. This methodology, introduced in the study by [Ciardo et al. (2025)](https://academic.oup.com/nar/article/53/3/gkaf007/7990345?login=false), enables the assessment of replication fork speeds and initiation rates through autocorrelation profiling.
The tool processes raw DNA fiber data, computes autocorrelation profiles, and applies statistical models to interpret replication dynamics.


## ⚙️ Installation

### Prerequisites

- MATLAB (version 2020b or later recommended)
- Required MATLAB toolboxes:
  - Global Optimization Toolbox
  - System Identification Toolbox
  - Statistics and Machine Learning Toolbox


### Steps

1. Clone the repository:

   ```bash
   git clone https://github.com/DidiCi/RepliCorr.git
   cd RepliCorr
   ```

2. Add the repository to your MATLAB path:

   ```matlab
   addpath(genpath('path_to_RepliCorr'));
   ```

3. Ensure all dependencies are met by checking the required toolboxes in MATLAB.

## Repository Structure

```
RepliCorr/
├── 1-Data_extraction/       # Scripts for importing and processing raw data
├── 2-Plot_autocorrelation/   # Functions to visualize autocorrelation profiles
├── 3-PCA_and_clustering/     # PCA and clustering analysis tools
├── 4-Fit_single_process/     # Model fitting for single-process replication
├── 5-Fit_two_processes/      # Model fitting for dual-process replication
├── 6-Analysis_fit_oneprocess/ # Analysis of single-process fit results
├── 7-Analysis_fit_twoprocesses/ # Analysis of dual-process fit results
├── 8-Analysis_of_fit_clustering/ # Analysis of clustering results
├── 9-Comparison_control_treatment/ # Scripts for comparing experimental conditions
├── Data_demo/                # Example datasets for demonstration
├── Functions/                # Utility functions used across scripts
├── LICENSE                   # MIT license file
├── README.md                 # Project documentation
```

## 📄 License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.
