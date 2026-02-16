# eeganalysis

*eeganalysis* is an R package for comprehensive EEG data analysis, specializing in BioSemi recordings. The package provides functions to import, preprocess, analyze, and visualize EEG data for research purposes.

Author: Christos Dalamarinis

Contact: [[dalamarinischristos\@gmail.com](mailto:dalamarinischristos@gmail.com)]

<!-- badges: start -->

[![R-CMD-check](https://github.com/ChristosDalamarinis/eeganalysis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ChristosDalamarinis/eeganalysis/actions/workflows/R-CMD-check.yaml) [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental) ![Version](https://img.shields.io/badge/version-0.0.0.9000-blue) ![License: MIT](https://img.shields.io/badge/License-MIT-green.svg) ![Last Commit](https://img.shields.io/github/last-commit/ChristosDalamarinis/eeganalysis)

<!-- badges: end -->

# Version: 0.0.0.9000 (under development)

Current code reads .bdf and .edf files. Scripts include functions for assisting the researcher in importing, down-sampling, and referencing the data.

## Interactive function include:

-   **identify_external_channels()**: Allows the user to interactively label external channels detected in imported EEG data.
-   **plot_electrode_3d()**: 3D visualization of electrode positions on a head model (helper function).
-   **plot_electrode_3d_spherical()**: 3D visualization of electrode positions on a spherical head model (helper function).

## Dependencies

-   "signal" -\> Signal processing and filtering

-   "ggplot2" -\> Data visualization

-   "dplyr" -\> Data manipulation

-   "plotly" -\> Interactive 3D plotting

## Quick Start

``` r
library(eeganalysis)
```

### Import BioSemi data

``` r
eeg_data <- read_biosemi("path/to/your/file.bdf")
```

### Downsample to 256 Hz

``` r
eeg_data <- downsample(eeg_data, target_srate = 256)
```

### Rereference to average

``` r
eeg_data <- eeg_rereference(eeg_data, ref_type = "average")
```

### Extract epochs around events

``` r
epochs <- epoch_eeg(eeg_data, events = c(1, 2), time_window = c(-0.2, 0.8))
```

# Analysis Pipeline

The analysis pipeline implemented in the package includes the following steps:

-   1.Data Import: Functions to import raw EEG data from various file formats (e.g., .edf, .bdf, .set).

-   2.Preprocessing: Functions for filtering, artifact removal, and epoching the EEG data.

-   3.Feature Extraction: Functions to extract relevant features from the EEG signals (e.g., power spectral density, event-related potentials).

-   4.Statistical Analysis: Functions to perform statistical tests on the extracted features.

-   5.Visualization: Functions to create visualizations of the EEG data and analysis results (e.g., topographic maps, time-frequency plots).

## Installation

You can install the *eeganalysis* package from GitHub using the following command:

``` r
# install.packages("devtools") # Uncomment if devtools is not installed
devtools::install_github("dalamarinischristos/eeganalysis")
```

## Contributing

Contributions to the *eeganalysis* package are welcome! If you would like to contribute, please fork the repository and submit a pull request with your changes. Please ensure that your code follows the existing style and includes appropriate documentation and tests. Note: currently focusing on Biosemi files, other formats will be added in future releases.

## Folder structure

``` r
eeganalysis
│
├── R/                                    ← Your R functions live here
│   ├── eeg_class.R                       ← EEG data structure & creation
│   │   ├── new_eeg()                     ← Create an eeg object
│   │   └── print.eeg()                   ← Display eeg object nicely
│   │
│   ├── read_bdf_native.R                 ← BioSemi file import
│   │   └── read_bdf_native()             ← Import .bdf files (main function)
│   │
│   ├── extract_bdf_events()
│   │   ├── extract_bdf_events()          ← Parse trigger codes from status channel
│   │   ├── summary_bdf_events()          ← Summary of extracted BDF events   
│   │   └── validate_bdf_events()         ← Validate and summarize BDF events
│   │
│   ├── channel_info2.R                   ← Electrode database & channel inspection
│   │   ├── get_electrode_database()      ← Access 64-ch BioSemi electrode database
│   │   ├── get_electrode_position()      ← Get coordinates for specific electrode
│   │   ├── inspect_bdf_channels()        ← Preview channels without full import
│   │   ├── classify_channel_naming()     ← Identify naming convention (10-20/10-10/BioSemi)
│   │   ├── plot_electrode_3d()           ← 3D electrode visualization (Cartesian)
│   │   └── plot_electrode_3d_spherical() ← 3D electrode visualization (Spherical)
|   |
|   ├── setexchannels.R                   ← External channel management
│   │   ├── identify_external_channels()  ← Interactive labeling (EOG, EMG, ECG, GSR)
│   │   ├── detect_external_channels()    ← Automated external channel detection
│   │   └── apply_external_labels()       ← Apply user-defined labels to data
│   │
│   ├── downsample.R                      ← Smart downsampling with anti-aliasing
│   │   └── downsample()                  ← Downsample EEG data with filters
|   |
|   ├── rereference.R                     ← Re-referencing utilities
│   │   └── eeg_rereference()             ← Change reference scheme (average/custom)
│   │
│   ├── epoch2.R                          ← Epoching functions
│   │   ├── inspect_triggers()            ← Inspect Event Triggers in EEG Data
│   │   ├── plot_epochs()                 ← Visualize  extracted epochs
│   │   └── epoch_eeg()                   ← Extract time-locked epochs around events
│   │
│   ├── preprocessing.R                   ← (Future) Data cleaning & preprocessing
│   │   ├── filter_eeg()                  ← Apply filters
│   │   ├── rereference_eeg()             ← Change reference scheme
│   │   └── ...                           ← Additional functions
│   │
│   └── feature_extraction.R             ← (Future) Feature computation
│       ├── compute_erp()                ← Calculate ERPs
│       ├── compute_power()              ← Band power analysis
│       └── ...                          ← Additional functions
│
├── man                                  ← Auto-generated help files
│   ├── new_eeg.Rd                       ← Help for new_eeg()
│   ├── print.eeg.Rd                     ← Help for print.eeg()
│   ├── read_biosemi.Rd                  ← Help for read_biosemi()
│   ├── extract_biosemi_events.Rd        ← Help for extract_biosemi_events()
│   ├── summarize_biosemi_import.Rd      ← Help for summarize_biosemi_import()
│   ├── get_electrode_database.Rd        ← Help for get_electrode_database()
│   ├── get_electrode_position.Rd        ← Help for get_electrode_position()
│   ├── inspect_bdf_channels.Rd          ← Help for inspect_bdf_channels()
│   ├── identify_external_channels.Rd    ← Help for identify_external_channels()
│   ├── detect_external_channels.Rd      ← Help for detect_external_channels()
│   ├── apply_external_labels.Rd         ← Help for apply_external_labels()
│   ├── downsample.Rd                    ← Help for downsample()
│   ├── eeg_rereference.Rd               ← Help for eeg_rereference()
│   ├── epoch_eeg.Rd                     ← Help for epoch_eeg()
│   └── ...Rd                            ← Additional help files
|
├── data                                 ← Example datasets
│   └── example_eeg.RData                ← Example eeg object for testing (coming)
│
├── NAMESPACE                            ← Exported functions (auto-generated)
├── DESCRIPTION                          ← Package metadata
├── LICENSE                              ← License information
├── README.md                            ← Package overview
├── .gitignore                           ← Git ignore rules
├── .Rbuildignore                        ← Build ignore rules
└── eeganalysis.Rproj                    ← RStudio project file
```

Note: This is structure is frequently expanded and updated with additional folders/files/functions.
