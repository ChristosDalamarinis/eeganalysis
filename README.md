# eeganalysis - Package under construction

*eeganalysis* is an (under constraction) R package for analyzing EEG data. The scope is to create functions to pre-processes, analyze and visualize EEG data for research purposes.

Author: Christos Dalamarinis

Contact: [[dalmarinischristos\@gmail.com](mailto:dalmarinischristos@gmail.com)]

| Version: 0.0.1

## Analysis Pipeline

The analysis pipeline implemented in the package includes the following steps:

-   1.Data Import: Functions to import raw EEG data from various file formats (e.g., .edf, .bdf, .set).

-   2.Preprocessing: Functions for filtering, artifact removal, and epoching the EEG data.

-   3.Feature Extraction: Functions to extract relevant features from the EEG signals (e.g., power spectral density, event-related potentials).

-   4.Statistical Analysis: Functions to perform statistical tests on the extracted features.

-   5.Visualization: Functions to create visualizations of the EEG data and analysis results (e.g., topographic maps, time-frequency plots).

## Installation You can install the *eeganalysis* package from GitHub using the following command:

``` r
# install.packages("devtools") # Uncomment if devtools is not installed
devtools::install_github("dalmarinischristos/eeganalysis")
```

## Contributing

Contributions to the *eeganalysis* package are welcome! If you would like to contribute, please fork the repository and submit a pull request with your changes. Please ensure that your code follows the existing style and includes appropriate documentation and tests. Note: currently focusing on Biosemi files, other formats will be added in future releases.

## Folder structure
``` r
eeganalysis/
│
├── R/                              ← Your R functions live here
│   ├── 01_eeg_class.R             ← EEG data structure & creation
│   │   ├── new_eeg()              ← Create an eeg object
│   │   └── print.eeg()            ← Display eeg object nicely
│   │
│   ├── 02_read_biosemi.R          ← BioSemi file import
│   │   ├── read_biosemi()         ← Import .bdf files (main function)
│   │   ├── extract_biosemi_events() ← Parse trigger codes from status channel
│   │   └── summarize_biosemi_import() ← Generate quality report
│   │
│   ├── 03_preprocessing.R         ← (Future) Data cleaning & preprocessing
│   │   ├── filter_eeg()           ← Apply filters
│   │   ├── rereference_eeg()      ← Change reference scheme
│   │   └── ...
│   │
│   └── 04_feature_extraction.R    ← (Future) Feature computation
│       ├── compute_erp()          ← Calculate ERPs
│       ├── compute_power()        ← Band power analysis
│       └── ...
│
├── man/                           ← Auto-generated help files
│   ├── new_eeg.Rd                ← Help for new_eeg()
│   ├── print.eeg.Rd              ← Help for print.eeg()
│   ├── read_biosemi.Rd           ← Help for read_biosemi()
│   ├── extract_biosemi_events.Rd ← Help for extract_biosemi_events()
│   ├── summarize_biosemi_import.Rd ← Help for summarize_biosemi_import()
│   └── ...
│
├── tests/                        ← Unit tests for your functions
│   └── testthat/
│       ├── test_eeg_class.R      ← Tests for eeg_class functions
│       ├── test_read_biosemi.R   ← Tests for biosemi reader
│       └── ...
│
├── data/                         ← Example datasets
│   └── example_eeg.RData         ← Example eeg object for testing
│
├── vignettes/                    ← User tutorials/guides
│   ├── getting_started.Rmd       ← How to get started
│   └── preprocessing_guide.Rmd   ← Preprocessing tutorial
│
├── NAMESPACE                     ← Exported functions (auto-generated)
├── DESCRIPTION                   ← Package metadata
├── LICENSE                       ← License information
├── README.md                     ← Package overview
├── .gitignore                    ← Git ignore rules
├── .Rbuildignore                 ← Build ignore rules
└── eeganalysis.Rproj            ← RStudio project file
``` r
Note: This is a simplified structure and can be expanded with additional folders/files as needed.