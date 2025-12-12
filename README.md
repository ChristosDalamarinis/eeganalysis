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
│
├── R/ ← Your functions live here
│ ├── eeg_class.R ← Defines eeg data structure
│ │ ├── new_eeg() ← Create EEG object 
│ │ └── print.eeg() ← Display EEG object 
│ │
│ └── read_biosemi.R ← Read BioSemi BDF files
│ ├── read_biosemi() ← Main function to import .bdf files 
│ ├── extract_biosemi_events() ← Parse event triggers 
│ └── summarize_biosemi_import() ← Quality report 
│ 
├── man/ ← Auto-generated help documentation
│ ├── new_eeg.Rd 
│ ├── print.eeg.Rd 
│ ├── read_biosemi.Rd
│ └── ... (more coming)
│ 
├── NAMESPACE ← What functions are exported 
├── DESCRIPTION ← Package metadata 
├── LICENSE ← License information 
└── README.md ← This file
