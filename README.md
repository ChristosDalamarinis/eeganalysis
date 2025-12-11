# eeganalysis - Package under construction

*eeganalysis* is an (under constraction) R package for analyzing EEG data. The scope is to create functions to pre-processes, analyze and visualize EEG data for research purposes.

Author: Christos Dalamarinis

Contact: [[dalmarinischristos\@gmail.com](mailto:dalmarinischristos@gmail.com)]

| Version: 0.1.0

## Analysis Pipeline
The analysis pipeline implemented in the package includes the following steps:
1. Data Import: Functions to import raw EEG data from various file formats (e.g., .edf, .bdf, .set).
2. Preprocessing: Functions for filtering, artifact removal, and epoching the EEG data.
3. Feature Extraction: Functions to extract relevant features from the EEG signals (e.g., power spectral density, event-related potentials).
4. Statistical Analysis: Functions to perform statistical tests on the extracted features.
5. Visualization: Functions to create visualizations of the EEG data and analysis results (e.g., topographic maps, time-frequency plots).
## Installation
You can install the *eeganalysis* package from GitHub using the following command:
```R
# install.packages("devtools") # Uncomment if devtools is not installed
devtools::install_github("dalmarinischristos/eeganalysis")