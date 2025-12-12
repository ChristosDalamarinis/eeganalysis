# eeganalysis - R Package for EEG Analysis

*eeganalysis* is an R package for reading, preprocessing, and analyzing EEG data from BioSemi recordings. The package provides tools to import raw EEG data, clean and preprocess it, extract features, and perform statistical analysis.

**Author:** Christos Dalamarinis  
**Contact:** [dalmarinischristos@gmail.com](mailto:dalmarinischristos@gmail.com)  
**Status:** Under Construction  
**Current Version:** 0.0.1

---

## Package Structure

```
eeganalysis (your package)
│
├── R/                           ← Your functions live here
│   ├── eeg_class.R              ← Defines eeg data structure
│   │   ├── new_eeg()            ← Create EEG object
│   │   └── print.eeg()          ← Display EEG object
│   │
│   └── read_biosemi.R           ← Read BioSemi BDF files
│       ├── read_biosemi()       ← Main function to import .bdf files
│       ├── extract_biosemi_events()  ← Parse event triggers
│       └── summarize_biosemi_import()  ← Quality report
│
├── man/                         ← Auto-generated help documentation
│   ├── new_eeg.Rd
│   ├── print.eeg.Rd
│   ├── read_biosemi.Rd
│   └── ...
│
├── NAMESPACE                    ← What functions are exported
├── DESCRIPTION                  ← Package metadata
├── LICENSE                      ← License information
└── README.md                    ← This file
```

---

## Analysis Pipeline

The package implements a complete EEG analysis pipeline:

### 1. **Data Import** (Currently Implemented)
   - Functions to import raw EEG data from BioSemi BDF files
   - Automatic extraction of channel names, sampling rate, and metadata
   - Event/trigger marker detection

### 2. **Preprocessing** (Coming Soon)
   - Filtering (high-pass, low-pass, notch filters)
   - Re-referencing (average, linked mastoids, etc.)
   - Artifact handling (ICA, automatic detection)
   - Epoching (segmentation around events)
   - Baseline correction

### 3. **Feature Extraction** (Coming Soon)
   - Event-Related Potentials (ERPs)
   - Time-frequency analysis (spectral power, wavelets)
   - Phase-based measures (connectivity, phase-locking)
   - Single-trial features

### 4. **Statistical Analysis** (Coming Soon)
   - Within-subject statistics
   - Group-level models (ANOVA, mixed-effects models)
   - Permutation tests
   - Multiple comparison corrections (FDR, cluster-based)

### 5. **Visualization** (Coming Soon)
   - ERP plots with confidence intervals
   - Topographic maps
   - Time-frequency plots
   - Connectivity networks

---

## Installation

You can install the *eeganalysis* package from GitHub using:

```r
# Install devtools if needed
install.packages("devtools")

# Install eeganalysis
devtools::install_github("ChristosDalamarinis/eeganalysis")
```

---

## Quick Start

### Load the package

```r
library(eeganalysis)
```

### Read a BioSemi BDF file

```r
# Import EEG data from BioSemi file
eeg <- read_biosemi("~/data/subject_01.bdf")

# View summary information
print(eeg)

# Get detailed quality report
summarize_biosemi_import(eeg)
```

### Access EEG data

```r
# Access different components of the eeg object
eeg$data              # Raw EEG signal matrix
eeg$channels          # Channel names (e.g., "Cz", "Pz", "Oz")
eeg$sampling_rate     # Sampling rate in Hz
eeg$times             # Time vector in seconds
eeg$events            # Event/trigger markers
eeg$metadata          # Recording information
```

---

## Current Features

✅ **Data Import**
- Read BioSemi BDF files
- Automatic channel labeling
- Event/trigger extraction
- Metadata preservation

✅ **Data Organization**
- EEG object (S3 class) for structured data storage
- Preprocessing history tracking
- Quality assessment tools

⏳ **Planned Features**
- Filtering and artifact correction
- ICA (Independent Component Analysis)
- ERP and time-frequency analysis
- Statistical testing and visualization

---

## Requirements

### Dependencies

The package requires:
- **edfReader** - Read BioSemi BDF files
- **signal** - Signal processing (filtering, spectral analysis)
- **pracma** - Mathematical functions (interpolation, etc.)
- **ggplot2** - Visualization
- **tidyverse** - Data manipulation
- **tibble** - Modern data frames

### R Version
- R >= 3.5.0

---

## Notes

- **Current Focus:** BioSemi BDF file format
- **Future Formats:** EDF, EGI, EEGLAB, and other common formats will be added in future releases
- **Reference Scheme:** BioSemi CMS/DRL (Common Mode Sense / Driven Right Leg) reference is preserved

---

## Contributing

Contributions to the *eeganalysis* package are welcome! To contribute:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Make your changes with documentation
4. Submit a pull request

Please ensure your code:
- Follows the existing style
- Includes roxygen2 documentation comments (`#'`)
- Has unit tests (in `tests/testthat/`)
- Includes meaningful commit messages

---

## License

This package is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Citation

If you use *eeganalysis* in your research, please cite it as:

```
Dalamarinis, C. (2025). eeganalysis: An R Package for EEG Analysis.
GitHub: https://github.com/ChristosDalamarinis/eeganalysis
```

---

## Acknowledgments

This package builds on the following R packages:
- [edfReader](https://cran.r-project.org/web/packages/edfReader/) - EDF/BDF file reading
- [ggplot2](https://ggplot2.tidyverse.org/) - Visualization
- [tidyverse](https://www.tidyverse.org/) - Data science tools

