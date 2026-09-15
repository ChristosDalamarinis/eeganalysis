# eeganalysis

<img src="man/figures/logo.png" align="right" height="139"/>

*eeganalysis* gives EEG researchers a straightforward R-native workflow — from raw BioSemi recordings to epoched, analysis-ready data — without leaving the R ecosystem.

Author: Christos Dalamarinis

Contact: [dalamarinischristos\@gmail.com](mailto:dalamarinischristos@gmail.com)

<!-- badges: start -->

[![R-CMD-check](https://github.com/ChristosDalamarinis/eeganalysis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ChristosDalamarinis/eeganalysis/actions/workflows/R-CMD-check.yaml) [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental) ![Version](https://img.shields.io/badge/version-0.0.0.9000-blue) ![License: MIT](https://img.shields.io/badge/License-MIT-green.svg) ![Last Commit](https://img.shields.io/github/last-commit/ChristosDalamarinis/eeganalysis)

<!-- badges: end -->

## Installation

You can install the *eeganalysis* package from GitHub using the following command:

``` r
# install.packages("devtools") # Uncomment if devtools is not installed
devtools::install_github("dalamarinischristos/eeganalysis")
```

## Status

Version 0.0.0.9000 — under active development. Currently reads BioSemi `.bdf` files; additional formats (`.edf`, `.set`) are planned.

## Features

| Stage | Function(s) | Status |
|------------------------|------------------------|------------------------|
| Import BioSemi `.bdf` | `read_bdf_native()` | ✅ |
| Channel inspection & labeling | `inspect_bdf_channels()`, `identify_external_channels()`, `detect_external_channels()` | ✅ |
| Downsampling / Filtering | `downsample()`, `eeg_bandpass()`, `eeg_notch()` | ✅ |
| Bad-channel detection & repair | `find_bad_channels()`, `interpolate_bads()` | ✅ |
| Re-referencing | `eeg_rereference()` | ✅ |
| ICA artifact removal | `fit_ica()`, `plot_ica_sources()`, `apply_ica()` | ✅ |
| Epoching & visualization | `epoch_eeg()`, `plot_epochs()` | ✅ |
| Spectral analysis | `eeg_fft()`, `eeg_psd_welch()`, `eeg_band_power()` | ✅ |
| Topography & montage | `plot_topography()`, `create_montage()` | ✅ |
| ERP averaging | `average_epochs()` | ⚠️ in development |

## Dependencies

- **signal** — filtering and downsampling
- **dbscan** — Local Outlier Factor bad-channel detection
- **akima** — spherical-spline interpolation for bad-channel repair and topography maps
- **fastICA** — ICA decomposition for artifact removal
- **MASS** — linear algebra support for ICA and interpolation
- **plotly** — interactive 3D electrode and topography plots
- **dplyr** — data manipulation
- **ggplot2** — plotting

## Quick Start

``` r
library(eeganalysis)
```

### Import and inspect

``` r
eeg_data <- read_bdf_native("path/to/your/file.bdf")
inspect_bdf_channels("path/to/your/file.bdf")
inspect_triggers(eeg_data)
```

### Downsample and filter

``` r
eeg_data <- downsample(eeg_data, target_srate = 256)
eeg_data <- eeg_bandpass(eeg_data, l_freq = 0.1, h_freq = 40)
eeg_data <- eeg_notch(eeg_data, freqs = 50)
```

### Find and repair bad channels

``` r
eeg_data <- find_bad_channels(eeg_data)
eeg_data <- interpolate_bads(eeg_data)
```

### Re-reference to average

``` r
eeg_data <- eeg_rereference(eeg_data, ref_type = "average")
```

### Remove artifacts with ICA

``` r
ica <- new_ica(n_components = 20)
ica <- fit_ica(ica, eeg_data)
plot_ica_sources(ica, eeg_data)                 # inspect components
ica <- set_exclude(ica, components = c(1, 3))   # mark eye/muscle artifacts
eeg_data <- apply_ica(ica, eeg_data)
```

### Extract epochs around events

``` r
epochs <- epoch_eeg(eeg_data, events = c(1, 2), time_window = c(-0.2, 0.8))
```

## Contributing

Contributions to the *eeganalysis* package are welcome! If you would like to contribute, please fork the repository and submit a pull request with your changes. Please ensure that your code follows the existing style and includes appropriate documentation and tests. Note: currently focusing on Biosemi files, other formats will be added in future releases.

## Folder structure

``` r
eeganalysis
├── R/            ← Package source code
├── man/          ← Auto-generated help files
├── tests/        ← testthat unit tests
├── data/         ← Example datasets
├── NAMESPACE     ← Exported functions (auto-generated)
├── DESCRIPTION   ← Package metadata
├── LICENSE       ← License information
└── README.md     ← Package overview
```
