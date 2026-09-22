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
| Channel inspection & labeling | `inspect_biosemi_file()`, `identify_external_channels()`, `detect_external_channels()` | ✅ |
| Downsampling / Filtering | `downsample()`, `eeg_bandpass()`, `eeg_notch()` | ✅ |
| Bad-channel detection & repair | `find_bad_channels()`, `interpolate_bads()` | ✅ |
| Time-range annotations (bad stretches) | `annotate_amplitude()`, `annotate_muscle()`, `annotate_nan()`, `annotate_break()` | ✅ |
| Re-referencing | `eeg_rereference()` | ✅ |
| ICA artifact removal (fit, auto-detect, exclude, apply) | `fit_ica()`, `find_bads_eog()`, `find_bads_ecg()`, `find_bads_muscle()`, `apply_ica()` | ✅ |
| EOG regression (fast alternative to ICA for eye artifacts; continuous and epoched data) | `fit_eog_regression()`, `apply_eog_regression()`, `subtract_evoked()` | ✅ |
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
inspect_biosemi_file("path/to/your/file.bdf")
inspect_triggers(eeg_data)
```

### Downsample and filter

``` r
eeg_data <- downsample(eeg_data, target_rate = 256)
eeg_data <- eeg_bandpass(eeg_data, l_freq = 0.1, h_freq = 40)
eeg_data <- eeg_notch(eeg_data, freqs = 50)
```

### Find and repair bad channels

``` r
eeg_data <- set_montage(eeg_data, create_montage())   # scalp positions, needed for repair
eeg_data <- find_bad_channels(eeg_data)
eeg_data <- interpolate_bads(eeg_data)
```

### Re-reference to average

``` r
eeg_data <- eeg_rereference(eeg_data, ref = "average")
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
epochs <- epoch_eeg(eeg_data, events = c(1, 2), tmin = -0.2, tmax = 0.8)
```

## Contributing

Contributions to the *eeganalysis* package are welcome! If you would like to contribute, please fork the repository and submit a pull request with your changes. Please ensure that your code follows the existing style and includes appropriate documentation and tests. Note: currently focusing on Biosemi files, other formats will be added in future releases.

## Folder structure

``` text
eeganalysis
│
├── R/                                       ← Package source code
│   ├── eeg_class.R                          ← Core EEG data structure
│   │   ├── new_eeg()                        ← Create an eeg object
│   │   └── print.eeg()                      ← Display eeg object nicely
│   │
│   ├── read_bdf_native.R                    ← BioSemi file import
│   │   └── read_bdf_native()                ← Import .bdf files
│   │
│   ├── extract_bdf_events.R                 ← Trigger/event extraction
│   │   ├── extract_bdf_events()             ← Parse trigger codes from status channel
│   │   ├── summary_bdf_events()             ← Summarize extracted events
│   │   └── validate_bdf_events()            ← Validate and report on events
│   │
│   ├── label_bdf_events.R                   ← Event labeling
│   │   ├── label_bdf_events()               ← Attach human-readable labels to trigger codes
│   │   └── apply_trigger_labels()           ← Apply a label scheme to events
│   │
│   ├── channel_info2.R                      ← Electrode database & channel inspection
│   │   ├── get_electrode_database()         ← Access 64-ch BioSemi electrode database
│   │   ├── get_electrode_position()         ← Get coordinates for a specific electrode
│   │   ├── scan_biosemi_channels()          ← Quick channel-name scan of a .bdf file
│   │   ├── inspect_biosemi_file()           ← Preview a .bdf file without full import
│   │   ├── detect_electrode_naming_system() ← Identify naming convention (10-20/10-10/BioSemi)
│   │   ├── plot_electrode_3d()              ← 3D electrode visualization (Cartesian)
│   │   └── plot_electrode_3d_spherical()    ← 3D electrode visualization (spherical)
│   │
│   ├── setexchannels.R                      ← External channel management
│   │   ├── identify_external_channels()     ← Interactive labeling (EOG, EMG, ECG, GSR)
│   │   ├── detect_external_channels()       ← Automated external channel detection
│   │   └── apply_external_labels()          ← Apply user-defined labels to data
│   │
│   ├── downsample.R                         ← Smart downsampling
│   │   └── downsample()                     ← Downsample with anti-aliasing filter
│   │
│   ├── filter1.R                            ← FIR filtering
│   │   ├── eeg_bandpass()                   ← Hamming-window FIR bandpass filter
│   │   └── eeg_notch()                      ← Multi-band notch filter for line noise
│   │
│   ├── bad_channels.R                       ← Bad-channel detection
│   │   └── find_bad_channels()              ← Flat/amplitude/outlier + spatial + LOF checks
│   │
│   ├── interpolate.R                        ← Bad-channel repair
│   │   └── interpolate_bads()               ← Spherical-spline interpolation (MNE-matched)
│   │
│   ├── annotations.R                        ← Time-range bad-data marking
│   │   ├── annotate_amplitude()             ← Flag flat/spiking stretches (or bad channels)
│   │   ├── annotate_muscle()                ← Flag EMG bursts (high-freq envelope z-score)
│   │   ├── annotate_nan()                   ← Flag amplifier dropouts (NA runs), per channel
│   │   └── annotate_break()                 ← Flag dead time between experimental blocks
│   │
│   ├── rereference.R                        ← Re-referencing utilities
│   │   └── eeg_rereference()                ← Change reference scheme (average/custom)
│   │
│   ├── ica1.R                               ← ICA artifact removal
│   │   ├── new_ica()                        ← Create an ICA container
│   │   ├── fit_ica()                        ← Decompose data (PCA + FastICA)
│   │   ├── get_sources()                    ← Get component time courses
│   │   ├── plot_ica_sources()               ← Plot component time courses
│   │   ├── get_component_topography()       ← Get a component's scalp topography
│   │   ├── plot_ica_topography()            ← Plot a component's scalp topography
│   │   ├── ica_component_summary()          ← Summary table of components
│   │   ├── set_exclude()                    ← Mark components for removal
│   │   ├── apply_ica()                      ← Reconstruct data without excluded components
│   │   ├── plot_ica_overlay()               ← Before/after overlay for a channel
│   │   └── print.eeg_ica()                  ← Display ICA object nicely
│   │
│   ├── ica_detect.R                         ← ICA Phase 2: automatic bad-component detection
│   │   ├── find_bads_eog()                  ← Flag eye-movement components (EOG correlation)
│   │   ├── find_bads_ecg()                  ← Flag heartbeat components (ECG correlation or CTPS)
│   │   ├── find_bads_muscle()               ← Flag muscle components (spectral slope + topography)
│   │   └── corrmap()                        ← Match a template topography across subjects' ICAs
│   │
│   ├── regression.R                         ← EOG regression for eye artifacts (continuous + epoched)
│   │   ├── new_eog_regression()             ← Create/validate an EOG regression model
│   │   ├── subtract_evoked()                ← Hide each trial's evoked response before fitting
│   │   ├── fit_eog_regression()             ← Learn per-channel weights from the EOG channel(s)
│   │   ├── apply_eog_regression()           ← Subtract weight × EOG from each EEG channel
│   │   └── print.eeg_eog_regression()       ← Display regression model nicely
│   │
│   ├── montage.R                            ← Electrode montage handling
│   │   ├── new_montage()                    ← Create a montage object
│   │   ├── create_montage()                 ← Build a montage from a template
│   │   ├── set_montage()                    ← Attach a montage to an eeg object
│   │   └── print.montage()                  ← Display montage object nicely
│   │
│   ├── topography.R                         ← Scalp topography plotting
│   │   └── plot_topography()                ← Interpolated scalp-map plot
│   │
│   ├── epoch2.R                             ← Epoching functions
│   │   ├── inspect_triggers()               ← Inspect event triggers before epoching
│   │   ├── epoch_eeg()                      ← Extract time-locked epochs around events
│   │   └── plot_epochs()                    ← Visualize extracted epochs
│   │
│   ├── erp_analysis.R                       ← ERP averaging (in development)
│   │   └── average_epochs()                 ← Average epochs into an ERP
│   │
│   ├── fourier.R                            ← Spectral analysis
│   │   ├── eeg_fft()                        ← One-sided FFT spectrum
│   │   ├── eeg_psd_welch()                  ← Welch's method PSD
│   │   ├── eeg_multitaper()                 ← Multitaper (DPSS) PSD
│   │   ├── eeg_band_power()                 ← Band-limited power
│   │   └── print.eeg_spectrum()             ← Display spectrum object nicely
│   │
│   ├── plot_signal.R                        ← Raw signal plotting
│   │   └── plot_eeg_signal()                ← Plot raw/continuous EEG traces
│   │
│   ├── eeg_summary.R                        ← Dataset diagnostics
│   │   └── eeg_summary()                    ← Print a diagnostic summary of an eeg object
│   │
│   └── imports.R                            ← Centralized @importFrom declarations
│
├── man/                                     ← Auto-generated help files (130, one per exported/internal function)
│   ├── new_eeg.Rd
│   ├── read_bdf_native.Rd
│   ├── find_bad_channels.Rd
│   ├── interpolate_bads.Rd
│   ├── fit_ica.Rd
│   ├── apply_ica.Rd
│   └── ...Rd                                ← Remaining help files
│
├── tests/testthat/                          ← Unit tests (testthat, one file per module, 20 total)
│   ├── test-read_bdf_native.R
│   ├── test-ica1.R
│   ├── test-bad_channels.R
│   ├── test-interpolate.R
│   └── ...R                                 ← Remaining test files
│
├── NAMESPACE                                ← Exported functions (auto-generated)
├── DESCRIPTION                              ← Package metadata
├── LICENSE                                  ← License information
├── README.md                                ← Package overview
├── .gitignore                               ← Git ignore rules
├── .Rbuildignore                            ← Build ignore rules
└── eeganalysis.Rproj                        ← RStudio project file
```

Note: this structure is expanded and updated as new modules are added.
