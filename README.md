# Electron_Fake_Rate
## Introduction
This repository contains scripts and plots for analyzing the electron-faking photon rate in the monophoton final state. The analysis is based on the 2017 datasets (41.53 fb<sup>-1</sup> luminosity and a center-of-mass energy of 13 TeV) collected by the Compact Muon Solenoid (CMS) Detector at the Large Hadron Collider (LHC).

The electron-faking photon ratio factor is utilized to quantify the rate of electron misidentification. The investigation focuses on different probe transverse momenta (p<sub>T</sub>) ranges and pseudorapidity (η) regions.

The scripts and corresponding plots are implemented using C++, ROOT, and CMS-Software (CMSSW). These tools enable comprehensive analysis and visualization of the electron-faking photon rate.

## Contents
**scripts/:** Contains the C++ scripts used for data processing and analysis. 


**Plots/:** Includes the generated plots and visualizations.

    Plots_eta: Include the generated plots for different pseudorapidity(η) regions

    Plots_pT: Include the generated plots for different transverse momenta
  
    Plots_fakerate: Include the generated plots for electron-faking photon rate for different pT and eta(η) values


**README.md:** Provides an overview of the repository and instructions for running the analysis.

## Usage
To reproduce the analysis or modify the scripts, follow these steps:

1. Clone the repository:              git clone https://github.com/SambitSatapathy9/Electron_Fake_Rate.git
2. Navigate to the cloned repository: cd Electron_Fake_Rate/
3. Execute the analysis script:       ./scripts/electron.C
4. Explore the generated plots in the Plots/.

## Dependencies
To run the analysis, you need the following dependencies installed:

1. C++ compiler
2. ROOT data analysis framework
3. CMS-Software (CMSSW) environment

Make sure to set up these dependencies properly before running the analysis.

## Contributing
Contributions to this repository are welcome. If you encounter any issues or have suggestions for improvements, please create an issue or submit a pull request.

  
