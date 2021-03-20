# Supplementary-material-for-Plasma-Turbulence-at-Comet-67P-Rosetta-Observations
This repository contains supplementary material to Suranga Ruhunusiri, "Plasma Turbulence at Comet 67P: Rosetta Observations", submitted to Journal of Geophysical Research: Space Physics, 2020.

The repository contains the MATLAB programs for reading Rosetta magnetometer data, perform various computations mentioned in the manuscript, and generate the figures shown:

1. A MATLAB program for reading Rosetta magnetometer data: RPCMAG_data_reader.m [3 kB]

2. A MATLAB program for computing magnetic field power spectra and spectral parameters (spectral-break frequency and spectral indices): RPCMAG_data_analysis.m [12 kB]

3. A MATLAB program for plotting monthly medians and uncertainties for spectral-break frequency, low- and high-frequency spectral indices:    spec_break_and_indices_plotter.m [7 kB] 

3. A MATLAB program for computing and plotting the monthly occurrence rates and their uncertainties for turbulent processes TP1-TP8: occurrence_rate_plotter_for_turbulent_processes.m [16 kB]

4. A MATLAB program for computing and plotting spatial occurrence rates for turbulent processes TP1-TP4: spatial_map_plotter_for_turbulent_process_occurrence_rate.m [13 kB]

5. A MATLAB program for computing autocorrelation function for magnetic field: magnetic_field_autocorrelation_function_calculator.m [5 kB]

6. A MATLAB program for computing and plotting monthly median correlation times and their uncertainties: median_correlation_value_plotter.m [4 kB] 

7. MATLAB data files containing the following data for Rosetta observations at comet 67P from September 2014 to September 2016 (total file size=41 MB):
    column 1: low-frequency spectral index
    column 2: uncertainty for low-frequency spectral index
    column 3: high-frequency spectral index 
    column 4: uncertainty for high-frequency spectral index
    column 5: spectral-break frequency (Hz)
    column 6: data validity indicator 
    column 7: mean X-CSEQ location for the spacecraft
    column 8: mean Y-CSEQ location for the spacecraft
    column 9: mean Z-CSEQ location for the spacecraft
