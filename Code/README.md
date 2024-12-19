This sub-directory contains the R-files that are needed to reproduce the simulation study and the real data analysis as reported in the preprint 'Scan statistics for the detection of anomalies in M-dependent random fields with applications to image data' by C. Kirch, P. Klein and M. Meyer available on arXiv:2311.09961. 

The file `functions.R` contains all necessary functions, while the simulation study can be reproduced using `simulation_study.R`. The real data analysis can be reproduced using `real_data.R`.

The thresholds for the simulation study and the real data analysis can be uploaded from `Threshold_simulation_study.RData` and `Thresholds_real_data.RData`, respectively.  The calculations to generate these thresholds for both the simulation study and the data analysis can be found in `Generate_Thresholds.R`.

The color palettes used for the coloring of the images can be uploaded from `Palettes.RData`.

All calculations were done under R version 3.4.4.
