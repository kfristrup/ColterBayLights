# ColterBayLights
R scripts and data for the paper: "National Park visitors perceive benefits for themselves and wildlife under blended red-white outdoor lighting."

This repository stores the code and data used to generate the results in the paper named above. These offer opportunities to replicate the results, and possibly explore other patterns in these data.

This R scripts use the R library _here_ to organize file access. Create a project directory containing three subdirectories: Data, Radiometry, Output, and Rcode (the last could be named anything). You must create an empty file ".here" in the top level (project) directory.
There are three R scripts:

GrteSocSciHeatmap4.R -- main script for processing

heatMkf6.R -- functions defined for the preceding file

ConsentChecks.r -- auxiliary analysis of differences between visitors who did or did not consent to participate in the survey.

06_AggBatDataWuv.csv,  GRTE lighting.xlsx, MooseWeatherHr.csv (yes, the Holy Grail reference was intentional), and MASTER_GRTE_NightSkies_Cleaned_NRB2.sav should be placed in the Data subdirectory.
AsphaltReflectance.RDat, LightMeasurementsLong.csv, and Vb2BrightnessCoefficients.csv should be placed in the Radiometry subdirectory.

Questions about the survey design and deployment should be directed to Prof. B. Derrick Taff at Penn State University -- bdt3@psu.edu. This repository is maintained by Kurt Fristrup, Kurt.Fristrup@colostate.edu
