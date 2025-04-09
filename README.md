# StructuredActionPrepVMDM
data and analysis/modeling code for all studies in "Mental graphs structure the storage and retrieval of visuomotor associations"

These scripts use R (4.2.1), RStudio (2022.07.1), and MATLAB (9.12.0, 2022a, Update 4). Links to installation help are below. Installation time can vary by computer.

We recommend using the precompiled data to run analyses in order to expedite run time of the scripts. 

NHB_hierarchy_Expt123.Rmd - (R version 4.2.1,https://www.r-project.org/ , RStudio version 2022.07.1, https://rstudio-education.github.io/hopr/starting.html, rmarkdown) main analysis script for Experiments 1-4. Script loads in csv files of usable data for each study (saved in compiledData) and conducts main analyses and figures. Run time is <5min. Experiment 1 = sonaHL, Experiment 2 = oneDFL, Experiment 3 = threeDFL, Modified structure control = fingerHL.

NHB_hierarchy_Expt45.Rmd - Main analysis script for Experiments 5 and 6. Same as previous except that it loads in different data and additionally has forced response analyses. Run time is ~5min. Experiment 4 = hierDyn, Experiment 5 = flatDyn.

NHB_SurveyAnalyses.Rmd - Script to plot survey results. Raw data included in strucureIntuitions.csv. Run time < 1min.

NHB_hierarchy_additionalAnalyses - Script for additional analyses added in revision. Run time <5 min. This script includes analyses for Experiments 1-5.

NHB_Expt6_multiday_learn - Main analysis script for the multiday experiment. Loads in compiled data for analyses. Analyses include all subjects that participated in the longitudinal experiment (not just participants with the forced response phase too.) Run time <5min.

NHB_Expt6_multiday_forced - Main analysis script for the forced response multi-day sessions. Also loads in data from Experiments 4 and 5 (other forced response experiments). Run time <5min.

FOLDERS
compiledData - csv files of experimental data after exclusions. Each file is labeled with the experiment number. Data from baseline RT task is in motorData_XXX.csv files. Additionally, data from Experiments 5 and 6 are tagged with whether it is data from the learning task or the forced response task.

modeling - This folder has scripts for the response time model for Experiments 5 and 6 in Matlab.
  .mat files are saved model fits
  fit_prepTime_model_final.m - (Matlab, version 9.12.0, 2022a, Update 4, https://www.mathworks.com/help/install/ug/install-products-with-internet-connection.html, licensed through Yale University) fits any of the models and saves output
  modelingForcedRT_simulations.m - Loads in model fit and plots simulations of response time functions. Also, saves    output of simulations for use in plotting scripts in NHB_hierarchy_Expt56.Rmd
  func_***.m - function used for each model in fit_prepTime_model_final.m
  presponseSims - folder where simulation csvs are saved (also saved in the StructuredActionPrepVMDM/simulations)

rawData - data files for each subject, organized by experiment.
  Experiment*_compiled.csv - Same as in StructuredActionPrepVMDM/compiledData
  compiledData_wExclusions.Rmd - Pulls in subjectwise data, prepares extra columns for analysis, performs exclusions, and motor correction. This script can take >5min per block of code.
