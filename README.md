# StructuredActionPrepVMDM
data and analysis/modeling code for all studies in "Structured action preparation during visuomotor decision-making"

NHB_hierarchy_Expt1234.Rmd - main analysis script for Experiments 1-4. Script loads in csv files of usable data for each study (saved in compiledData) and conducts main analyses and figures.

NHB_hierarchy_Expt56.Rmd - main analysis script for Experiments 5 and 6. Same as previous except that it loads in different data and additionally has forced response analyses.

compiledData - csv files of experimental data after exclusions. Each file is labeled with the experiment number. Data from baseline RT task is in motorData_XXX.csv files. Additionally, data from Experiments 5 and 6 are tagged with whether it is data from the learning task or the forced response task.

modeling - This folder has scripts for the response time model for Experiments 5 and 6.
  .mat files are saved model fits
  fit_prepTime_model_final.m - fits any of the models and saves output
  modelingForcedRT_simulations.m - Loads in model fit and plots simulations of response time functions. Also, saves    output of simulations for use in plotting scripts in NHB_hierarchy_Expt56.Rmd
  func_***.m - function used for each model in fit_prepTime_model_final.m
  presponseSims - folder where simulation csvs are saved (also saved in the StructuredActionPrepVMDM/simulations)

rawData - data files for each subject, organized by experiment.
  Experiment*_compiled.csv - Same as in StructuredActionPrepVMDM/compiledData
  compiledData_wExclusions.Rmd - Pulls in subjectwise data, prepares extra columns for analysis, performs exclusions, and motor correction.
