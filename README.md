# Repository for time-consistent innovation policy with creative destruction
### Jonathan Goldberg, David Lopez-Salido, Nicholas von Turkovich
### 2023

## Overview
Included are the necessary codes and data files to produce the results in the paper. The main blocks of code can be generally divided into three categories: **Model**, **Calibration**, and **Output**. The model code contains the core scripts needed to solve the model and perform the exercises that generate the data underlying the figures and tables in the paper. The calibration codes are used to solve for the baseline parameters of the model given a set of target moments. The output codes are those that actually generate the tables and figures in the paper. For a detailed description of the organization of the codebase and how to replicate the results, see sections Organization and Replication Instructions below. The codes are built primarily in Matlab though a few scripts are written in R (see Dependencies below).

For ease of reproducing our figures and tables, the data files generated by the code are contained in a zip folder. Many of the exercises are resource intensive and can take several hours even with numerous additional processors and additional memory. Considerable effort has been made to reduce runtime including efficient coding practices and parallelization wherever possible. Still, the vast majority of the exercises in the paper require access to a computing cluster to be completed in a reasonable amount of time. Please see Replication Instructions for more details.

## Organization

The root folder for the project is **innov_commit**. Within this folder, when the project is first initiated, there will be a single folder (**Code**) and a few core scripts (*main.m*, *slurm_job_matlab.sh*, *slurm_job_short.sh*). The script *main.m*, when run the first time, will instantiate the remaining folder structure needed for the project including the folder **Data** in the root directory as well as subfolders within **Data** and **Output**. The folders, **Code** contains all the scripts for calibrating the model, running experiments, and producing output. Note that the user should not need to interact with these files directly as *main.m* aggregates all the actions one would need to carry out. The folder **Data** will contain staging folders for any raw data that is downloaded (**Data/Raw**) as well as an intermediate folder that primarily contains the *.mat* files output from the exercises managed by *main.m.* The folder **Output** has four subfolders (**Figures**, **Logs**, **Tables**, **Tracked**). **Output/Figures** will contain the figures of the paper generated from the code. Likewise, **Output/Logs** and **Output/Tables** will contain logs generated by the code as well as the tables in the paper. 

The folder **Output/Tracked** is of particular importance as this contains static data files, including shell *.tex* files for the table output, the compiled output from the two calibration exercises, and the data files produced by the routines that are necessary to generate all the figures and tables. As mentioned before, these files are included as it is time/resource-intensive to run some of the exercises and this allows the user to quickly replicate the figures/tables.

However, should one want to run the exercises from scratch, that is also handled by the *main.m* script.

## Replication Instructions

### About *main.m*

This file contains all of the configurations needed to replicate the results in the paper. The user should need only to interact with the code through this script or the SLURM wrappers (*slurm_job_matlab.sh*, *slurm_job_short.sh*). 

Immediately, the user should adjust the variable *nlopt_path* which specifies where the code can find the NLOPT routines. This isn't strictly necessary unless running the calibration code. 

The script is divided into three sections:
- Setup (instantiating folders, setting up filepaths, detecting number of worker threads, etc.)
- Switches (this is where the user makes choices about what to run, notes indicate runtime expectations)
- Code execution (these blocks of code translate the switches and run the code needed to either produce the output or data associated with the switch)

### Running using the static files in **Output/Tracked**

To save time, the user can produce the charts/figures by following these steps.
- Run main.m once without any switches, this will instantiate the file structure
- Copy *Output/Tracked/Data_Intermediate.zip* into **Data/Intermediate** and unzip
- Copy *Output/Tracked/Output_Logs_Main.zip* into **Output/Logs/Main** and unzip
- Set *switch_plots = true;* and all others to false
- Run either *main.m* or use *slurm_job_short.sh* to generate a SLURM job (some parameters in this file would need to be adjusted depending on your own cluster)

The output will be found in **Output/Figures/Main** and **Output/Tables**.

### Running exercises from scratch

For all intents and purposes, how one interacts with the code is the same as in the prior section. By turning other switches to true, the code is executed in the same way. The difference being is that it will write new files to **Data/Intermediate** (over the static files if they were unzipped). 

One should take care to setup the SLURM script to accommodate the demands of the exercise. Generally, requesting 50 CPUs and ~25 GB per CPU has been sufficient. The most intensive exercise is the *switch_triple* which can take a considerable amount of time and resources given that it nests multiple dependent optimization routines.

## Data

Very little in the way of raw data is required by the code. For the calibration exercises, to re-run the moment targets, it accesses static files from **Output/Tracked**.

## Dependencies
- NLOPT (Matlab, v 2.6.2 though v 2.7.1 appears to differ very little)
- R packages:
	- tidyverse
	- readxl
	- quantmod
	- zoo
	- purrr
	- lubridate
	- here

## Contact
For any questions related to the codebase, please reach out to Nick von Turkovich (nvonturk@mit.edu).
