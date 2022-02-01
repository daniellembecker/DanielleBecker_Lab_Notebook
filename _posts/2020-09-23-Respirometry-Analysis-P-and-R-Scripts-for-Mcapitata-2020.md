---
layout: post
title: Respirometry Analysis: P and R Scripts for Mcapitata 2020 
date: '2020-08-30' 
categories: Mcapitata_EarlyLifeHistory_2020 Scripts_Analysis 
tags: Analysis Respirometry 
--- 

This entry references scripts and analysis for respirometry data obtained from Oxy-10 SDR with microplate. 

# Overview 

Scripts were modified from K. Wong and S. Gurr scripts in the LoLinR package in R.  

All scripts and data can be found in my EarlyLifeHistory_Energetics Github repository here: https://github.com/AHuffmyer/EarlyLifeHistory_Energetics.

To analyze and obtain slopes for photosynthesis and respiration in the same run, I modified previous respiration scripts with the following changes:

[Slope Extraction Script](https://github.com/AHuffmyer/EarlyLifeHistory_Energetics/blob/master/Mcap2020/Scripts/Respirometry_Extraction_Mcap2020.Rmd)
* Split data by start and end times specified in dataframe to separate P and R runs.  
* Calculate blanks for each section of P and R runs and subtract from slopes obtained.  
* Output corrected slopes into .csv file.  

[Analysis and Plotting Script](https://github.com/AHuffmyer/EarlyLifeHistory_Energetics/blob/master/Mcap2020/Scripts/Respirometry_Plotting_Stats_Mcap2020.Rmd)
* Plotting R and P values by net and gross metrics for each life stage measured in this experiment.  
* Application of statistics (parametric and non-parametric).  
