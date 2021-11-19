---
layout: post
title: E5 Time Series Data QC
date: '2021-11-17'
categories: Analysis E5
tags: R
---
This posts describes the QC process from November 2021 for the E5 Rules of Life timeseries data. 

# Overview

The purpose of this round of QC is to go through each time point and check that data is complete for each colony as expected and to account for any missing data by correcting data entry errors or tracking down reasons why data is missing.  

I use the `completeness.Rmd` scripts in the [E5 repository](https://github.com/urol-e5/timeseries). These scripts generate tables that populate each colony with `TRUE` or `FALSE` statements for each data type. From this, I can track those that have `FALSE` statements and account for the missing data.  

The most common reasons for `FALSE` statements are: 

- Colonies were not sampled or were not found in the field  
- Colony tag ID was incorrectly entered in either the data sets or in the metadata file  
- Sample was not able to be processed due to spilling or errors in processing  
- Typos in the data frame for colony ID 

All of these issues can be traced through looking at the `completeness.Rmd` outputs as well as the notebook photographs.

My notes and tracking of completeness as well as "answers" to missing data can be seen in [this excel document in our Git repo](https://github.com/urol-e5/timeseries/blob/master/metadata/E5_QC.xlsx). 

### Timepoint 1  

All missing data has been resolved for timepoint 1. There were 19 colonies listed as having missing data. These were due to either typos in data entry (corrected) or colonies that were not collected/found/measured at this time point.  

In some cases, colony data was not found in individual data sets. This was confirmed by checking that the colony was not analyzed according to the notebook. 

### Timepoint 2  

All missing data has been resolved for timepoint 2. This timepoint was the most complete with only 7 colonies listed with missing data. Corrected a typo for incorrect labeling fo colony names. Other missing data was from colonies that were not collected/found.  

### Timepoint 3  

One factor contributing to missing data in this timepoint is the longer period of time from time point 2, so several colonies were not found or had to be relabeled. In total, 47 colonies had missing data.  

There were several instances of POR or POC colonies being labeled as the opposite. Corrected these in data entry.  

I also noticed we have very low n for all metrics for site 3 Acropora. This is due to mortality and inability to sample. Fewer were found in TP3 than in TP4, some colonies were later recovered.  

### Timepoint 4  

This timepoint had several major problems:   
- Colonies that were not sampled, largely Acropora  
- Typos and errors in data entry. There were several instances of incorrect data entry in chlorophyll, molecular, surface area, and tac. These were corrected.  
- In some cases, notes were made about samples being put in whirlpacks with other colony labels. I tracked down these colonies and corrected to the correct colony id.  
- Incomplete notebooks; some assays that were conducted after time in Moorea are not recorded in these notebooks.  
- The protein and TAC platemaps were switched, leading to lots of errors in protein/tac.  
- There were several colonies for each metric that were lost or not processed for some of the metrics.  

There were 33 colonies that had entry errors.  

