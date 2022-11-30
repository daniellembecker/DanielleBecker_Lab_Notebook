---
layout: post
title: Analyzing Acropora larval metabolism experiment respirometry data
date: '2022-11-29'
categories: Acropora_Laraval_Metabolism_Moorea_2022
tags: Acropora R Respirometry
---
This post details analysis of respirometry data for the *Acropora pulchra* larval metabolism experiment in Moorea in Oct 2022.  

[Notebook posts for this project can be found here](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/Moorea-2022-Coral-Spawning-and-Field-Expedition-Daily-Entry/).  

# Dataset overview  

We exposed *Acropora pulchra* larvae to 4 days of thermal stress in larval tank systems in Moorea, French Polynesia in October 2022. Larvae were expoesd to either high (+2-2.5°C) or ambient temperature. Within each temperature, tanks were either dosed with symbionts or not. However, larvae in the dosed tanks did not take up symbionts. Therefore, this is primarily a temperature experiment. 

Throughout the experiment, we sampled for gene expression, metabolomics, and respirometry on days 0, 2, and 4 of temperature exposure.  

Here, I am conducting the initial analysis of respirometry data.  

# Extracting respiration rates   

Oxygen evolution was measured first in the light (20 min at ~500 PAR). The purpose of this was to measure photosynthesis in anticipation of symbiont uptake. This light period will also standardize the light environment to account for any differences in light enhanced respiration.  

Oxygen consumption (respiration) was calculating using LoLinR using localized linear regressions. The script for rate extraction can be [found on GitHub](https://github.com/AHuffmyer/acropora_larval_metabolism/blob/main/scripts/Respirometry_Extraction.Rmd).  

For each sample, the calculations looked like this:  

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/respiration/respiration_rates_dark.png)  

This calculation shows the slope calculated of oxygen consumption during the dark (respiration) period.  

However, because no symbionts were taken up, photosynthesis and respiration rates were similar, indicating no oxygen production (described later below). I also calculated respiration rates using teh entire measurement period.  

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/respiration/respiration_rates_all.png)

In this example, you can see the start of the measurement period is highly variable. This is the light measurement period. Since there was no photosynthesis, respiration rates can be calculated across the entire measurement period. 

Rates were calculated as nmol oxygen per larva per minute and are on [GitHub here](https://github.com/AHuffmyer/acropora_larval_metabolism/blob/main/output/respirometry/oxygen_P_R_calc.csv).  

# Visualizing respiration rates  

First, I visualized respiration rates across run number colored by measurement temperature (either 27 or 30°C).  

![](https://raw.githubusercontent.com/AHuffmyer/acropora_larval_metabolism/main/figures/respirometry/resp_all_runs_boxplot.png)  

This plot raised several red flags right away. Runs 1-4 were taken on the first day of the experiment (day 0 of exposure), runs 5-12 were taken on day 2 of exposure, and run 13-16 were taken on the last day (day 4) of the experiment. Runs 5-8 and runs 9-12 were duplicate sets of runs of the same larval groups. Therefore, the values should be the same. There is a clear change in both the absolute values and the variability starting in run 9. 

The first thing I did was check the blank values. The blank values were also highly variable on these runs which further indicate a technical problem.  

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/respiration/blanks.png)  

There was no change in the software or equipment during this time and each day a new calibration was applied. It is still unclear why these runs were so different within this same day. So, I removed runs 9-12 from the dataset and then looked at the data.  

![](https://raw.githubusercontent.com/AHuffmyer/acropora_larval_metabolism/main/figures/respirometry/resp_filtered_runs_boxplot.png)  
 
In this figure, we can see that runs 1-4 and 5-8 are similar, which makes sense. Runs 13-16 are still different than the previous data in absolute value and variability. Due to the issues that showed up with runs 9-12, I am unsure if this is a biological difference or a technical artifact.

I then looked at the data at the treatment level first with all runs included:  

![](https://raw.githubusercontent.com/AHuffmyer/acropora_larval_metabolism/main/figures/respirometry/resp_all_runs_treatment.png) 
 
Here we can see values are higher and very variable on days 2 and 4. 

When I remove runs 9-12, the data looks like this:  

![](https://raw.githubusercontent.com/AHuffmyer/acropora_larval_metabolism/main/figures/respirometry/resp_filtered_runs_treatment.png) 

This suggests that the differences in respiration values are an artifact not a biolgical difference.  

This raises uncertainty in if we can trust data from runs >9.  

Photosynthesis rates show the same result as respiration, indicating no production of oxygen. Note here photosynthesis is shown as negative values (oxygen consumption). In the plots above respiration is shown as absolute values, so they are positive in the charts (but are negative due to oxygen consumption).  

With all runs photosynthesis looks like this:  

![](https://raw.githubusercontent.com/AHuffmyer/acropora_larval_metabolism/main/figures/respirometry/photo_all_runs_treatment.png) 

With runs 9-12 the data looks like this:  

With all runs photosynthesis looks like this:  

![](https://raw.githubusercontent.com/AHuffmyer/acropora_larval_metabolism/main/figures/respirometry/photo_filtered_runs_treatment.png)

In both photosynthesis and respiration we have the same issue with runs >9.  

# Next steps  

Next, I will try to dig deeper to find out why runs >9 were different and if we will be able to use the data.   

