---
layout: post
title: Revising E5 Time Series Multivariate Plasticity Analysis
date: '2023-03-16'
categories: E5 Analysis
tags: Physiology R Multivariate
---
This post details revisions to calculations for multivariate physiological plasticity analysis for the [E5 Rules of Life project](https://urol-e5.github.io/). 

All code used in this post can be found in the [E5 GitHub repository](https://github.com/urol-e5/timeseries/blob/master/time_series_analysis/5_plasticity_analysis.Rmd).  

The original analysis and description of calculations can be found [in this post](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/E5-Multivariate-Physiological-Plasticity-Analysis/).  

# Overview 

The E5 physiological time series project includes data on coral and symbiont physiology in colonies of *Acropora*, *Pocillopora*, and *Porites* corals in Moorea, French Polynesia. Samples of individual colonies were collected at four timepoints in the year of 2020: Jan, Mar, Sep, and Nov.  

For the physiological responses collected (e.g., biomass, symbiont density, calcification, protein, antioxidant capacity), we have analyzed to investigate the influence of species, site, and time on physiology. We have found that physiological response of corals is influenced by environmental variability both across space (site) and seasons (time point).  

This can be visualized using trajectory principle component analyses, which visually show the movement in physiological trajectories of each species across space and time. 

![trajectory pca](https://raw.githubusercontent.com/urol-e5/timeseries/master/time_series_analysis/Figures/Multivariate/Full_Panel_PCAs_All.png)  

The next step in our analysis is to quantify the "movement" or physiological acclimatization in corals between species and sites. To do this, we conducted an analysis to quantify the differences in physiological profiles for each colony across time in multivariate space. [This post detials our original analysis to do this and quanify physiological plasticity](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/E5-Multivariate-Physiological-Plasticity-Analysis/).  

# Revisions to Analysis  

Previously, we utilized Euclidean distance in multivariate physiological profiles for each colony across time. I calculated the Euclidean distance between the PCA coordinates (using all PC's) for each colony between Jan - Mar, Jan - Sept, and Jan - Nov. 

Because these calculations were conducted from the change from the first time point to subsequent time points, we are missing some of the information on the total distance traveled (i.e., dispersion). To do this, I am going to calculate the euclidean distances for each colony as: Jan to Mar + Mar to Sept + Sept to Nov = total dispersion. This will provide a more accurate calculation of physiological plasticity across the time series. 

The previous analysis produced the following figure, which we will compare to at the end of this analysis.  

![plasticity](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/E5_physiology/Plasticity_Figure.png)   

# Revising Calculations  

See the [previous notebook post](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/E5-Multivariate-Physiological-Plasticity-Analysis/) for code for previous calculations as described above. In this analysis, I am using the same approach, but changing the equation for calculating the distances between time points for each colony. 

All scripts in this post can be found [in the E5 GitHub here](https://github.com/urol-e5/timeseries/blob/master/time_series_analysis/5_plasticity_analysis.Rmd).  

First, I generated a dataframe that calculated the euclidian distance between all timepoints for each colony in a pairwise comparison.  

```
data_all<-master%>%
```

The data looked like this: 

```
                 row                col    value
```

Second, I selected only the comparisons that we want to calculate (Jan-March, March-Sept, and Sept-Nov) for each colony and removed the duplicate reciprocal comparisons (Jan-March distance and March-Jan distance, which are the same). 

I then generated rules to calculate distances. This was complicated because we were often missing data from 1-3 time points for some colonies. This meant that I needed to make rules for calculating distances. For example, if a colony was missing data from March, I wanted to calculate the distance from January to September instead. Because we are looking at the sum of all distances, I only kept colonies that we had at least 3 time points of information. This resulted in 36 colonies out of 123 without plasticity calculations.   

The code for these rules is below:  

```
mutate(sum_1=rowSums(.[ , c("T1_T2", "T2_T3", "T3_T4")]))%>% #add all values for corals with all four timepoints recorded by summing T1-T2, T2-T3 and T3-T4 
```

The full code for calculations and data manipulation are below here:  

```
#calculate distances between subsequent time points (T1-T2, T2-T3, T3-T4) and sum them together

```

This produced data that looked like this:  

```
  Colony1 colony_dispersion site_code   species 
```

Where `colony_dispersion` is the metric of interest to calculate plasticity.  

This calculation process was repeated for Holobiont (all metrics), Host, and Symbiont metrics. 

# Plotting 

I generated plots for these plasticity values at the Holobiont, Host, and Symbiont levels and displayed anova analyses on the effects of site and species on these values.  

![](https://raw.githubusercontent.com/urol-e5/timeseries/master/time_series_analysis/Figures/Plasticity/Plasticity_Figure.png)

This figure shows multivariate dispersion, our metric of plasticity, on the y-axis across site and species. Each point is an individual colony with the group mean displayed in the dark lines and points. There are some interesting take aways from this analysis. 

- Holobiont physiological plasticity was significantly different by species. Pocillopora had the lowest plasticity levels across site with Acropora and Porites more equal. There was a trend for reduced plasticity in Acropora at Hilton, but the site x species interaction was not significant. 
- Host physiology was different by species as expected with lowest values in Pocillopora. The site x species interaction was significant. This is driven by higher plasticity in Acropora than Porites at Mahana, but lower in Acropora at Hilton. Plasiticty for each species was dependent on site. 
- For both holobiont and host physiology, site was not significant, but there was a non-significant trend for variation by site. 
- Symbiont physiology, interestingly, was not different by species or site. This suggests that holobiont plasticity was driven more by the host than the symbiont. 

This interpretation differs from the previous version, where host and symbiont plasticity were significantly influenced by site and species.   

