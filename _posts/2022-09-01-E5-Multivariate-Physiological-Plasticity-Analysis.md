---
layout: post
title: E5 Multivariate Physiological Plasticity Analysis
date: '2022-09-01'
categories: E5 Analysis
tags: Physiology R Multivariate
---

This post details multivariate physiological plasticity analysis for the [E5 Rules of Life project](https://urol-e5.github.io/).  

# Overview  

The E5 physiological time series project includes data on coral and symbiont physiology in colonies of **Acropora**, **Pocillopora**, and **Porites** corals in Moorea, French Polynesia. Samples of individual colonies were collected at four timepoints in the year of 2020: Jan, Mar, Sep, and Nov.  

For the physiological responses collected (e.g., biomass, symbiont density, calcification, protein, antioxidant capacity), we have analyzed to investigate the influence of species, site, and time on physiology. We have found that physiological response of corals is influenced by environmental variability both across space (site) and seasons (time point).  

This can be visualized using trajectory principle component analyses, which visually show the movement in physiological trajectories of each species across space and time.  

First, we viewed physiological responses in all responses measured:  

![trajectory PCAs](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/E5_physiology/Full_Panel_PCAs.png)

We also viewed responses separated by holobiont (e.g., calcification, host protein) and symbiont responses (e.g., symbiont density, chlorophyll):  

![trajectory PCAs holobiont](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/E5_physiology/Full_Panel_PCAs_Holobiont.png)

![trajectory PCAs symbiont](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/E5_physiology/Full_Panel_PCAs_Symbiont.png)

This analysis shows us that physiology varies across species, space, and time and we have additionally explored univariate drivers of these shifts. 

The next step in our analysis is to quantify the "movement" or physiological acclimatization in corals between species and sites. To do this, we need to conduct an analysis to quantify the differences in physiological profiles for each colony across time in multivariate space. This post detials our analysis to do this and quanify physiological plasticity.  

# Analysis  

This analysis utilizes Euclidean distance in multivariate physiological profiles for each colony across time. In brief, we calculate the Euclidean distance between the PCA coordinates (using all PC's) for each colony between Jan - Mar, Jan - Sept, and Jan - Nov.  These PCA coordinates are generated from the PCA analyses shown above. 

The scripts to generate these calculations is as follows:  



This process is repeated for physiological profiles for 1) all responses, 2) holobiont responses, and 3) symbiont responses. 

# Results  

# Next Steps  

 

