---
layout: post
title: Mcapitata Development 16S Analysis - Mothur pipeline visualizations
date: '2022-12-13'
categories: Mcapitata_EarlyLifeHistory_2020
tags: 16S Mcapitata Molecular R
---
This post details recent updates to the *Montipora capitata* developmental time series 16S analysis. 

See these notebook posts for [the mothur bioinformatic pipeline that we used](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/16S-Analysis-in-Mothr-Part-1/) and [initial analyses in R](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/Analyzing-16S-Mothur-Results-in-R/).   

# Overview  

The dataset we are analyzing was collected across early life stages of *M. capitata* in Hawaii in 2020 including egg, embryo, larval, metamorphosed polyp, and attached recruit stages. See [my notebook page](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/categoryview/#mcapitata-earlylifehistory-2020) for more information on the project.    

Following bioinformatics steps above in the mothur pipeline, we lost several samples due to low sequence counts during subsampling. For this analysis, we are working with the following lifestages:  

- Embryos (38 hpf)    
- Embryos (65 hpf)     
- Larvae (93 hpf)  
- Larvae (163 hpf)  
- Larvae (231 hpf)  
- Metamorphosed Polyps (231 hpf)  
- Attached Recruits (183 hpf)  
- Attached Recruits (231 hpf)  

In this analysis, my goals were to vizualize:  

1. Alpha diversity (diversity metrics)   
2. Beta diversity (NMDS)  
3. Relative abundance of microbial taxa    

Working scripts for this analysis can be found on [GitHub here](https://github.com/AHuffmyer/EarlyLifeHistory_Energetics/blob/master/Mcap2020/Scripts/16S/16S_Mothur.Rmd).  

# 1. Alpha diversity  

Data for this project were analyzed using the `phyloseq` package in R (v1.42.0).  I first visualized alpha diveristy calculated using the Shannon and Inverse Simpson indices.  

![](https://raw.githubusercontent.com/AHuffmyer/EarlyLifeHistory_Energetics/master/Mcap2020/Figures/16S/alpha_rare.png)  

These data indicate that alpha diversity increases with development and is highest in the attached recruit stages. 

An ANOVA analysis shows significant variation in alpha diversity by lifestage. 

```
summary(aov(InvSimpson~lifestage, data=alpha_rare))

```

```
summary(aov(Shannon~lifestage, data=alpha_rare))

```

Next, I plotted beta diversity using a non-metric multidimensional scaling analysis (NMDS).  

# 2. Beta diversity/NMDS  

I used an NMDS to visualize beta diversity across lifestages.  

![](https://raw.githubusercontent.com/AHuffmyer/EarlyLifeHistory_Energetics/master/Mcap2020/Figures/16S/nmds_16S.jpeg)  

```
Square root transformation
```

Stress values are approximately 0.09-0.1 indicating a good fit.  

I also ran a PERMANOVA analysis to test for statistical variation across life stages.  

```
perm1<-adonis2(phyloseq::distance(physeq1_rare, method="bray") ~ lifestage, + data = samples)

```

The NMDS plot and PERMANOVA analysis show that there is significant variation across lifestages in 16S profiles. The variation over time is primarily driven across the x-axis of the NMDS plot. Further, there is higher dispersion in 16S microbial profiles in later life stages (larvae, metamorphosed polyps, attached recruits after 183 hpf) as compared to earlier life stages.  

Experimentally, this timepoint marks the start of settlement trials in different tank conditions than outdoors. Therefore, we have to consider that some of this variation could be due to different water conditions. We therefore will focus on comparisons within these two time periods (e.g., larvae vs metamorphosed polyps after 183 hpf).   

To investigate this further, I plotted relative abundance of taxa in our samples.  

# 3. Relative abundance   

I vizualized relative abundance (abundance / sum(abundance) * 100) using heatmaps in the package `complexHeatmap` in R.  

First, I calculated the mean relative abundance of each taxa for each lifestage group. I then filtered taxa to include only those that comprised >1% relative abundance within each lifestage. 

I first vizualized relative abundance at the level of phylum.  

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/16S/heatmap_phylum.png) 

We see that samples are dominated by Proteobacteria. It appears that Bacteriodes and unclassified Bacteria increase in later life stages.  

Then I vizualized relative abundance at the level of Genus:Family. Note that after filtering for taxa <1% relative abundance, we are left with 16 taxa.  

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/16S/heatmap_genus.png) 

From this dataset, we can specifically look at comparisons of lifestages within the pre-183 hpf and comparison of lifestages within the time period after 183 hpf.   

I noticed that in earlier lifestages there is a higher prevalence of Vibrio spp. There is also a peak of Oceanosprilliacaea (unclassified) in mid points of early life history. Alteromonas peaks in early stages of embryos. Pseudoalteromonas, Rhodobacteraceae, and Kordiimonadaceae are more prevalent in later time points.  

# Next Steps  

The next steps will be to look at the taxa present in samples more closely and consult with collaborators on this data set.     


 


