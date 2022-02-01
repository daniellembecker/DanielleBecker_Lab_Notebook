---
layout: post
title: Analyzing 16S Mothur Results in R
date: '2022-01-31'
categories: Analysis Mcapitata_EarlyLifeHistory_2020
tags: 16S Mcapitata Molecular Protocol R
---

This post details analyzing relative abundance and diversity metrics in R output from mothur analysis of 16S data for the 2020 early life history *Montipora capitata* time series.  

# Prepare Data  

We previously transferred out the following files from the [mothur pipline described here](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/16S-Analysis-in-Mothr-Part-1/).  

This includes the following, 

Bray-Curtis matrices (lt.dist file is not subsampled/rarefied; lt.ave.dist file is subsampled/rarefied):     

```
mcap.opti_mcc.braycurtis.0.03.lt.dist
mcap.opti_mcc.braycurtis.0.03.lt.ave.dist
```

Taxonomy file:  

```
mcap.taxonomy
```

Shared files (subsample.shared is the file with rarefied; .shared does not have subsampling):  

```
mcap.opti_mcc.0.03.subsample.shared
mcap.opti_mcc.shared
```

From these files we will be comparing results from rarefying and therefore removing eggs/embryo samples with low sequences, or not conducting rarefication and keeping all samples.  

I opened these files in Excel and resaved as .txt tab delimited files. 


