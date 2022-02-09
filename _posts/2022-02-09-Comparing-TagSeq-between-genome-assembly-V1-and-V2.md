---
layout: post
title: Comparing TagSeq between genome assembly V1 and V2
date: '2022-02-09'
categories: Mcapitata_EarlyLifeHistory_2020 Analysis
tags: Mcapitata Molecular R GeneExpression
---

This post details analysis of TagSeq data from previous alignment to *M. capitata* [version 1 assembly](http://cyanophora.rutgers.edu/montipora/) and [version 2 assembly](http://cyanophora.rutgers.edu/montipora/). 

# Comparing TagSeq results between *M. capitata* genome assembly versions 1 and 2 

Recently, a new version of the *M. capitata* genome was made available from Rutgers University at http://cyanophora.rutgers.edu/montipora/. To see how this new version impacts our alignment with our TagSeq data, I am running the [bioinformatics pipeline used previously](https://github.com/AHuffmyer/EarlyLifeHistory_Energetics/blob/master/Mcap2020/Scripts/TagSeq/TagSeq_BioInf.md) with the new version. The purpose of this analysis is to see whether the clear lifestage comparisons found in previous analysis are still present when working with the new assembly.  

Previously, our analysis shows clear lifestage separation in our gene expression data obtained through TagSeq sequencing and WGCNA analysis.  

![genes pca](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/GeneExpression/pca_genes.png)

![genes wgcna](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/GeneExpression/wgcna_genes.png)