---
layout: post
title: Correcting subsampling error in Mcap2020 ITS2 code and revise results
date: '2024-03-24'
categories: Mcapitata_EarlyLifeHistory_2020
tags: ITS2 R
---

This post describes correcting code and revising analyses for *Montipora capitata* 2020 developmental time series project.  

# Overview 

While conductings some tests for subsampling thresholds, I found an error in my ITS2 code. The error was that sequence abundance was being transformed to relative abundance prior to subsampling, and therefore it was not subsampling at the desired sequence count threshold. I have corrected this error and have summarised the revised results below.  

Prior to this revision, we saw no change in symbiont communities across development. We can now see significant changes across development, but these changes are very minor and the major ITS2 profiles show only minimal shifts between eggs and later life stages.   

# DIV-level analysis 







# Profile-level analysis 