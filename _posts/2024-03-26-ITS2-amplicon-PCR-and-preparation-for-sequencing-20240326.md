---
layout: post
title: ITS2 amplicon PCR and preparation for sequencing 20240326
date: '2024-03-26'
categories: Larval_Symbiont_TPC_2023
tags: PCR ITS2 Molecular Mcapitata
---

# Overview 

# Protocol 

Base protocol
List of modifications 

# Equipment and supplies 

xxx

# 20240325

## Sample dilution 

All gDNA has previously been extracted and quantified by Jill detailed in these posts: [xxx](xxx), [xxx](xxx), [xxx](xxx).

Today we diluted gDNA samples in preparation for ITS2 amplification and preparation for sequencing. Dilution is performed to standardize the amount of input in each reaction. Samples were diluted to a concentration of 4ng/uL with the exception of two samples, which had concentrations lower than 4ng/uL in the original gDNA. For these samples, we diluted to 3ng/uL and we will increase the input volume to have a total of 4ng gDNA input for every reaction.  

We diluted to a total volume of 10uL. For several samples, we were unable to obtain a Qubit concentration reading (see previous posts linked above), so for these we based dilutions off of the Nanodrop concentration.  

Dilutions were generated in individually labeled 0.6mL tubes and not in 96-well plates as described in Emma's base protocol. 

We will use Ultrapure H20 as our negative control. Our positive control is an Mcap2020 project sample D18, which was [successfully amplified and sequenced previously]().  

| Sample | Group       | Treatment | Qubit Conc (ng/uL) | Nanodrop Conc (ng/uL) | DNA for dilution (uL) | Water | Total volume | Final Concentration (ng/ul) | Notes                                                                    |
|--------|-------------|-----------|--------------------|-----------------------|-----------------------|-------|--------------|-----------------------------|--------------------------------------------------------------------------|
| R55    | Cladocopium | Ambient   |               14.3 | NA                    |                   2.8 |   7.2 |           10 |                           4 |                                                                          |
| R56    | Cladocopium | Ambient   | NA                 |                     6 |                   6.7 |   3.3 |           10 |                           4 | Off Nanodrop concentration                                               |
| R57    | Cladocopium | Ambient   |              14.45 |                  15.8 |                   2.8 |   7.2 |           10 |                           4 |                                                                          |
| R58    | Cladocopium | Ambient   |              15.05 |                   8.4 |                   2.7 |   7.3 |           10 |                           4 |                                                                          |
| R59    | Cladocopium | Ambient   |               6.66 |                  10.7 |                   6.0 |   4.0 |           10 |                           4 |                                                                          |
| R60    | Cladocopium | Ambient   | NA                 |                   6.3 |                   6.3 |   3.7 |           10 |                           4 | Off Nanodrop concentration                                               |
| R61    | Mixed       | Ambient   |               8.89 |                   8.4 |                   4.5 |   5.5 |           10 |                           4 |                                                                          |
| R62    | Mixed       | Ambient   |              13.15 |                  14.2 |                   3.0 |   7.0 |           10 |                           4 |                                                                          |
| R63    | Mixed       | Ambient   |              12.65 |                   6.3 |                   3.2 |   6.8 |           10 |                           4 |                                                                          |
| R64    | Mixed       | Ambient   |               8.69 |                     9 |                   4.6 |   5.4 |           10 |                           4 |                                                                          |
| R65    | Mixed       | Ambient   |               8.24 |                   8.3 |                   4.9 |   5.1 |           10 |                           4 |                                                                          |
| R66    | Mixed       | Ambient   |              20.05 |                  17.3 |                   2.0 |   8.0 |           10 |                           4 |                                                                          |
| R67    | Wildtype    | Ambient   |              11.65 |                  13.6 |                   3.4 |   6.6 |           10 |                           4 |                                                                          |
| R68    | Wildtype    | Ambient   |               4.34 |                   9.7 |                   9.2 |   0.8 |           10 |                           4 |                                                                          |
| R69    | Wildtype    | Ambient   | NA                 |                     3 |                  10.0 |   0.0 |           10 |                           3 | Off Nanodrop concentration/max out at 10 uL; reduce water in PCR; 3ng/ul |
| R70    | Wildtype    | Ambient   |               14.3 |                     3 |                   2.8 |   7.2 |           10 |                           4 |                                                                          |
| R71    | Wildtype    | Ambient   |               14.1 |                    12 |                   2.8 |   7.2 |           10 |                           4 |                                                                          |
| R72    | Wildtype    | Ambient   |               14.1 |                  13.2 |                   2.8 |   7.2 |           10 |                           4 |                                                                          |
| POS    | Mcap2020    | Ambient   |               24.2 | NA                    |                   1.7 |   8.3 |           10 |                           4 | MCap 2020 (below)                                                        |
| M60    | JA          | Ambient   |               3.56 | NA                    |                   8.4 |   1.6 |           10 |                      2.9904 | Reduce water in PCR, 3ng/uL                                              |
| M72    | JA          | Ambient   |               4.62 | NA                    |                   8.7 |   1.3 |           10 |                           4 |                                                                          |
| M80    | JA          | Ambient   |               10.8 | NA                    |                   3.7 |   6.3 |           10 |                           4 |                                                                          |

After dilutions were made, samples were stored at -20°C until the next step.  

# 20240328

Today we will conduct the PCR and gel QC steps of this protocol.  








## Sample and master mix calculations  

**Number of reactions with sufficient concentration (4ng/uL)** 

| Reactions:           |    |
|----------------------|----|
| Samples (AH)         | 17 |
| Samples (JA)         |  2 |
| Samples (POC)        |  0 |
| Positive (Mcap2020)  |  1 |
| Negative (1 uL H20)  |  1 |
| TOTAL SAMPLES        | 21 |
| TOTAL REACTIONS      | 63 |
|                      |    |
| Added for 10% error  | 69 |


**Master mix for samples with sufficient concentration (4ng/ul)**  

| Component            | Per Rxn (uL) | Rxns | Total Volume (uL) |
|----------------------|--------------|------|-------------------|
| 2X Phusion Mastermix |         12.5 |   69 |             866.3 |
| F primer (10uM)      |          0.5 |   69 |              34.7 |
| R primer (10uM)      |          0.5 |   69 |              34.7 |
| Ultra Pure H20       |         10.5 |   69 |             727.7 |
| DNA (4ng)            |            1 |      |                   |
| TOTAL RXN VOLUME     |           25 |      |                   |


**Number of reactions with low concentration (4ng/uL)** 

| Reactions:           |   |
|----------------------|---|
| Samples (AH)         | 1 |
| Samples (JA)         | 1 |
| Samples (POC)        | 0 |
| Positive             | 0 |
| Negative             | 0 |
| TOTAL SAMPLES        | 2 |
| TOTAL REACTIONS      | 6 |
|                      |   |
| Added for 10% error  | 7 |

**Master mix for samples with low concentration (3ng/ul)**  

| Component            | Per Rxn (uL) | Rxns | Total Volume (uL) |
|----------------------|--------------|------|-------------------|
| 2X Phusion Mastermix |         12.5 |    7 |              82.5 |
| F primer (10uM)      |          0.5 |    7 |               3.3 |
| R primer (10uM)      |          0.5 |    7 |               3.3 |
| Ultra Pure H20       |        10.17 |    7 |            67.122 |
| DNA (3ng)            |         1.33 |      |                   |
| TOTAL RXN VOLUME     |           25 |      |                   |

**Primer dilutions** 

| Primer Dilution | Desired Volume | Current Concentration (uM) | Desired Concentration (uM) | Stock Primer Vol | H20  | Total reactions |
|-----------------|----------------|----------------------------|----------------------------|------------------|------|-----------------|
| F primer        |           38.0 |                        100 |                         10 |             3.80 | 34.2 |              76 |
| R primer        |           38.0 |                        100 |                         10 |             3.80 | 34.2 |              76 |

## PCR amplification 


xxxx


## Gel for QC 

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/uri_2024/20240328_ITS2_gel_1.png?raw=true)

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/uri_2024/20240328_ITS2_gel_1.png?raw=true)


## PCR amplification to QC for contamination using new water source  

# 20240329

## Gel for QC of yesterday's PCR 