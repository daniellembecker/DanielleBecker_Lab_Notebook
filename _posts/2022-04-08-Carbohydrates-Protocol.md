---
layout: post
title: Carbohydrates Protocol
date: '2022-04-08'
categories: Protocol
tags: Protocol Carbohydrates Physiology
---

This protocol details quantification of carbohydrate content in coral larvae in 96-well plates modified from the [Bove Baumann protocol](https://www.protocols.io/view/coral-carbohydrate-assay-for-96-well-plates-j8nlk4ro1g5r/v1) and the [Putnam Lab carbohydrate protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/20210713-Carbohydrate-Assay-Test/).  

# Materials  

## Reagents: 

Concentrated Sulphuric Acid (95% certified ACS grade)  
Phenol (certified ACS)  
MilliQ Water  
L-(-)-Glucose  

## Equipment:  

96-well plastic plates  
Water bath (room temperature)  
Vortex  
Fume hood  
Pipettes  
Plate reader (can read absorbance at 485 nm)  

# Reagent and Standard preparations  

### 10 mM Glucose stock solution (50mL)

Calculation to get to 10 mM L-(-)-Glucose:  

- Molecular Weight of L-(-)-Glucose: 180.16g  
- 1 M = 180.16g/1L Water  
- 1 mM (or 0.001 M) = 0.18g/L of water  
- 10 mM = 1.8g/1L of water  
- g of L-(-)-Glucose needed for 10mM in 50mL = 1.8g/200 = 0.009g (9mg)  

Add 50mL of milliQ water to a 50 mL falcon tube  
- label: 10 mM Glucose stock solution, Date, Initials  

Add 0.009g of L-(-)-Glucose  

Vortex and store in 4°C fridge.  

### 1 mM Glucose stock solution (50mL)

Add 45mL if milliQ water to a 50mL falcon tube  
- label: 1 mM Glucose stock solution, Date, Initials  

Add 5mL of 10 mM Glucose stock solution  

Vortex and store in 4°C fridge.  

# Protocol  

1. Homogenize and separate fractions of each sample [as described here](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/Montipora-capitata-2021-larval-temperature-experiment-sample-processing/).  
2. Take out samples and thaw at room temperature. Samples are run in triplicates, so you can do 22 samples per plate in addition to the standards. You can process samples for multiple plates at once and then read out the plates at the end of the protocol. I recommend doing no more than 2 plates at a time due to the number of samples and time to pipette.  
3. Label one tube per sample (e.g., one tube for host fraction, one tube for holobiont fraction, and one tube for symbiont fraction for each biological sample if desired).  
4. Label tubes for standard curve (see table below).  
5. While samples are thawing, create standards and blanks as shown below.   

| Tube ID | Concentration (mg/mL) | Vol water (uL) | Vol 1 mM Glucose (ul) | Vol 10mM Glucose (uL) |
|:-------:|:---------------------:|:--------------:|:---------------------:|:---------------------:|
|    B    |         0.0000        |      1000      |           0           |           0           |
|    1    |        0.00901        |       950      |           50          |           0           |
|    2    |        0.01802        |       900      |          100          |           0           |
|    3    |        0.02703        |       850      |          150          |           0           |
|    4    |        0.03604        |       800      |          200          |           0           |
|    5    |        0.05406        |       700      |          300          |           0           |
|    6    |         0.0901        |       500      |          500          |           0           |
|    7    |         0.1802        |        0       |          1000         |           0           |
|    8    |         0.3604        |       800      |           10          |          200          |
|    9    |         0.901         |       500      |           0           |          500          |

6. Vortex samples after thawing
7. Add 25 uL of homogenized larval tissue for the desired fraction (host, holobiont, or symbiont) and 975 uL milliQ water to pre-labelled test tube for all samples
  * **Note:** This was determined by testing a gradient of concentrations to generate concentrations within the standard curve. This may need to be altered depending on the sample type.  
8. Set up a room temperature water bath in the fume hood with test tube rack (i.e. DI water in a plastic bin)
  * This is to prevent the tubes from over heating from the addition of sulphuric acid
9. Add 25 uL of phenol to first sample
10. Vortex (in the hood) briefly  
11. Immediately add 2.5 mL sulphuric acid to the sample
12. Transfer to water bath
13. **Repeat steps addition of phenol and sulfuric acid for all tubes**
14. When the last sample is placed in the water bath, incubate all samples for 30 minutes (set a timer)    
15. Pipette 200 uL of all standards and samples into the bottom of the wells in a 96-well plate
  * Have a printed plate layout with well numbers and sample IDs prepared during the incubation step. Run all samples and standards in triplicates.  
16. After incubation, read on spectrophotometer at 485 nm and save data.  

# Calculations  

1. Create standard curve with known standard concentrations and absorbance values (y = mx + b)  
2. Using the resulting equation, convert sample absorbance to concentrations (mg/mL)  
3. Multiply sample concentration (mg/mL) by total slurry volume (mL) and dilution factor (1000/v of sample, usually 100 mL), then divide by surface area (cm2) for resulting units: mg/cm2  

R Scripts for these calculations can be found in my [GitHub repository](https://github.com/AHuffmyer).    

