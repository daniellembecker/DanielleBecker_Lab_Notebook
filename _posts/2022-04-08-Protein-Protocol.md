---
layout: post
title: Protein Protocol
date: '2022-04-08'
categories: Protocol
tags: Protocol Protein Physiology
---

This protocol details quantification of soluble protein from coral larval samples. This protocol is modified from the [E5 Protein Protocol](https://github.com/urol-e5/protocols/blob/master/2020-01-01-Total-Protein-Protocol.md).  

# Soluble Protein Protocol

Quantification of soluble protein is described here. For soluble + insoluble protein, see the full protocol linked above.  

Contents  
- [**Materials**](#Materials)    
- [**Protocol**](#Protocol)  
- [**Standard Table**](#Table)  
- [**References**](#References)  
 
1. <a name="Materials"></a> **Materials**
    - 	[Pierce BCA Protein Assay Kit from Thermo Scientific](https://www.thermofisher.com/order/catalog/product/23225?SID=srch-srp-23225).  
    - 	Plastic 96 Well plate
    - 	Incubator (37°C), preferably with rocker  
    - 	Plate reader Spectrophotometer
    -  Pipettes P10, P200, P1000 and tips
    -  1.5ml microfuge tubes
    -  DI water

2. <a name="Protocol"></a> **Protocol** 

**Sample Preparation for Soluble Protein**  

Separate fractions of sample (host, symbiont, holobiont fractions) using the separation protocols [described here](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/Montipora-capitata-2021-larval-temperature-experiment-sample-processing/).  

After separations are complete, thaw desired sample fractions.  

1. Thaw tissue homogenate.  
2. Vortex to re-suspend the sample. 
3. <a name="Table"></a> **Standard Table**  
**Preparation of Diluted Albumin (BSA) Standards**    
*Prepare these standards while samples are thawing*  

4. Use the following table as a guide to prepare a set of protein standards. For this protocol we will use the microplate procedure. Diluent is DI water Type II. Each vial will be a sterile 1.5 mL microcentrifuge tube. Label the cap of the microcentrifuge tube with the Vial ID ("A", "B", etc.).  

| Vial | Volume of Diluent (μL) | Volume of Source of BSA (μL) | Final BSA Concentration (μg/mL) |
|------|------------------------|------------------------------|---------------------------------|
| A    | 0                      | 300 of Stock                 | 2000                            |
| B    | 125                    | 375 of Stock                 | 1500                            |
| C    | 325                    | 325 of Stock                 | 1000                            |
| D    | 175                    | 175 of vial B dilution       | 750                             |
| E    | 325                    | 325 of vial C dilution       | 500                             |
| F    | 325                    | 325 of vial E dilution       | 250                             |
| G    | 325                    | 325 of vial F dilution       | 125                             |
| H    | 400                    | 100 of vial G dilution       | 25                              |
| I    | 400                    | 0 (Blank)                    | 0                               |

**Preparation of the BCA Working Reagent (WR)**   

1. Use the following formula to determine the total volume of WR required depending on the number of samples:      
(# standards + # samples) x (# replicates) x (volume of WR per sample) = total volume WR required
  
For this protocol, we will use 9 standards and 200 μL of WR is required for each sample in the microplate procedure. Samples and standards will be run in triplicate.    

> *(9 standards + # samples) x (3 replicates) x (200 μL of WR) = total volume WR required*  

2. Prepare WR by mixing 50 parts of BCA Reagent A with 1 part of BCA Reagent B (50:1, Reagent A:B) in a clean protein-free container of the appropriate size, based on how many samples are going to be run.  

> *Total volume WR required / 51 = volume of Reagent B needed.*  
> *Volume of Reagent B * 50 = volume of Reagent A needed.*  

Using these volumes, create the working reagent. It is recommended to round up the amount needed (~1-2 mL extra).  

**Microplate Procedure (Sample to WR ratio = 1:8) from Pierce BCA Protein Assay Kit:**  

*Note that in this protocol for larval tissue I am using diluted sample by adding DI water to each well. It is recommended to run dilution gradients of each sample to determine the optimal sample concentration for measurements within the standard curve.  

1. Pipette 25 μL of each standard and 10 μL of each sample into triplicate microplate wells. Cover plate after adding standards and samples.  
2. Using a multi-channel pipette, add 15 μL of DI water to each well with samples. Note that the desired total concentration of each sample and standard is 25 μL.    
2. Add 200 μL of the working reagent (WR) to each well and mix by gently pipetting up and down in the well.  
3. Cover the plate and incubate at 37&deg;C for 30 minutes.  
4. Measure the plates at 562 nm absorbance.  
4. Subtract the average 562 nm absorbance measurement of the Blank standard replicates from the 562 nm measurements of all other individual standard and unknown sample replicates.  
5. Calculate the standard curve by plotting the average Blank-corrected 562nm measurement for each BSA standard vs. its concentration in μg/mL. Use the standard curve equation to determine the protein concentration of each unknown sample.  

4. <a name="References"></a> **References**  
[Pierce BCA Protein Assay](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/MAN0011430_Pierce_BCA_Protein_Asy_UG.pdf)
