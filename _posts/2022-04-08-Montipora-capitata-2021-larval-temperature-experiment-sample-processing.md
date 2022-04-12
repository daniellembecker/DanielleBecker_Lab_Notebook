---
layout: post
title: Montipora capitata 2021 larval temperature experiment sample processing
date: '2022-04-08'
categories: Mcapitata_LarvalTemp_2021
tags: Mcapitata Processing Protocol Metabolomics Physiology Carbohydrates Protein StableIsotopes
---
This post details sample processing completed in April 2022 for the 2021 *M. capitata* larval thermal exposure experiment conducted in Hawaii.  

# Sample Processing Overview  

We processed samples in the Putnam Lab at the University of Rhode Island collected from June 2021 thermal exposure of *M. capitata* larvae in Hawaii. During this round of processing, we conducted metabolomic extractions and physiological assays including protein and carbohydrate measurements. Additional assays will be conducted by members of the Putnam Lab. 

Putnam Lab PhD students and candidates Jill Ashey, Kevin Wong, and Emma Strand collaborated in this sample processing. 

Protocols are linked in the sections below. This post contains notes and information on optimizing processes in April 2022. 

Notebook posts detailing the experimental design for these samples can be found [here](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/categoryview/#mcapitata-larvaltemp-2021).    

# Separations of host, holobiont, and symbiont fractions for physiological assays  

Samples were collected in Hawaii 2021 for a variety of assays in separate tubes. In order to obtain physiological data from the same replicate, we are processing samples with the highest biomass ("E" tubes). The goal was to process protein, carbohydrate, chlorophyll, and symbiont cell density data from each replicate to allow for normalization and correlation across metrics.  

These larval samples were collected after 3 days of temperature exposure during this experiment. Samples processed include two replicate samples from each larval conical (n=6 in ambient and n=6 in high) for physiological metrics. Physiology samples therefore included a total of 24 samples. 

In order to process samples, we first homogenize larvae and separate into host, holobiont, and symbiont fractions, tracking volumes of each fraction for calculation and normalization. 

These samples are larvae in tubes with water removed. This generates a concentrated larval tissue sample for use in physiological assays.  
 
**Separation protocol was conducted as follows:**   
1. Thaw sample on ice  
2. Add 600uL of DI water to sample  
3. Homogenize using immersion tip homogenizer for 30-45 sec  
4. Transfer 500uL of homogenate into a new tube for the *holobiont* fraction ("ID-Holo")  
5. Transfer the remaining homogenate (200-300uL) to a new tube labeled for the *symbiont* fraction ("ID-Sym"). Track the volume.  
6. Centrifuge the sample at 2000g for 90 sec to pellet the symbiont.  
7. Remove the supernatant for the *host* fraction and move to a new tube ("ID-Host"). Record the volume.  
8. There are now 3 tubes per sample for each fraction.  
9. Wash homogenizer with 10% bleach, 70% ethanol, and DI water in between samples.  

Using this protocol, we noticed that the symbiont pellet contained a layer of lipid-type material due to the high lipid content of the larvae. Due to this poor separation, we will only analyze the *host* and *holobiont* fractions. This will allow us to obtain data on the host and provide comparisons to the holobiont fraction, which can provide information on the symbiont contribution to each response.  

# Carbohydrate Assays  

We processed a host and holobiont fraction sample for each of the 24 samples (N=48) for total carbohydrates. 

The protocol was modified from [C. Bove protocol HERE](). 

Briefly, 25 uL of respective fractioned samples were added to 975 uL of DI water in 5mL tubes. We also generated a gradient of standards including a blank using L-(-)-Glucose as a standard. 44uL of phenol and 2.5mL of sulfuric acid was then added to each tube and allowed to incubate for 30 minutes. Following the incubation, 200uL of samples and standards were then added to a 96-well plate in triplicates. Absorbance was then read on a plate reader at 485nm.  

**Initial Plate (6 April 2022)**:  
We first ran a plate with samples to optimize sample volume. First, we added 100uL of test samples ("B" tubes) to 900uL of DI water at the initial step of the protocol, which is the volume used for adult coral tissue. We also ran a set of test samples at 25uL of sample and 975uL of DI water as a test for dilution. We then read the plate. We found that the diluted samples were within the range of the standard curve while the 100uL samples were far greater than the standards. Moving forward we selected 25uL as our volume of sample for both host and holobiont fractions. 

**Plates 1-3 (6 April 2022)**:  
We then ran 3 plates for Mcap2021 samples with a standard curve run on each plate. Following our [modified protocol](). We also included n=4 samples from Jill Ashey's *Astrangia* samples.  

Assays were completed with all values within the standard curve. Preliminary data analysis shows that holobiont fraction carbohydrate content was higher than host, as expected.  

These data will be normalized to soluble protein for analysis.  

![INSERT CARB PICTURES HERE]()

Data and scripts for carbohydrates can be found [on the GitHub repository here]().  

# Protein Assays  

We measured soluble protein content in the host and holobiont fractions of all samples using the Pierce BCA protein assay kit (Thermo Fischer) following The [Putnam Lab and E5 Coral protein protocol here]().  

As conducted for the carbohydrate assay above, we modified the protocol for the sample volume used to test whether dilutions are necessary. The protocol specifies that 25uL of sample is used for adult coral tissues. We tested 25uL of sample as well as diluted 10uL sample with 15uL of DI water.  

Briefly, we generated a standard curve using bovine albumin as specified in the kit protocol. We then added the standards and specified concentration of larval sample (host and holobiont fractions) in a 96-well plate in triplicates. We then added 200uL of the working reagent to each well using a micro-channel pipette. The working reagent was a 50:1 solution of Reagent A:Reagent B using the kit protocols. Plates were then incubated at 37°C on a rocking table for 30 minutes prior to reading absorbance at 562nm.  

**Initial Plate (7 April 2022)**:  
We ran test samples on an initial plate to check whether dilution of samples was required. Samples that were diluted (10uL sample with 15uL DI water) were within the range of the standard curve, while samples that were not diluted (25uL of sample) were outside the range of the curve. Sample volume was modified for the protocol with no other modifications made.  

**Plates 1-3 (7 April 2022)**:  
We then ran 3 plates for Mcap2021 samples with a standard curve run on each plate. Following our [modified protocol](). We also included samples from Jill Ashey's *Astrangia* samples.  

These data will be used to normalize other physiological assays.  

![INSERT PROTEIN PICTURES HERE]()

Data and scripts for protein can be found [on the GitHub repository here]().  

# Stable Isotope Metabolomic Extractions  

Stable isotope metabolomic samples processed here included n=10 samples collected in *M. capitata* larval stable isotope incubations over a 24 hr time series. These samples will allow us to better understand the timing of isotopic label incorporation. We also processed n=8 samples collected from incubations of *P. acuta* larvae from Moorea in December 2021. These samples will include host, symbiont, and holobiont fractions to examine our capacity to detect the isotopic label in each fraction and ensure our separations and extractions are working correctly. 

We followed a modified protocol from Kevin Wong's [metabolomic extraction protocol here](). Briefly, we separated fractions as specified above, ensuring that all materials were washed with bleach, ethanol, and water between each use to avoid isotope contamination. Samples were then added to extraction buffer in a glass dounce, with samples kept on ice throughout the entire protocol. After extraction, samples were centrifuged and the supernatant was removed, buffered with ammonium bicarbonate, and transferred to autosampler vials and stored at -80°C.  

These vials will be sent to Rutgers Metabolomics Shared Resource for analysis of isotopic labeling in metabolomic profiles at positive and negative polarities.  

Following results of these test samples, Mcap2021 and Moorea Pocillopora samples will be further processed by Kevin Wong and Jill Ashey.   

![INSERT METABOLOMICS PICTURES HERE]()

# Additional Processing 

Physiology samples will be further processed for (1) chlorophyll, (2) symbiont cell densities, and (3) larval size led by Emma Strand.  

