---
layout: post
title: E5 POC and POR Species ID PCR and RFLPs
date: '2024-03-26'
categories: E5
tags: PCR RFLP Molecular Sanger Protocol
---

This post details protocls and results from identifying species of *Porites* and *Pocillopora* samples for the E5 RoL project. 

Danielle Becker  
Ariana Huffmyer  

# Overview 

Previously, our lab has conducted Sanger sequencing of genetic markers to identify species within *Porites* and *Pocillopora* samples from the E5 RoL project using genetic markers.  

*Pocillopora* species are identified using Sanger sequencing of the mtORF marker and a gel based RFLP of the POC Histone marker. *Porites* species are identified using Sanger sequencing of the H2 marker.  

Protocols were developed by Zoe and Hollie based on previous work from [Johnston et al. 2018](https://peerj.com/articles/4355/) and [Tisthammer et al. 2020](https://peerj.com/articles/8550/).  

Zoe and Hollie have previously processed samples using these protocols:  
   
- Zoe's PCR protocol [is here](https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/PCR-Protocol/) and Species ID protocol [is here](https://github.com/zdellaert/ZD_Putnam_Lab_Notebook/blob/master/protocols/SpeciesID-via-PCR-Sanger-Sequencing.md). 
- Zoe's previous notebook posts using this protocol are [here](https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/E5-PCR-POC-POR-Sites2and3-SpeciesID/) and [here](https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/E5-PCR-POC-POR-Site1-SpeciesID/)
- Hollie's previous notebook posts using this protocol are [here](https://github.com/hputnam/Putnam_Lab_Notebook/blob/master/_posts/2023-07-23-PocID.md) and [protocol here](https://github.com/hputnam/Putnam_Lab_Notebook/blob/master/_posts/Species%20ID%20via%20PCR%20and%20Sanger%20Sequencing.md). 

After the processing detailed in Zoe's notebooks, we need to run an RFLP assay on the *Pocillopora* samples to distinguish between *P. meandrina* and *P. grandis* in the samples that were distinguished from *P. verrucosa* through [Sanger sequencing of the mtORF marker](https://github.com/hputnam/Putnam_Lab_Notebook/blob/master/_posts/Species%20ID%20via%20PCR%20and%20Sanger%20Sequencing.md).  

These samples are:  

| Sample  | Marker      |
|---------|-------------|
| POC-200 | POC Histone |
| POC-205 | POC Histone |
| POC-207 | POC Histone |
| POC-215 | POC Histone |
| POC-217 | POC Histone |
| POC-238 | POC Histone |
| POC-239 | POC Histone |
| POC-248 | POC Histone |
| POC-254 | POC Histone |
| POC-366 | POC Histone |
| POC-386 | POC Histone |
| POC-395 | POC Histone |
| POC-41  | POC Histone |
| POC-45  | POC Histone |
| POC-55  | POC Histone |

We also need to re run two samples of *Porites* that failed sequencing and/or amplification to distinguish between *P. lutea/lobata* and *P. evermanni*.  

These samples are:  

| Sample  | Marker                   |
|---------|--------------------------|
| POR-362 | Reamplify and resequence for H2 |
| POR-83  | Reamplify and resequence for H2 |

We will conduct PCR and RFLP for the *Pocillopora* samples and will conduct PCR and Sanger sequencing for the *Porites* samples.   

# Protocol

We are following Zoe's Species ID PCR protocol [here](https://github.com/zdellaert/ZD_Putnam_Lab_Notebook/blob/master/protocols/SpeciesID-via-PCR-Sanger-Sequencing.md).  

We will note below where we deviate from the standard protocol.  

# Equipment and supplies 

- Sample Preservation and Lysis [Zymo Research DNA/RNA Shield R1100-250](https://github.com/hputnam/Putnam_Lab_Notebook/blob/master/images/Zymo_r1100-250_dna_rna_shield.pdf)
- DNA Extraction [Zymo Research Quick-DNA™ Miniprep Plus Kit D4069](https://github.com/hputnam/Putnam_Lab_Notebook/blob/master/images/d4068_d4069_quick-dna_miniprep_plus_kit.pdf)    
- Foward primer [FatP6.1 200µM Stock IDT)](https://github.com/hputnam/Putnam_Lab_Notebook/blob/master/images/Fatp6.1_IDT_Spec_328104852.pdf) 
- Reverse primer [RORF 200µM Stock IDT](https://github.com/hputnam/Putnam_Lab_Notebook/blob/master/images/RORF_IDT_Spec_328104853.pdf)         
- Master Mix [EmeraldAmp GT PCR Master Mix](https://github.com/hputnam/Putnam_Lab_Notebook/blob/master/images/TaKaRa_Emerald_RR320A_DS.pdf)
- EmeraldAmp GT PCR Master Mix is a loading-dye-added version of EmeraldAmp MAX PCR Master Mix that is optimized for great performance and convenience in both standard and high-throughput PCR applications.
- Loading Dye [NEB 6X Purple Loading Dye NEB Cat # B7024S]()        
- Gel Stain [Biotium GelGreen Nucleic Acid Gel Stain, 10,000X in Water Fisher Cat NC9728313](https://www.fishersci.com/shop/products/gel-green-stain-5ml/NC9728313#?keyword=NC9728313)
- DNA Ladder [1kb Gel ladder](https://github.com/hputnam/Putnam_Lab_Notebook/blob/master/images/NEB_1kb_Ladder_N3232S.png?raw=true) 
- Gel ladder [GeneRuler 100 bp DNA Ladder (Thermo FIsher Catalog SM0241)](https://github.com/hputnam/Putnam_Lab_Notebook/blob/master/images/SM0241_GeneRuler_100bp_DNALadder.pdf) 
- Restriction Enzyme [XhoI R0146Sk](https://github.com/hputnam/Putnam_Lab_Notebook/blob/master/images/XhoI%20_%20NEB_R0146S.pdf) 
- Restriction Enzyme Buffer [rCutSmart B6004S](https://github.com/hputnam/Putnam_Lab_Notebook/blob/master/images/rCutSmart_Buffer%20_%20NEB_B6004S.pdf) 
 
# 20240325 

## Pocillopora 

*Note that we ended up re doing this run on 20240326 because we added a positive control using a sample that was successfully amplified previously. The post for today details what we did, but the protocol did change slightly on subsequent days.*  

Danielle prepared samples and pulled the samples we will need to use from the freezer. 

### Sample metadata 

| Sample ID | Tube No | Extracted Date |
|-----------|---------|----------------|
| POC 200   | 273     | 20211007       |
| POC 205   | 289     | 20210924       |
| POC 207   | 309     | 20211105       |
| POC 215   | 335     | 20211115       |
| POC 217   | 81      | 20211015       |
| POC 238   | 321     | 20211018       |
| POC 239   | 317     | 20211104       |
| POC 248   | 325     | 20211007       |
| POC 254   | 311     | 20211019       |
| POC 366   | 363     | 20210916       |
| POC 386   | 465     | 20211005       |
| POC 395   | 447     | 20211122       |
| POC 41    | 407     | 20211015       |
| POC 45    | 383     | 20211001       |
| POC 55    | 391     | 20211118       |
| POR 362   | 707     | 20211020       |
| POR 83    | 469     | 20220208       |

15 POC samples  
2 POR samples  

1 negative control 

Today we are just doing the POC samples, so there are 16 total to run.  

### POC Master Mix 

We are planning the master mix for the POC samples for 16 samples ± 5-10% for error. We are planning for 18 total reactions to have enough master mix volume. 

| Component                      | Concentration | 1 rxn (uL) | 18 rxns (uL) |
|--------------------------------|---------------|------------|--------------|
| EmeraldAmp 2X GT PCR Master Mix |               | 12.55      | 225.9        |
| FatP6.1 F primer               | 10 uM         | 0.32       | 5.76         |
| RORF R primer                  | 10 uM         | 0.32       | 5.76         |
| Ultra Pure H20                 |               | 10.8       | 194.4        |
| DNA                            |               | 1          |            |
| Total volume                   |               | 25uL       |              |

Our primers are:  

FatP6.1 5′-TTTGGGSATTCGTTTAGCAG-3′  
RORF 5′-SCCAATATGTTAAACASCATGTCA-3′  

### POC PCR amplification  

Master mix was made (without template DNA) by Danielle by first adding water, then primers, then the master mix.  

Labeled PCR tubes with tube number in ascending order, with negative control as the last tube. Used PCR strip tubes. 

Added 24uL of master mix made above to each tube. Added 1uL Ultrapure H20 as the negative control. Added 1 uL template DNA to each sample tube from each DNA extraction tube.  

Samples added to the thermocycler at 11:30 for PCR, completed at 13:00.  

The PCR protocol was as follows:  

- [94°C 60 secondes] 1 cycle
- [94°C 30 sec,53°C 30 sec, 72°C 60 sec] 30 cycles
- [72°C 5 minutes] 1 cycle
- [4°C infinity]
- stored at 4°C until gel was run and we proceeded with restriction enzyme incubations 

This PCR program is already programmed in the PPP lab thermocyclers.  

### POC Gel to QC PCR product 

We next ran a gel to QC the PCR product. We prepared a 1% agarose gel in a medium box (20 wells) with 1XTAE buffer and 1uL of Gel Green added. Gel was run at 80V for 30 minutes from 14:56-15:26. 4uL of sample were added to each well with 2 ladders on either side (1KB Purple ladder 100uM).  

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/uri_2024/20240325_POC_QC_gel.jpg?raw=true)

The gel shows a consistent band length at 669bp as expected and the negative control is clear. 

### POC Restriction enzyme RFLP incubation 

We then prepared the restriction enzyme master mix. We first diluted 10X rCutSmart to 1X (1uL rCutSmart and 9 uL Ultrapure H20). We planned for 20 reactions (18 samples + 10% error).   

| Component | Concentration | 1 rxn | 20 rxns |
|-----------|---------------|-------|---------|
| rCutSmart | 1X            | 0.5uL | 10uL    |
| Xho1 RE   |               | 0.5uL | 10uL    |
| Total     |               |       | 20uL    |

We then added 1uL of master mix and 15uL of PCR product into new strip tubes for each sample.  

Tubes were loaded into a thermocycler and run for 1 h at 37°C followed by 20 min at 65°C with hold at 10°C.  

### POC Restriction enzyme RFLP gel

After the incubation, samples were removed and loaded into a 2% gel and ran at 70V for 1 h from 16:39-17:39.  

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/uri_2024/20240325_POC_RFLP_gel.jpg?raw=true)  

There were no samples that were cut in this RFLP, which could mean they were all *P. meandrina*. However, without a positive control, we cannot be sure that the restriction enzymes worked, which could generating false identification of *P. meandrina*. We decided that tomorrow we will include previously amplified and RFLP'd samples from [Hollie's previous RFLP](https://hputnam.github.io/Putnam_Lab_Notebook/PocID/).  

# 20240326  

## Pocillopora 

### POC Histone PCR 

Today we re did samples ran yesterday but with a positive control. We are using sample POC-198 that successfully cut and was therefore identified as *P. grandis* in [Hollie's previous RFLP](https://hputnam.github.io/Putnam_Lab_Notebook/PocID/).  

All other preparations an dprotocols were the same as described above, but with increased master mix volume for the one additional sample.  

We started the PCR at 10:15 and it was completed by 11:50. 

### POC Gel PCR QC 

We ran a 2% gel to QC POC Histone PCR product. We ran the gel at 80 V for 1h with 4uL of each sample and two ladders (same ladder as described above). Gel ran from 13:02-14:02.  

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/uri_2024/20240326_POC_QC_gel.jpeg?raw=true)  

We noticed that we have a band in our negative control. However, we think that this is a loading dye artifact because the same band is not seen in our other samples. We will try to fix this in our next round of PCR. The positive control sample amplified as did our other samples as expected.  

### POC Restriction enzyme RFLP incubation

We conducted the restriction enzyme incubation as done yesterday.  

| Component | Concentration | 1 rxn | 20 rxns |
|-----------|---------------|-------|---------|
| rCutSmart | 1X            | 0.5uL | 10uL    |
| Xho1 RE   |               | 0.5uL | 10uL    |
| Total     |               |       | 20uL    |

We then added 1uL of master mix and 15uL of PCR product into new strip tubes for each sample.  

Tubes were loaded into a thermocycler and run for 1 h at 37°C followed by 20 min at 65°C with hold at 10°C from 13:32-15:00.  

We also added another negative control (NEG2) that was the restriction enzyme master mix with 15 uL of H20.  

### POC RFLP gel

We ran a 2% gel at 80V for 1 h to view RFLP results. Gel ran from 13:45-14:45. We imaged the gel at 14:07. 1uL loading dye was added to NEG2 before adding to the gel. 

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/uri_2024/20240326_POC_RFLP_gel.jpeg?raw=true)  

We again saw no cutting in any of our samples and no cutting in the positive control. This suggests that that either one of two things is happening:  

1. The POC-198 colony is actually not *P. grandis* and is therefore not cutting. 
2. The restriction enzymes are bad or the protocol is not correct.  

We checked the protocol multiple times against Hollie's posts, Zoe's posts, and the Johnston et al. 2018 paper and could not find any differences in our approach and other approaches. It is more likely that there is an issue with the restriction enzymes.  

## Porites 

### POR H2 PCR

Today we also started the process to PCR and prepare POR samples for Sanger sequencing using the H2 marker (see resources above).  

We are planning the master mix for the POR samples for 2 samples, 2 positive controls, and 1 negative ± 5-10% for error. We are planning for 6 total reactions to have enough master mix volume. 

| Component                      | Concentration | 1 rxn (uL) | 6 rxns (uL) |
|--------------------------------|---------------|------------|--------------|
| EmeraldAmp 2X GT PCR Master Mix |               | 12.55      | 75.3        |
| zH2AH4f             				 | 10 uM         | 0.32       | 1.92         |
| zH4Fr R primer                  | 10 uM         | 0.32       | 1.92         |
| Ultra Pure H20                 |               | 10.8       | 64.8        |
| DNA                            |               | 1          |            |
| Total volume                   |               | 25uL       |              |

Our primers are:  

zH2AH4f 5′-GTGTACTTGGCTGCYGTRCT-3′
zH4Fr 5′-GACAACCGAGAATGTCCGGT-3′

Our positive controls are: 
 
- POR 209: *P. lobata/lutea* from the E5 time series previously sequenced (POS1)
- POR 72: *P. evermanni* from the E5 time series previously sequenced (POS2)

Our samples are:  

- POR 362 (tube 707)
- POR 83 (tube 469)  

We made the master mix and added 1uL template DNA to each tube, with 1uL H20 added to the negative control.  

Samples were added to the thermocycler at 10:54-12:41. PCR was run using the same protocol as described above.  

The PCR protocol was as follows:  

- [94°C 60 secondes] 1 cycle
- [94°C 30 sec,53°C 30 sec, 72°C 60 sec] 30 cycles
- [72°C 5 minutes] 1 cycle
- [4°C infinity]
- stored at 4°C until prepared for sequencing  

### POR H2 QC Gel 

We ran a 2% gel as described above for the POR H2 PCR product.  

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/uri_2024/20240326_POR_QC_gel.jpeg?raw=true)  

- We found that the POR 707 sample did not amplify (POR-362). This sample also did not amplify on previous attempts to sequence this sample. 
- Both positive controls amplified and negative control was clear. 
- 469 amplified from POR-83. In previous attempts for E5 processing this sample did amplify but failed sequencing. We will have to see what the sequencing result is to know if it was successful.  

POR-362 was taken from tube 707 from time point 3. Notes in the E5 notebook say that this sample was taken from a physiology fragment, rather than a molecular tube. We will do another round of PCR where we take a different tube from this same colony from a different time point.  

It could also be that DNA concentrations were too high or there are inhibitors present. To investigate this, we will also run diluted samples in our next round of PCR. 

### POR H2 PCR Attempt #2 today 

We ran other PCR today and made two changes: 

1. We diluted sample 707 by 1/10 and 1/100
2. We selected a new sample for colony POR-362 and added a sample of the full concentration, 1/10, and 1/100 dilutions. 

We chose tube 837 as our new sample for POR-362. This sample is from time point 4 rather than time point 3.  

We have 8 samples for this PCR:  

- 707 diluted 1/10
- 707 diluted 1/100
- 837 full concentration
- 837 diluted 1/10
- 837 diluted 1/100
- POS1 (*P. lobata/lutea* sample 283)
- POS2 (*P. evermanni* sample 493)
- Negative control  

We made a master mix, loaded samples into strip tubes for PCR as described above, and ran protocol from 16:46-18:31. PCR product was stored at 4°C overnight until, checking on a gel tomorrow morning.  

# 20240327 

## Pocillopora 

### POC PCR 

We ran a new round of PCR and RFLP today using three new positive control samples that were identified as *P. grandis* from POC spawning work in Moorea.  

- POC-150 (POS1)
- POC-156 (POS2)
- POC-159 (POS3)
- POC-198 (POS4)

We chose to include POC-198 because it *should* be *P. grandis*, and we want to resolve this discrepancy.  

We prepared the PCR as described above for 22 reactions and ran the PCR from 10:15-11:48.    

### POC PCR QC Gel 

We ran another 2% gel as described above.  

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/uri_2024/20240327_POC_QC_gel.jpeg?raw=true)

We had amplification as expected from our samples and the POC-198 positive control (POS4). Our negatives are clear. Unexpectedly, there was not much amplification of the POS1-3 controls. The band also appears to be shorter. It is possible that this specific region is a shorter sequence in *P. grandis* and that is why the bands are shorter. But if POC-198 is also *P. grandis*, it should have the same banding pattern and it does not - it is the same as our other samples. 

### POC RFLP incubation 

We then took the PCR product, added the restriction enzyme master mix as described above, and incubated with restriction enzymes. We have 6 uL of PCR product left over after using 4uL for QC and 15uL for incubations. Left over PCR product was kept at 4°C.  

Incubations took place from 12:45-14:05.  

After this incubation, we ran another gel. 

### POC RFLP Gel QC 

We ran a 2% gel from 14:44-15:44 for 1h at 80V as described above.  

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/uri_2024/20240327_POC_RFLP_gel.jpeg?raw=true)

There was no cutting in any sample, even the positive controls. We concluded that the restriction enzyme is the problem. Either the reagents are degraded (which is unlikely with an expiration of 2025), or something is wrong with the master mix.  

We looked back at Hollie's notebooks and in [her post](https://hputnam.github.io/Putnam_Lab_Notebook/PocID/) she uses the same master mix and protocol that we did. But after talking with Hollie, she believes that she likely did not perform the dilution from 10X to 1X of the CutSmart buffer. In the case that a higher concentration is required of this buffer, we decided to do another RFLP incubation with 10X CutSmart rather than 1X to test this.  

### POC RFLP incubation with remaining PCR product using higher CutSmart concentration  

We then prepared the following restriction enzyme master mix:  

- 12 uL of Xhol
- 12 uL of 10X rCutSmart buffer 

1uL of this mix was added to each remaining PCR product (6uL) along with 9uL ultrapure water (to bring sample volume to 15 uL). This resulted in a diluted PCR product, but because our bands are so bright, dilution should not reduce the concentration too low that we will not be able to see the banding pattern.  

These were added to the thermocycler for incubation protocol at 17:20-18:40.  

We then prepared a 2% gel and ran the gel at 80V for 50 min.  

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/uri_2024/20240327_POC_RFLP_gel_2.jpeg?raw=true)

We did not see cutting in this gel. This means that the restriction enzyme is the problem.  

We will run a new PCR tomorrow at a larger volume and with new enzymes that are scheduled to arrive tomorrow. If the enzymes don't work again, we will prepare PCR product for sequencing as the alternative option.  

## Porites 

### POR PCR gel QC  

We ran a 2% gel on yesterday's PCR product at 70V for 1h at 10:20-11:20 as described above.  

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/uri_2024/20240327_POR_QC_gel.jpeg?raw=true)

There was no amplification of sample 707 again, even with dilution. This sample may just be degraded. There was, however, amplification of sample 837 at all dilutions!   

We will move forward with preparing the full concentration of 837 and 469 for Sanger sequencing.    

### PCR product clean up 

xxx

### Preparation for Sanger sequencing 

Danielle prepared the samples for Sanger sequencing for Janet at RI-INBRE.  

First, dilute the two primer stock for sequencing (both F and R).  

10 uM primer stock diluted to 3.2 uM:  

1.6 uL of 10 uM primer + 3.4 uL H2) = 5 uL 

Samples were diluted 1:1 with ultra pure H20 with 2uL H20 and 2uL PCR product 

| Sample ID | Well (GSC use only) | Template type | A. Template size (bases) | B. Template stock conc. (ng/ul) | C. PCR template: ng needed = ((A/100)*1.25)*2 | D. PCR template: Volume=(C/B)ul | E. Plasmid template: volume= (200ng/B)*2 | F. Volume PCR H20 needed (10-D or E) | G. Volume primer needed 1ul per rxn |                                            |             |   |
|:---------:|:-------------------:|:-------------:|:------------------------:|:-------------------------------:|:---------------------------------------------:|:-------------------------------:|:----------------------------------------:|:------------------------------------:|:-----------------------------------:|--------------------------------------------|-------------|---|
|           |                     |               |                          |                                 |                                               |                                 |                                          |                                      |                                     | Actual   PCR product nanodop concentration |             |   |
|           |                     |               |                          |                                 |                                               |                                 |                                          |                                      |                                     |                                            |             |   |
|           |                     |               |                          |                                 |                                               |                                 |                                          |                                      |                                     |                                            | Colony name |   |
|     1     |                     |      PCR      |           1500           | 10.85                           |                      37.5                     |               3.5               |                                          |                  6.5                 | 2                                   | 21.7                                       | POR-469     | F |
|     2     |                     |      PCR      |           1500           | 10.85                           |                      37.5                     |               3.5               |                                          |                  6.5                 | 2                                   | 21.7                                       | POR-469     | R |
|     3     |                     |      PCR      |           1500           | 19.9                            |                      37.5                     |               1.9               |                                          |                  8.1                 | 2                                   | 39.8                                       | POR-837     | F |
|     4     |                     |      PCR      |           1500           | 19.9                            |                      37.5                     |               1.9               |                                          |                  8.1                 | 2                                   | 39.8                                       | POR-837     | R |

Danielle labeled the tubes and we will take them to Janet for sequencing tomorrow.  

# 20240328  

## Pocillopora 

Today we did a new round of PCR and RFLP with new restriction enzymes.  

Danielle ran this protocol today and Ariana ran the final gel.   

### POC PCR 

xxx

### POC PCR QC gel

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/uri_2024/20240328_POC_QC_gel.png?raw=true)

### POC RFLP incubation 

xxx

### POC RFLP gel 

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/uri_2024/20240328_POC_RFLP_gel.png?raw=true)
 
xxx

## Porites

Danielle and Jill delivered the Porites prepared product for sequencing at RI-INBRE!  

# 20240329

## Pocillopora

### POC RFLP incubation without diluted rCutSmart

### POC RFLP gel 
