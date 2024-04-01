---
layout: post
title: ITS2 amplicon PCR and preparation for sequencing 20240326
date: '2024-03-30'
categories: Larval_Symbiont_TPC_2023
tags: PCR ITS2 Molecular Mcapitata
---


**REPLACE XXX with links throughout** 




The post details ITS2 amplicon PCR amplificaiton and preparation for sequencing.  

# Overview 

We are conducting ITS2 amplicon sequencing to characterize the symbiont community in *M. capitata* larvae from the Hawaii 2023 project. See posts about this project and samples [in my previous posts here](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/categoryview/#larval-symbiont-tpc-2023).  

# Protocol 

I am following my revised and updated [ITS2 amplification protocol - IN PROGRESS](XXX) that is based on Emma's [ITS2 sequencing protocol](https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2020-01-31-ITS2-Sequencing-Protocol.md) with modifications described in this post. 

Other protocols used in this post:  

- [Gel electrophoresis](https://github.com/Putnam-Lab/Lab_Management/blob/master/Lab_Resources/DNA_RNA-protocols/Agarose-Gel-Protocol.md) 

# Sample status 

We are using samples that have already gone through DNA and RNA extraction from Jill's previous work [detailed in her notebook](https://jillashey.github.io/JillAshey_Putnam_Lab_Notebook/). gDNA has been previously quantified with data [available here](https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/data/rna_seq/rna-extraction-workbook.xlsx).   

# Equipment and supplies 

- [2X Phusion Mastermix](https://www.neb.com/en-us/products/m0531-phusion-high-fidelity-pcr-master-mix-with-hf-buffer)
- Forward primer (*ITSintfor2*) at 100uM: TCG TCG GCA GCG TCA GAT GTG TAT AAG AGA CAG GAA TTG CAG AAC TCC GTG
- Reverse primer (*ITS2_Reverse*) at 100uM: GTC TCG TGG GCT CGG AGA TGT GTA TAA GAG ACA GGG GAT CCA TAT GCT TAA GTT CAG CGG GT
- Loading Dye [NEB 6X Purple Loading Dye NEB Cat # B7024S](https://www.neb.com/en-us/products/b7024-gel-loading-dye-purple-6x)        
- Gel Stain [Biotium GelGreen Nucleic Acid Gel Stain, 10,000X in Water Fisher Cat NC9728313](https://www.fishersci.com/shop/products/gel-green-stain-5ml/NC9728313#?keyword=NC9728313)
- DNA Ladder [1kb Gel ladder](https://github.com/hputnam/Putnam_Lab_Notebook/blob/master/images/NEB_1kb_Ladder_N3232S.png?raw=true) and [100bp gel ladder](https://www.neb.com/en-us/products/n3231-100-bp-dna-ladder)
- Agarose and 1XTAE buffer for gel electrophoresis 
- Ultrapure water 
- PCR strip tubes and 1.5 mL tubes 

# 20240325

## Sample dilution 

All gDNA has previously been extracted and quantified by Jill detailed in these [posts stored on GitHub](https://github.com/AHuffmyer/larval_symbiont_TPC/tree/main/notebooks).

Today we diluted gDNA samples in preparation for ITS2 amplification and preparation for sequencing. Dilution is performed to standardize the amount of input in each reaction. Samples were diluted to a concentration of 4ng/uL with the exception of two samples, which had concentrations lower than 4ng/uL in the original gDNA. For these samples, we diluted to 3ng/uL and we will increase the input volume to have a total of 4ng gDNA input for every reaction.  

We diluted to a total volume of 10uL. For several samples, we were unable to obtain a Qubit concentration reading (see previous posts linked above), so for these we based dilutions off of the Nanodrop concentration.  

Dilutions were generated in individually labeled 0.6mL tubes and not in 96-well plates as described in Emma's base protocol. 

We will use Ultrapure H20 as our negative control. Our positive control is an Mcap2020 project sample D18, which was [successfully amplified and sequenced previously](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/ITS2-and-16s-Extractions/).  

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

## Master mix calculations  

I generated master mixes for samples that have sufficient concentration (4ng/uL) and those that had lower concentration (3ng/uL). These calculations can be found below.  

We labeled strip PCR tubes for triplicate reactions for each sample. We have a total of 69 reactions but will prepare enough master mix for a total of 76 reactions to account for pipetting error. Negative controls included the addition of Ultrapure rather in place of gDNA.  

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


**Number of reactions with low concentration (3ng/uL)** 

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

Master mixes were prepared with diluted primers as described below. Note that Phusion master mix is brand new, received today. Primers are from 2021, so these primers are older than we would like, but we will use them for this protocol. For any future ITS2 amplifications, we need to order new primers.   

**Primer dilutions** 

- Forward primer (*ITSintfor2*) at 100uM: TCG TCG GCA GCG TCA GAT GTG TAT AAG AGA CAG GAA TTG CAG AAC TCC GTG
- Reverse primer (*ITS2_Reverse*) at 100uM: GTC TCG TGG GCT CGG AGA TGT GTA TAA GAG ACA GGG GAT CCA TAT GCT TAA GTT CAG CGG GT

| Primer Dilution | Desired Volume | Current Concentration (uM) | Desired Concentration (uM) | Stock Primer Vol | H20  | Total reactions |
|-----------------|----------------|----------------------------|----------------------------|------------------|------|-----------------|
| F primer        |           38.0 |                        100 |                         10 |             3.80 | 34.2 |              76 |
| R primer        |           38.0 |                        100 |                         10 |             3.80 | 34.2 |              76 |

After preparing the master mix and thawing samples on ice, we added master mix to each tube and took independent technical replicate samples of each sample diluted gDNA into the three replicate tubes per sample.  

Tubes were spun down briefly and kept on ice.  

## PCR amplification 

Tubes were added to a thermocycler and run on the pre-programmed "ITS2" PCR protocol:  

- 1 cycle at 95°C for 3 min 
- 35 cycles at 95°C for 30 sec, 52°C for 30 sec, and 72°C for 30 sec
- 1 cycle at 72°C for 2 min
- Hold at 4°C

PCR was started at 10:48 and samples were moved to a 4°C fridge at 12:17 after PCR finished.  

## Gel for QC 

We ran two medium 2% agarose gels at 80 V for 1 h following the lab [gel protocol](https://github.com/Putnam-Lab/Lab_Management/blob/master/Lab_Resources/DNA_RNA-protocols/Agarose-Gel-Protocol.md)). 

For each gel, we mixed 100 mL 1X TAE buffer and 2 g agarose in a flask and heated for 1 min in the microwave and swirled every 20 sec. After warming, we added 1 µL GelGreen to each gel. Gel was set in a medium box with two rows of combs (40 wells in each gel total).  

Gels were allowed to cool until set. Combs were removed and trays were rotated with the top of the gel towards the black and bottom of the gel towards the red electrodes. Boxes were filled with 1X TAE buffer until a thin layer was over the gel and wells were filled.  

Samples were loaded by mixing 3 µL of each sample with 1 µL purple loading dye on parafilm and then transferring into gel wells. Pipette tips were changed between each transfer. We used 1Kb purple plus DNA ladder in our gels. We loaded all reactions with sample triplicates next to each other.  

Gen ran from 15:04-16:04 at 80V.    

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/uri_2024/20240328_ITS2_gel_1.PNG?raw=true)

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/uri_2024/20240328_ITS2_gel_2.JPG?raw=true)

Unfortunately there are double bands in every reaction, the positives, and the negatives. The most likely source of this is primer dimers, but we cannot rule out contamination. We need to run some additional testing to determine if it is due to primer dimer or contamination. We are going to rum another small set of amplifications using a new water source to start evaluating whether we have a contamination problem. We are also going to run another gel with a 100bp ladder to see if we can determine whether the bands are primer dimers (~100bp).  

In addition, samples R58 and M85 may not have amplified well. The double band may be obscuring amplified band, so it is possible there was low amplification in these samples. We will check this again after we finish the QC steps and any potential clean up.  

## PCR amplification to QC for contamination and primer dimers

To test for contamination, we will conduct another PCR with new water, which would be the most likely source of contamination since the master mix is new. If contamination is still present, we may have contamination in primers or lab equipment.  

I selected two samples that amplified well in our first PCR (R55 and R56) and the positive control (D18). We will do 1 replicate per sample for this test. We will also have another negative control with new water added. We will use a brand new nuclease free water for this PCR.  

I will also test whether we have contamination from the original DNA dilution. I will include one sample with the original dilution ("-o") and one with a new dilution at the same concentrations in the table above for these samples ("-n").  

I also used tubes that were individually sealed for all of these preps and re did the dilutions for the primers with new water. This should eliminate any source of water contamination if it is present. 

We have 7 total reactions for this QC PCR.     

**Master Mix**   

| Component            | Per Rxn (uL) | Rxns | Total Volume (uL) |
|----------------------|--------------|------|-------------------|
| 2X Phusion Mastermix |         12.5 |   7 |             87.5 |
| F primer (10uM)      |          0.5 |   7 |              3.5 |
| R primer (10uM)      |          0.5 |   7 |              3.5 |
| Ultra Pure H20       |         10.5 |   7 |             727.7 |
| DNA (4ng)            |            1 |      |                   |
| TOTAL RXN VOLUME     |           25 |      |                   |

**Sample key**  

| Sample            | Sample | Dilution | 
|----------------------|--------------|------|
| R55-o |         R55 |   original |
| R55-n      |          R55 |   new |
| R56-o |         R56 |   original |
| R56-n      |          R56 |   new |
| POS-o |         POS |   original |
| POS-n      |          POS |   new |
| NEG |         NA |   NA |

Added 24uL and 1uL template DNA (or water in the negative sample) to strip tubes.  

PCR started at 17:40 and samples were kept in 4°C fridge overnight for running the gel tomorrow.  

# 20240329

## Gel for QC of yesterday's PCR 

I ran a 2% gel at 80V for 1.5 hours and with a 100bp ladder to check for primer dimers. Samples with "-o" indicate samples with the original DNA dilution and "-n" indicate samples with new DNA dilution.  

1Kb ladder was diluted with 4 uL ladder to 16 uL Ultrapure water. This ladder stock does not include purple dye and was at a concentration of 500ng. This dilution generated a ladder stock at 100ng. 

Gel was ran from 10:30-11:00.  

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/uri_2024/20240329_ITS2_gel_1.jpeg?raw=true)

In this gel, we saw that the second band seen in our gel yesterday is at ~100bp. This further suggests that we have primer dimers, rather than contamination. Since this PCR product was made with a new water source, we are now going to run another gel including some of the original PCR product from yesterday. If bands are ~100bp in our original product, this will indicate that our double bands from our first gel are indeed primer dimers, not contamination. We may not have been able to see this in the last gel because I did not include a smaller ladder and it needed more time to separate the small band sizes. If this is the case, we will be able to proceed using this PCR product after a clean up step to remove the primer dimers. This would also indicate that we do not have a source of contamination in the lab.  

I ran another 2% gel at 80V for 1.5 hours with a 100 bp ladder. I included the same samples as above and samples from our original PCR. Samples with "-o" indicate samples with the original DNA dilution and "-n" indicate samples with new DNA dilution. Samples without these designations indicate PCR product from our first PCR.   

| Sample            | Sample | Dilution | PCR |
|----------------------|--------------|------|------|
| R55 |         R55 |   original |	first	|
| R55-o |         R55 |   original |	second	|
| R55-n      |          R55 |   new |	second	|
| R56 |         R56 |   original |	first	|
| R56-o |         R56 |   original |	second	|
| R56-n      |          R56 |   new | second		|
| POS |         POS |   original |	first	|
| POS-o |         POS |   original |	second	|
| POS-n      |          POS |   new | second		|
| NEG-1 |         NA |   NA | first		|
| NEG-2 |         NA |   NA | second		|

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/uri_2024/20240329_ITS2_gel_2.jpeg?raw=true)

In this gel we found that all samples from both PCRs and new and original dilutions had bands at ~100bp and are likely primer dimers. Negatives are clear from amplification of DNA. This is good news, we likely do not have contamination and can potentially remove the primer dimer with a size selection clean up.  

## Pooling PCR product from first PCR 

After determining that we have primer dimers and will proceed with a size selection clean up, we pooled the triplicate PCR product for each sample into new 1.5mL tubes labeled with "PCR product" and stored them in a box labeled "Huffmyer/Ashey ITS2 amplicon PCR product pooled 202403290" in Freezer D at -20°C on the shelf Putnam-2 on the right side.  

There is about 70 uL of product for each sample - plenty to work with!  

These samples are now ready for clean up, QC gel, and submission for sequencing.  

gDNA was put back in original storage boxes for Mcap larval gDNA boxes 1 and 2 in -20°C storage.  

Sample D18 was put back into the Mcap2020 larval extraction box at -20°C. 

Diluted gDNA was stored in a box in the -20°C freezer F on shelf 2 on the right side.   

# Next steps 

Next, Jill will conduct a clean up to remove primer dimers. We will use a KAPA bead clean up protocol and I will next write a template protocol and Jill will conduct the clean up, run a QC gel to confirm primer dimers are removed, and submit for sequencing at RI-INBRE. These steps will be detailed in a future post.  