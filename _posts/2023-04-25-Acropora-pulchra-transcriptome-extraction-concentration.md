---
layout: post
title: Concentrating Gametogenesis Acropora pulchra RNA for De Novo Transcriptome Assembly
date: '2023-04-25'
categories: Protocols
tags: RNA, DNA, Protocol, Gametogenesis, Transcriptome
projects: Putnam Lab
---

#### Goal:
Use *Acropora pulchra* samples from Mo'orea, French Polynesia, sampled on January 15th 2022 from the north shore backreef site Mahana (17°29'13.9"S 149°53'14.7"W) to create a de novo transcriptome for *A. pulchra*. Extracted samples (August 23rd 2022) where concentrated using the Zymo [RNA Clean and Concentrate Kit](https://www.zymoresearch.com/products/rna-clean-concentrator-5) on April 25th 2023 at URI.

### **Sample Collection**

Samples were collected and preserved in DNA-RNA shield in Mo'orea, French Polynesia following the [Sample Collection Protocol](https://github.com/daniellembecker/Gametogenesis/blob/main/protocols/2021-12-26-Sample_Same_Day_Processing_Protocol.md). Initial samples had 1-2 mL of DNA RNA shield added.

### **Sample Information**

January 2022 samples were brought back with CITES permit 20220819 and stored at URI. Samples DNA/RNA was extracted on August 23rd 2022 following this [exact protocol](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2022-08-23-Gametogenesis-January%202022-DNA-RNA-Extractions.md).

# **Sample IDs**

## **11 January 2022 samples used for transcriptome**

| timepoint    | tube.number | collection.date | colony.ID | uL taken |
|--------------|-------------|-----------------|-----------|----------|
| January 2022 |     418     | 20220115        | ACR-418   | 300      |
| January 2022 |     422     | 20220115        | ACR-422   | 300      |
| January 2022 |     428     | 20220115        | ACR-428   | 300      |
| January 2022 |     432     | 20220115        | ACR-432   | 300      |
| January 2022 |     438     | 20220115        | ACR-438   | 300      |
| January 2022 |     457     | 20220115        | ACR-457   | 300      |
| January 2022 |     458     | 20220115        | ACR-458   | 300      |
| January 2022 |     459     | 20220115        | ACR-459   | 300      |
| January 2022 |     460     | 20220115        | ACR-460   | 300      |
| January 2022 |     464     | 20220115        | ACR-464   | 300      |
| January 2022 |     465     | 20220115        | ACR-465   | 300      |

![image](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/images/20220823_extractions.jpeg)


### DNA/RNA Quantity - 11 January 2022 Samples

| tube.number    | colony.ID | qubit.reading.1 | qubit.reading.2 | quibit.average | sample.type | timepoint |
|----------------|-----------|-----------------|-----------------|----------------|-------------|-----------|
| DNA_standard_1 | NA        |          200.58 | NA              |         200.58 | DNA.STD     | NA        |
| DNA_standard_2 | NA        |           24151 | NA              |          24151 | DNA.STD     | NA        |
|              1 | ACR-418   |            25.6 |            25.4 |           25.5 | DNA         | JANUARY   |
|              2 | ACR-422   |              24 |            23.4 |           23.7 | DNA         |  JANUARY   |
|              3 | ACR-428   |            19.8 |            19.3 |          19.55 | DNA         | JANUARY   |
|              5 | ACR-432   |            15.1 |            14.8 |          14.95 | DNA         | JANUARY   |
|              6 | ACR-438   |              25 |            24.4 |           24.7 | DNA         | JANUARY   |
|              7 | ACR-457   |            27.6 |              27 |           27.3 | DNA         | JANUARY   |
|              8 | ACR-458   |            29.8 |              29 |           29.4 | DNA         | JANUARY   |
|              9 | ACR-459   |            13.1 |              13 |          13.05 | DNA         | JANUARY   |
|             10 | ACR-460   |            15.6 |            15.1 |          15.35 | DNA         | JANUARY   |
|             11 | ACR-464   |            19.2 |            18.8 |             19 | DNA         | JANUARY   |
|             12 | ACR-465   |              14 |            13.5 |          13.75 | DNA         | JANUARY   |
| RNA_standard_1 | NA        |          422.43 | NA              |         422.43 | RNA.STD     | NA        |
| RNA_standard_2 | NA        |         9892.36 | NA              |        9892.36 | RNA.STD     | NA        |
|              1 | ACR-418   |            31.6 |            31.6 |           31.6 | RNA         | JANUARY   |
|              2 | ACR-422   |            21.6 |              22 |           21.8 | RNA         | JANUARY   |
|              3 | ACR-428   |            27.6 |            27.6 |           27.6 | RNA         |  JANUARY   |
|              5 | ACR-432   |            16.8 |            16.4 |           16.6 | RNA         | JANUARY   |
|              6 | ACR-438   |            15.2 |            14.4 |           14.8 | RNA         | JANUARY   |
|              7 | ACR-457   |            18.4 |            18.2 |           18.3 | RNA         | JANUARY   |
|              8 | ACR-458   |            15.2 |            15.2 |           15.2 | RNA         | JANUARY   |
|              9 | ACR-459   |            17.4 |            16.6 |             17 | RNA         | JANUARY   |
|             10 | ACR-460   |            20.4 |            20.6 |           20.5 | RNA         | JANUARY   |
|             11 | ACR-464   |            13.8 |            13.8 |           13.8 | RNA         | JANUARY   |
|             12 | ACR-465   |            13.2 |              13 |           13.1 | RNA         | JANUARY   |

### **Clean and Concentrate Protocol**

For our desired sequence volume to submit to Genewiz for RNAseq for transcriptome analyses, we needed one sample with >50ng/uL. Using the 11 samples above to get a diverse range of genotypes and RNA form multiple individuals, we followed the [Zymo RNA Clean Concentrate Protocol](https://github.com/zdellaert/ZD_Putnam_Lab_Notebook/blob/master/protocols/Zymo_RNA_Clean_Concentrate.pdf) **based on the volume for each tube** for the Zymo Clean and Concentrate kit, eluting in 25 µL of DNAse/RNAse free water.

10 uL of each sample (10 uL * 11 samples = 110 uL) was added to a 1.5 mL microcentrifuge tube 

Base volume: 110 uL

1. Since each sample is **110** uL, add **220** uL (2 x **110**) of RNA binding buffer.
2. Add an equal volume (**330** uL (3 x **110**)) of 100% ethanol and mix.
3. Transfer the sample to the Zymo-Spin™ IC Column in a Collection Tube and centrifuge. Discard the flow-through.
4. Add 400 µl RNA Prep Buffer to the column and centrifuge. Discard the flow-through.
5. Add 700 µl RNA Wash Buffer to the column and centrifuge. Discard the flow-through.
6. Add 400 µl RNA Wash Buffer to the column and centrifuge for 1 minute to ensure complete removal of the wash buffer. Carefully, transfer the column into a RNase-free tube (not provided).
7. Add 25 µl DNase/RNase-Free Water directly to the column matrix and centrifuge.

## Qubit Results

- Used High Sensitivity range dsDNA and RNA Qubit [Protocol](https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/Qubit-Protocol/)
- The concentrated sample was read twice, the standards were only read once

 RNA Standards: 418.23 (S1) & 8244.21 (S2)

### Concentrating all sample

| colony_id | RNA_QBIT_AVG | Volume in Tube currently | Amount RNA in whole tube |
|-----------|--------------|--------------------------|--------------------------|
| ACRP   | 43.2        | 24                       | 1036.8                  |
