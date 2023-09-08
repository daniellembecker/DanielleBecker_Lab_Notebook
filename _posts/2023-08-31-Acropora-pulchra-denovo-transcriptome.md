---
layout: post
title: Workflow for Acropora pulchra de novo transcriptome
category: [ de novo transcriptome , DNA]
tag: [ Acropora pulchra, de novo transcriptome ]
---
## Designing a workflow to create a de novo transcriptome for *Acropora pulchra*

#### Goal:
Use one *Acropora pulchra* [concentrated sequence sample](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2023-04-25-Acropora-pulchra-transcriptome-extraction-concentration.md) from 11 colonies collected in Mo'orea, French Polynesia on January 15th 2022 from the north shore backreef site Mahana (17°29'13.9"S 149°53'14.7"W) part of a 12-month [Gametogenesis timeseries project](https://github.com/daniellembecker/Gametogenesis) and five sequence samples also collected from Mo'orea, French Polynesia part of the [E5 Rules of Life project](https://github.com/urol-e5) to create a de novo transcriptome for *A. pulchra*. Literature review of current *Acropora* de novo transcriptomes and genomes completed already.

**Important notes about de novo transcriptomes**

- [Raghavan et al. 2022](https://academic.oup.com/bib/article/23/2/bbab563/6514404#)
  - "A simple guide to de novo transcriptome assembly and annotation"
  - De novo transcriptome assembly, in contrast, is ‘reference-free’. The process is de novo (Latin for ‘from the beginning’) as there is no external information available to guide the reconstruction process. It must be accomplished using the information contained in the reads alone.
  - This approach is useful when a genome is unavailable, or when a reference-guided assembly is undesirable.
  - For instance, in opposition to a de novo assembler successfully producing a transcript, a reference-guided approach might not be able to reconstruct it correctly if it were to correspond to a region on the reference containing sequencing or assembly gaps [15, 16].

![figure1](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/images/transcriptome.png)

#### *Acropora* genomes

- [Shinzato et al. 2011](https://www.nature.com/articles/nature10249)
  - First *Acropora* genome sequenced was *Acropora digitifera* in "Using the Acropora digitifera genome to understand coral responses to environmental change"

- [Shinzato et al. 2020](https://academic.oup.com/mbe/article/38/1/16/5900672)
  - Sequenced genomes of 15 *Acropora* species (*A. acuminata*, *A. awi*, *A. cytherea*, *A. digitifera*, *A. echinata*, *A. florida*, *A. gemmifera*, *A. hyacinthus*, *A. intermedia*, *A. microphthalma*, *A. muricata*, *A. nasta*, *A. selago*, *A. tenuis*, and *A. yongei*)

#### *Acropora* de novo assemblies and notes

- [Oldach and Vize 2018](https://www.sciencedirect.com/science/article/pii/S1874778717303422?via%3Dihub)
  - De novo assembly and annotation of the *Acropora gemmifera* transcriptome
  - Used Trintiy to assemble 31.6 million combined raw reads and built into 104,000 contigs

  - [Kitchen et al. 2015](https://academic.oup.com/g3journal/article/5/11/2441/6025398)
    - De Novo Assembly and Characterization of Four Anthozoan (Phylum Cnidaria) Transcriptomes

#### Other de novo transcriptome resources

- [Wong and Putnam *Porites astreoides* genome](https://gigabytejournal.com/articles/65)
  - [GitHub](https://github.com/hputnam/Past_Genome)
  - Structural annotation of the P. astreoides genome was completed on the University of Rhode Island High Performance Computer ‘Andromeda’. As input for MAKER v3.01.03 (RRID:SCR_005309) [64] we used an existing P. astreoides transcriptome from samples collected in the Florida Keys, USA [32] and existing congener P. lutea peptide sequences from a sample collected in Australia [57﻿].

- [Chui et al. 2020](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-07113-9#Sec16)
  - "De novo transcriptome assembly from the gonads of a scleractinian coral, *Euphyllia ancora*: molecular mechanisms underlying scleractinian gametogenesis"

#### Trinity Resources and Vignette

[Full-length transcriptome assembly from RNA-seq data without a reference genome. Grabherr et al. 2011](https://www.nature.com/articles/nbt.1883)

[De novo transcript sequence reconstruction from RNA-seq using the Trinity platform for reference generation and analysis. Haas et al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3875132/)
  - [Software, documentation, and demonstrations](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3875132/)

#### *Acropora pulchra* Gametogenesis transcriptome data files on URI andromeda:

  Location on Andromeda, the HPC server for URI:

  ```
  cd /data/putnamlab/KITT/hputnam/20230825_Bermuda_Reference_Transcriptomes/

  ACRP_R1_001.fastq.gz
  ACRP_R2_001.fastq.gz
  ACRP_R2_001.fastq.gz.md5
  ```

#### *Acropora pulchra* E5 Rules of Life project transcriptome data files on University of Washington OWL data storage HPC database

  Location on OWL, the HPC server for UW:

  ```
  https://owl.fish.washington.edu/nightingales/A_pulchra/30-789513166/

  ACR-140_TP2
  ACR-145_TP2
  ACR-150_TP2
  ACR-173_TP2
  ACR-178_TP2
  ```

#### Samples QC and Qubit Results

  | Tube Label  | RNA_ng_µl | RNA_µl | RNA µg | Link to notebook post1                                                                                                                                  |
  |-------------|-----------|--------|--------|---------------------------------------------------------------------------------------------------------------------------------------------------------|
  | ACR-140_TP2 | 12        | 90     | 1.08   | https://kterpis.github.io/Putnam_Lab_Notebook/20211012-RNA-DNA-extractions-from-E5-project/                                                             |
  | ACR-145_TP2 | 20.8      | 87     | 1.8096 | https://kterpis.github.io/Putnam_Lab_Notebook/20211012-RNA-DNA-extractions-from-E5-project/                                                             |
  | ACR-150_TP2 | 13        | 87     | 1.131  | https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-09-03-20210903-RNA-DNA-extractions-from-E5-project.md                            |
  | ACR-173_TP2 | 11.4      | 87     | 0.9918 | https://kterpis.github.io/Putnam_Lab_Notebook/20211102-RNA-DNA-extractions-from-E5-project/                                                             |
  | ACR-178_TP2 | 12.2      | 87     | 1.0614 | https://github.com/Kterpis/Putnam_Lab_Notebook/blob/master/_posts/2021-09-02-20210902-RNA-DNA-extractions-from-E5-project.md                            |
  | ACRP-CON    | 43.2      | 24     | 1.0368 | https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2023-04-25-Acropora-pulchra-transcriptome-extraction-concentration.md |

#### Workflow Steps
