---
layout: post
title: Workflow for Acropora pulchra de novo transcriptome
category: [ de novo transcriptome , DNA]
tag: [ Acropora pulchra, de novo transcriptome ]
---
## Designing a workflow to create a de novo transcriptome for Acropora pulchra

#### Goal:
Use one *Acropora pulchra* [concentrated sequence sample]() from 11 colonies collected in Mo'orea, French Polynesia on January 15th 2022 from the north shore backreef site Mahana (17°29'13.9"S 149°53'14.7"W) part of a 12-month [Gametogenesis timeseries project]() and five sequence samples also collected from Mo'orea, French Polynesia part of the [E5 Rules of Life project]() to create a de novo transcriptome for *A. pulchra*. Literature review of current *Acropora* de novo transcriptomes and genomes completed already.

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
  - Sperm DNA obtained from a single colony of the coral Acropora digitifera was used for genome sequencing by Roche 454 GS-FLX10 and Illumina Genome Analyser IIx (GAIIx)11.
  - The 454 shotgun and paired-end reads were assembled de novo by GS De novo Assembler version 2.3 (Newbler, Roche)10, and subsequent scaffolding was performed by SOPRA27 and SSPACE28 using the Illumina mate-pair information.
  - A set of gene model predictions (the A. digitifera Gene Model v. 1) was generated mainly by AUGUSTUS29, and a genome browser has been established using the Generic Genome Browser (GBrowser) 2.17.

- [Shinzato et al. 2020](https://academic.oup.com/mbe/article/38/1/16/5900672)
  - Sequenced genomes of 15 *Acropora* species (*A. acuminata*, *A. awi*, *A. cytherea*, *A. digitifera*, *A. echinata*, *A. florida*, *A. gemmifera*, *A. hyacinthus*, *A. intermedia*, *A. microphthalma*, *A. muricata*, *A. nasta*, *A. selago*, *A. tenuis*, and *A. yongei*)
  -  We assessed completeness of genome assembly with BUSCO ver. 3.0.2 (Simao et al. 2015; Waterhouse et al. 2018) and the Metazoan set (978 genes).
  - Repetitive elements in the draft genomes of Acroporidae (Acropora, Montipora, and Astreopora) were identified de novo with RepeatScout v.1.0.5 (Price et al. 2005) and annotated with BlastN and BlastX searches against RepeatMasker.lib and RepeatPeps.lib bundled with RepeatMasker v.4.0.6 (Smit et al. 1996–2010), as reported in Luo et al. (2015, 2018).
  - For nonannotated or putative novel repeats, one additional class “Novel” was introduced.
  - The expansion history of repetitive elements were calculated and visualized using perl scripts from RepeatMasker package (calcDivergenceFromAlign.pl and createRepeatLandscape.pl) as reported in Khalturin et al. (2019).
  - gene prediction of Acropora genome assemblies, Augustus version 3.2.3 (Stanke et al. 2006) was first trained using 2,000 high-quality A. digitifera-assembled transcriptome sequences (Shinzato et al. 2011) selected by PASA (Haas et al. 2003).
  - Then, the trained Augustus was used for gene prediction from repeat-masked genome assemblies produced by RepeatMasker (Smit et al. 1996–2010) together with the A. tenuis and A. digitifera RNA-Seq data as gene structure hints.
  - We also assessed completeness of repertoires of predicted genes (mRNA) using BUSCO with the “transcriptome” setting.
  - All proteomes were BLASTed against the Uniprot/Swissprot (UniProt Consortium 2018) database and were analyzed with InterProScan 5 (Jones et al. 2014).
  - Genome browsers for the 18 acroporid genomes are available from the Marine Genomics Unit web site (https://marinegenomics.oist.jp/gallery).

#### *Acropora* de novo assemblies and notes

- [Oldach and Vize 2018](https://www.sciencedirect.com/science/article/pii/S1874778717303422?via%3Dihub)
  - De novo assembly and annotation of the *Acropora gemmifera* transcriptome
  - Used Trintiy to assemble 31.6 million combined raw reads and built into 104,000 contigs
  - Functional gene annotation was performed using dammit, Gene Ontology (GO), KOG (WebMGA) and KEGG pathway analyses (Kaas)
  - Raw reads were deposited at the NCBI Sequence Read Archive (SRA Accession SRX1629970 and SRX1629969, BioProject PRJNA293100). Unannotated transcript contigs were deposited to the Transcriptome Shotgun Assembly (TSA) database in GenBank (Accession no. GFQE00000000). The annotated transcriptome (FASTA and GFF3) is available on Mendelay.

  - [Kitchen et al. 2015](https://academic.oup.com/g3journal/article/5/11/2441/6025398)
    - De Novo Assembly and Characterization of Four Anthozoan (Phylum Cnidaria) Transcriptomes
    - Useed putative gene names and functional categories [Gene Ontology (GO) and Kyoto Encyclopedia of Genes and Genomes (KEGG)] to assembled transcripts based on sequence comparisons with online databases.
    - All sequence comparisons were conducted using BLAST+ from National Center for Biotechnology Information (NCBI) (Package version 2.2.29) (Altschul et al. 1990).
    - Gene names were assigned by comparing transcript sequences against UniProt protein sequence databases (SwissProt and TREMBL) using BLASTx with an expect value (E value) cutoff of 10−4. Each transcript was assigned a gene name based on its best match, excluding matches with uninformative names (e.g., uncharacterized, unknown, or hypothetical).

#### Other de novo transcriptome resources

- [Wong and Putnam *Porites astreoides* genome](https://gigabytejournal.com/articles/65)
  - [GitHub](https://github.com/hputnam/Past_Genome)
  - Structural annotation of the P. astreoides genome was completed on the University of Rhode Island High Performance Computer ‘Andromeda’. As input for MAKER v3.01.03 (RRID:SCR_005309) [64] we used an existing P. astreoides transcriptome from samples collected in the Florida Keys, USA [32] and existing congener P. lutea peptide sequences from a sample collected in Australia [57﻿].
  - These files were used as input for an initial round of MAKER to predict gene models directly from this transcriptomic and protein data, respectively. The first round of MAKER included RepeatMasker (RRID:SCR_012954) [65].
  - The output from the initial MAKER round was used to train ab initio gene predictors SNAP (RRID:SCR_002127) [66] and AUGUSTUS (RRID:SCR_008417) [67] using BUSCO v5.2.2(RRID:SCR_015008) [63].
  - A second round of MAKER was performed using the general feature format (GFF) output of the first round containing the information on the locations of repetitive elements for masking, as well as the locations of expressed sequence tags (ESTs) and proteins and ab initio gene prediction from the SNAP and AUGUSTUS outputs.
  - Another round of ab initio gene prediction was performed on the output from the second round of MAKER.
  - A third and final round of MAKER was conducted, including training information from the GFF files generated from the second round of MAKER and ab initio gene prediction as input. Genome structural annotations were compared against all four Porites congeners using AGAT v0.8.1 (Table 2).

- [Chui et al. 2020](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-07113-9#Sec16)
  - "De novo transcriptome assembly from the gonads of a scleractinian coral, *Euphyllia ancora*: molecular mechanisms underlying scleractinian gametogenesis"
  - Clean reads of 24 individual samples were assembled de novo using Trinity v2.0.6 software [137] (parameter settings: –min_contig_length 150 –CPU 8 –min_kmer_cov 3 –min_glue 3 –bfly_opts ‘-V 5 –edge-thr = 0.1) and assembled sequences were clustered using Tgicl v2.0.6 software [138] (parameter settings: -l 40 -c 10 -v 25 -O ‘-repeat_stringency 0.95 -minmatch 35 -minscore 35′)
  - All assembled sequences were aligned to available genomic databases of 4 scleractinian corals and 6 Symbiodiniaceae transcriptomic databases using BLASTN (−evalue 1e-3). These databases included A. digitifera [31, 139], P. damicornis [35, 140], S. pistillata [34], and O. faveolata [141], Symbiodinium sp. A1 [142], Symbiodinium sp. A2 [143], Breviolum sp. B2 [143], Breviolum muscatinei [144], Uncultured Cladocopium sp. [145] and uncultured Durusdinium sp. [145] (For more detailed information on the databases, see Additional file 6).
  - Contigs aligned exclusively to the coral genome database were annotated as “E. ancora contigs”, while those that aligned only to Symbiodiniaceae transcriptome databases were annotated as “Symbiodiniaceae contigs”.
  - To separate contigs aligned to both the coral genome and Symbiodiniaceae transcriptomic databases, contigs were re-aligned BLASTN (−evalue 1e-3) using a combined database of coral genomes and Symbiodiniaceae transcriptomes.
  - Based on the top hit results of BLASTN (corals or Symbiodiniaceae), contigs were annotated as “E. ancora contigs” or “Symbiodiniaceae contigs” (Fig. 2)
  - Nucleotide sequences were again clustered using CD-HIT [146] with 97% identity for removing sequences possibly originating from different individuals or haplotypes in a single individual.
  - Finally, contigs were translated into amino acid sequences using the longorf script [147] and clustered using CD-HIT with 95% identity. Completeness using the assembled sequences was assessed using BUSCO (bench-marking universal single-copy orthologs) version 3 [148, 149] in transcriptome mode. Reference E. ancora gonadal transcriptome contigs were annotated as follows: 1) BLAST searches against public protein databases: SWISS-PROT database (−evalue 1e-5) (Consortium 2011) (3/18/2019), 2) Identification of conserved protein domains with the Pfam database (−evalue 1e-5) [150, 151].


#### *Acropora pulchra* Gametogenesis transcriptome data files on URI andromeda:

'/data/putnamlab/KITT/hputnam/20230825_Bermuda_Reference_Transcriptomes/'

'ACRP_R1_001.fastq.gz
ACRP_R2_001.fastq.gz
ACRP_R2_001.fastq.gz.md5'

#### *Acropora pulchra* E5 Rules of Life project transcriptome data files on
