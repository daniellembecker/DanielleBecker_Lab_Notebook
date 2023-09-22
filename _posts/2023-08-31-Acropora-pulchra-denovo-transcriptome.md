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
  ACRP_R1_001.fastq.gz.md5
  ACRP_R2_001.fastq.gz
  ACRP_R2_001.fastq.gz.md5
  ```

Copied all data files to new location on Andromeda

  ```

  cp -r /data/putnamlab/KITT/hputnam/20230825_Bermuda_Reference_Transcriptomes/ACRP* /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/raw

  ACRP_R1_001.fastq.gz
  ACRP_R1_001.fastq.gz.md5
  ACRP_R2_001.fastq.gz
  ACRP_R2_001.fastq.gz.md5
  ```

#### *Acropora pulchra* E5 Rules of Life project transcriptome data files on University of Washington OWL data storage HPC database

  Location on OWL, the HPC server for UW:

  ```
  https://owl.fish.washington.edu/nightingales/A_pulchra/30-789513166/

  ```

Downloaded all files to Andromeda URI HPC location

  ```
  cd /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/raw/

  wget -r -nd -A .fastq.gz https://owl.fish.washington.edu/nightingales/A_pulchra/30-789513166/
  wget -r -nd -A .md5 https://owl.fish.washington.edu/nightingales/A_pulchra/30-789513166/

  ACR-140_TP2_R1_001.fastq.gz
  ACR-140_TP2_R2_001.fastq.gz
  ACR-145_TP2_R1_001.fastq.gz
  ACR-145_TP2_R2_001.fastq.gz
  ACR-150_TP2_R1_001.fastq.gz
  ACR-150_TP2_R2_001.fastq.gz
  ACR-173_TP2_R1_001.fastq.gz
  ACR-173_TP2_R2_001.fastq.gz
  ACR-178_TP2_R1_001.fastq.gz
  ACR-178_TP2_R2_001.fastq.gz
  ACR-140_TP2_R1_001.fastq.gz.md5
  ACR-140_TP2_R2_001.fastq.gz.md5
  ACR-145_TP2_R1_001.fastq.gz.md5
  ACR-145_TP2_R2_001.fastq.gz.md5
  ACR-150_TP2_R1_001.fastq.gz.md5
  ACR-150_TP2_R2_001.fastq.gz.md5
  ACR-173_TP2_R1_001.fastq.gz.md5
  ACR-173_TP2_R2_001.fastq.gz.md5
  ACR-178_TP2_R1_001.fastq.gz.md5
  ACR-178_TP2_R2_001.fastq.gz.md5

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

Trim ACRP-CON reads.

Reads for samples ACR-140, ACR-145, ACR-150, ACR-173, and ACR-178 were trimmed using the built-in version of Trimmomatic with the default settings, following the 9FastQ QC and Trimming - E5 Coral RNA-seq Data for A.pulchra protocol)[https://robertslab.github.io/sams-notebook/2023/05/19/FastQ-QC-and-Trimming-E5-Coral-RNA-seq-Data-for-A.pulchra-P.evermanni-and-P.meandrina-Using-FastQC-fastp-and-MultiQC-on-Mox.html].

Downloaded all files to Andromeda URI HPC location

  ```
cd /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed

wget -r -nd -A .fastq.gz https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/A_pulchra/trimmed/

RNA-ACR-140-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
RNA-ACR-140-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
RNA-ACR-145-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
RNA-ACR-145-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
RNA-ACR-150-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
RNA-ACR-150-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
RNA-ACR-173-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
RNA-ACR-173-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
RNA-ACR-178-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
RNA-ACR-178-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz

  ```

# 1) Run FastQC

a) Make folders for raw FastQC results and scripts

```
cd /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/

mkdir fastqc_results

```

b) Write script for checking quality with FastQC and submit as job on Andromeda

```
nano /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/fastqc_raw.sh
```

```  
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=danielle_becker@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/raw
#SBATCH --error="script_error" #if your job fails, the error report will be put in this file
#SBATCH --output="output_script" #once your job is completed, any final job report comments will be put in this file

module load FastQC/0.11.9-Java-11

for file in /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/raw/*ACRP
do
fastqc $file --outdir /data/putnamlab/dbecks/DeNovo_transcriptxome/2023_A.pul/data/fastqc_results/
done
```

```
sbatch /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/fastqc_raw.sh

Submitted batch job 281440
```

## Combined QC output into 1 file with MultiQC, do not need a script due to fast computational time

```
module load MultiQC/1.9-intel-2020a-Python-3.8.2

multiqc /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/fastqc_results/*fastqc.zip -o /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/fastqc_results/multiqc/

```

c) Copy MultiQC and FastQC files to local computer

```
scp -r danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/fastqc_results/multiqc/multiqc_report.html /Users/Danielle/Desktop/Putnam_Lab/Gametogenesis/bioinformatics/transcriptome/original_fastqc

scp -r danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/fastqc_results/*.html /Users/Danielle/Desktop/Putnam_Lab/Gametogenesis/bioinformatics/transcriptome/original_fastqc

```

# 4) Trim and clean reads

a) Make trimmed reads folder in all other results folders

```
mkdir data/trimmed
cd trimmed

```

c) Write script for Trimming and run on Andromeda

#Run fastp on files
#Trims 20bp from 5' end of all reads
#Trims poly G, if present

```
nano /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/trim.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=danielle_becker@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/raw
#SBATCH --error="script_error" #if your job fails, the error report will be put in this file
#SBATCH --output="output_script" #once your job is completed, any final job report comments will be put in this file

module load fastp/0.19.7-foss-2018b

array1=($(ls *R1*.fastq.gz)) #Make an array of sequences to trim
for i in ${array1[@]}; do
fastp --in1 ${i} --in2 $(echo ${i}|sed s/_R1/_R2/) --detect_adapter_for_pe --trim_poly_g --trim_front1 20 --trim_front2 20 --out1 ../trimmed/${i} --out2 ../trimmed/$(echo ${i}|sed s/_R1/_R2/)  
done

```
```
sbatch /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/trim.sh


```


# 5) Check quality of trimmed files

a) Check number of files in /trimmed directory

```
ls -1 | wc -l
#64
```


b) Check number of reads in /trimmed directory

```
zgrep -c "@GWNJ" *.gz > trimmed_seq_counts

```

c) Run FastQC on trimmed data
```
mkdir trimmed_qc

```
```
nano /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/fastqc_trimmed.sh
```

```  
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trimmed_qc
#SBATCH --error="script_error"
#SBATCH --output="output_script"

module load FastQC/0.11.8-Java-1.8

for file in /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/*.gz
do
fastqc $file --outdir /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trimmed_qc
done
```

```
sbatch /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/fastqc_trimmed.sh



```


d) Run MultiQC on trimmed data, Combined QC output into 1 file with MultiQC, do not need a script due to fast computational time

```
module load MultiQC/1.9-intel-2020a-Python-3.8.2

multiqc /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trimmed_qc/*fastqc.zip -o /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trimmed_qc/trimmed_multiqc
```

e) Copy multiqc and fastqc to computer, use terminal window fro desktop not in server

```
scp -r danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trimmed_qc/trimmed_multiqc/multiqc_report.html /Users/Danielle/Desktop/Putnam_Lab/Gametogenesis/bioinformatics/transcriptome/trimmed_fastqc

scp -r danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trimmed_qc/*.html /Users/Danielle/Desktop/Putnam_Lab/Gametogenesis/bioinformatics/transcriptome/trimmed_fastqc

```

Run Trinity with forward and reverse sequences

  ```

#!/bin/bash
#SBATCH --job-name=20230923_trinity
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=danielle_becker@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed
#SBATCH --error="script_error" #if your job fails, the error report will be put in this file
#SBATCH --output="output_script" #once your job is completed, any final job report comments will be put in this file


  # Load Trinity module

  module load Trinity/2.15.1-foss-2022a


  # Document programs in PATH (primarily for program version ID)

  date >> system_path.log
  echo "" >> system_path.log
  echo "System PATH for $SLURM_JOB_ID" >> system_path.log
  echo "" >> system_path.log
  printf "%0.s-" {1..10} >> system_path.log
  echo ${PATH} | tr : \\n >> system_path.log


  # Run Trinity
  /gscratch/srlab/programs/trinityrnaseq-Trinity-v2.8.3/Trinity \
  --trimmomatic \
  --seqType fq \
  --max_memory 500G \
  --CPU 28 \
  --left \
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-15_S7_L004_R1_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-15_S7_L004_R1_0348.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-16_S8_L004_R1_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-16_S8_L004_R1_0348.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-17_S9_L004_R1_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-17_S9_L004_R1_0348.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-18_S10_L004_R1_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-18_S10_L004_R1_0348.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-1_S1_L004_R1_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-1_S1_L004_R1_0348.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-2_S2_L004_R1_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-2_S2_L004_R1_0348.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-3_S3_L004_R1_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-3_S3_L004_R1_0348.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-4_S4_L004_R1_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-4_S4_L004_R1_0348.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-4Spl_S11_L004_R1_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-4Spl_S11_L004_R1_0348.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-5_S5_L004_R1_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-5_S5_L004_R1_0348.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-6_S6_L004_R1_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-6_S6_L004_R1_0348.fastq.gz \
  --right \
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-15_S7_L004_R2_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-15_S7_L004_R2_0348.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-16_S8_L004_R2_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-16_S8_L004_R2_0348.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-17_S9_L004_R2_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-17_S9_L004_R2_0348.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-18_S10_L004_R2_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-18_S10_L004_R2_0348.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-1_S1_L004_R2_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-1_S1_L004_R2_0348.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-2_S2_L004_R2_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-2_S2_L004_R2_0348.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-3_S3_L004_R2_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-3_S3_L004_R2_0348.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-4_S4_L004_R2_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-4_S4_L004_R2_0348.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-4Spl_S11_L004_R2_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-4Spl_S11_L004_R2_0348.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-5_S5_L004_R2_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-5_S5_L004_R2_0348.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-6_S6_L004_R2_0343.fastq.gz,\
  /gscratch/scrubbed/samwhite/data/O_lurida/RNAseq/CP-6_S6_L004_R2_0348.fastq.gz
  ```
