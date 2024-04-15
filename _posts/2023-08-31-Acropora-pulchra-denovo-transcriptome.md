---
layout: post
title: Workflow for Acropora pulchra de novo transcriptome
category: [ de novo transcriptome , DNA]
tag: [ Acropora pulchra, de novo transcriptome ]
---
## Designing a workflow to create a de novo transcriptome for *Acropora pulchra*

#### Goal:
Use samples from 11 *Acropora pulchra* colonies collected in Mo'orea, French Polynesia on January 15th 2022 from the north shore backreef site Mahana (17°29'13.9"S 149°53'14.7"W) part of a 12-month [Gametogenesis timeseries project](https://github.com/daniellembecker/Gametogenesis) that were than [concentrated into one sequenced sample](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2023-04-25-Acropora-pulchra-transcriptome-extraction-concentration.md) using the [Zymo RNA Clean Concentrate Protocol](https://github.com/zdellaert/ZD_Putnam_Lab_Notebook/blob/master/protocols/Zymo_RNA_Clean_Concentrate.pdf) and five sequence samples also collected from Mo'orea, French Polynesia part of the [E5 Rules of Life project](https://github.com/urol-e5) to create a de novo transcriptome for *A. pulchra*. Literature review of current *Acropora* de novo transcriptomes and genomes completed already.

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

Reads for samples ACR-140, ACR-145, ACR-150, ACR-173, and ACR-178 were trimmed using the built-in version of Trimmomatic with the default settings, following the 9FastQ QC and Trimming - [E5 Coral RNA-seq Data for A.pulchra protocol](https://robertslab.github.io/sams-notebook/posts/2023/2023-05-19-FastQ-QC-and-Trimming---E5-Coral-RNA-seq-Data-for-A.pulchra-P.evermanni-and-P.meandrina-Using-FastQC-fastp-and-MultiQC-on-Mox/). All fast

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

# 1) Run FastQC on ACRP_R1 and ACRP_R2

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

array1=($(ls *ACRP_R1*.fastq.gz)) #Make an array of sequences to trim
for i in ${array1[@]}; do
fastp --in1 ${i} --in2 $(echo ${i}|sed s/_R1/_R2/) --detect_adapter_for_pe --trim_poly_g --trim_front1 20 --trim_front2 20 --out1 ../trimmed/${i} --out2 ../trimmed/$(echo ${i}|sed s/_R1/_R2/)  
done

```
```
sbatch /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/trim.sh
Submitted batch job 281915

```


# 5) Check quality of trimmed files

a) Check number of files in /trimmed directory

```
ls -1 | wc -l
#12
```

b) Check number of reads in /trimmed directory and in /raw directory

```
# look at raw reads

zgrep -c "@A01587" ACRP* > seq_counts

ACRP_R1_001.fastq.gz:178057412
ACRP_R2_001.fastq.gz:178057412

# look at trimmed reads

zgrep -c "@A01587" ACRP* > trimmed_seq_counts

ACRP_R1_001.fastq.gz:173194930
ACRP_R2_001.fastq.gz:173194930

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

for file in /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/ACRP*
do
fastqc $file --outdir /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trimmed_qc
done
```

```
sbatch /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/fastqc_trimmed.sh

Submitted batch job 281933

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

# 6) Run Trinity with forward and reverse sequences

Trinity commands as found on [Trinity vignette](https://www.broadinstitute.org/videos/introduction-de-novo-rna-seq-assembly-using-trinity) and [Roberts Lab Guide](https://robertslab.github.io/resources/bio_Transcriptome-assembly/):
 -    Trinity  = runs trinity program
 -    seqType fq  = specifys FASTQ file format
 -    SS_lib_type RF  = Trinity has the option (--SS_lib_type) to specify whether or not the sequences you're assembly are      "stranded". User should specify this in the following fashion as on option in the Trinity command (example specifies typical stranded libraries): --SS_lib_type RF
 -    max_memory 100G  = max memory for trinity 100G should not be changed, per communications with the developer
 -    CPU 36  = match however many CPUs your computing cluster has access to, maximum number of parallel processes (for andromeda the maximum number of cores is 36.  On Unity, the nodes in the uri-cpu partition can go up to 64, and in the cpu and cpu-preempt partitions it can go to 128)
 -    left  = paired end reads filenames
 -    right = paired end reads filenames


a) Write Trinity script

 ```
 nano /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/trinity.sh
 ```

```

#!/bin/bash
#SBATCH --job-name=20230925_trinity
#SBATCH --time=30-00:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --exclusive
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

# Run Trinity

Trinity \
--seqType fq \
--SS_lib_type RF \
--max_memory 100G \
--CPU 36 \
--left \
/data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/ACRP_R1_001.fastq.gz,\
/data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/RNA-ACR-140-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz,\
/data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/RNA-ACR-145-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz,\
/data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/RNA-ACR-150-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz,\
/data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/RNA-ACR-173-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz,\
/data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/RNA-ACR-178-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz \
--right \
/data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/ACRP_R2_001.fastq.gz,\
/data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/RNA-ACR-140-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz,\
/data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/RNA-ACR-145-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz,\
/data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/RNA-ACR-150-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz,\
/data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/RNA-ACR-173-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz,\
/data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/RNA-ACR-178-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz

  ```

  b) Run trinity script

  ```
  sbatch /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/trinity.sh

Submitted batch job 282369 on 20230927 at 15:56

  ```

Trinity finished 20231011 with a FailedCommands.txt file. I submitted a GitHub issues with my output script to the [Trinity GitHub Issues](https://github.com/trinityrnaseq/trinityrnaseq/issues/1342#issuecomment-1760111445).

```
Response to issue:

It looks like there's 18 failed commands.  I'm seeing the errors like this:

/glfs/brick01/gv0/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trinity_out_dir/read_partitions/Fb_8/CBin_8455/c845648.trinity.reads.fa.out/inchworm.fa.tmp':


If you find the failed_commands.txt file (or similar name) in the trinity
output directory, it should have the list of the commands it's having a
hard time with.   Each should have an output directory, like in the above
example, it would be:

/glfs/brick01/gv0/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trinity_out_dir/read_partitions/Fb_8/CBin_8455/c845648.trinity.reads.fa.out

If you can remove those output directories and then rerun your Trinity
command, it should hopefully process everything just fine (unless you ran
out of disk space...)

```

I then removed all the corrupt files (18) listed in the FailedCommands.txt file and re-ran Trinity on 20231012 at 16:04 after changing the amount of CPUs to 20 because a number of solutions in the developer's GitHub Issues suggest that sometimes the CPUs process things too quickly for the computer's file system to keep up, leading to problems

```
rm -rf /glfs/brick01/gv0/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trinity_out_dir/read_partitions/Fb_8/CBin_8455/c845648.trinity.reads.fa.out

rm -rf /glfs/brick01/gv0/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trinity_out_dir/read_partitions/Fb_6/CBin_6008/c600976.trinity.reads.fa.out

rm -rf /glfs/brick01/gv0/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trinity_out_dir/read_partitions/Fb_8/CBin_8295/c829640.trinity.reads.fa.out

rm -rf /glfs/brick01/gv0/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trinity_out_dir/read_partitions/Fb_7/CBin_7785/c778615.trinity.reads.fa.out

rm -rf /glfs/brick01/gv0/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trinity_out_dir/read_partitions/Fb_0/CBin_139/c13942.trinity.reads.fa.out

rm -rf /glfs/brick01/gv0/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trinity_out_dir/read_partitions/Fb_7/CBin_7922/c792325.trinity.reads.fa.out

rm -rf /glfs/brick01/gv0/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trinity_out_dir/read_partitions/Fb_4/CBin_4006/c400692.trinity.reads.fa.out

rm -rf /glfs/brick01/gv0/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trinity_out_dir/read_partitions/Fb_3/CBin_3881/c388169.trinity.reads.fa.out

rm -rf /glfs/brick01/gv0/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trinity_out_dir/read_partitions/Fb_8/CBin_8455/c845648.trinity.reads.fa.out

rm -rf /glfs/brick01/gv0/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trinity_out_dir/read_partitions/Fb_5/CBin_5650/c565150.trinity.reads.fa.out

rm -rf /glfs/brick01/gv0/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trinity_out_dir/read_partitions/Fb_5/CBin_5326/c532719.trinity.reads.fa.out

rm -rf /glfs/brick01/gv0/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trinity_out_dir/read_partitions/Fb_6/CBin_6986/c698736.trinity.reads.fa.out

rm -rf /glfs/brick01/gv0/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trinity_out_dir/read_partitions/Fb_6/CBin_6598/c659988.trinity.reads.fa.out

rm -rf /glfs/brick01/gv0/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trinity_out_dir/read_partitions/Fb_6/CBin_6003/c600479.trinity.reads.fa.out

rm -rf /glfs/brick01/gv0/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trinity_out_dir/read_partitions/Fb_4/CBin_4526/c452745.trinity.reads.fa.out

rm -rf /glfs/brick01/gv0/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trinity_out_dir/read_partitions/Fb_2/CBin_2558/c255895.trinity.reads.fa.out

rm -rf /glfs/brick01/gv0/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trinity_out_dir/read_partitions/Fb_2/CBin_2424/c242423.trinity.reads.fa.out

rm -rf /glfs/brick01/gv0/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trimmed/trinity_out_dir/read_partitions/Fb_7/CBin_7230/c723106.trinity.reads.fa.out

```

```
sbatch /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/trinity.sh

Submitted batch job 285355
```

Trinity completed successfully!!!! When Trinity completes, it creates a 'trinity_out_dir.Trinity.fasta' output file (or prefix based on the output directory you specify).

Trinity groups transcripts into clusters based on shared sequence content. Such a transcript cluster is very loosely referred to as a 'gene'. This information is encoded in the Trinity fasta accession. An example Fasta entry for one of the transcripts is formatted like so:

```
>TRINITY_DN1000_c115_g5_i1 len=247 path=[31015:0-148 23018:149-246]
 AATCTTTTTTGGTATTGGCAGTACTGTGCTCTGGGTAGTGATTAGGGCAAAAGAAGACAC
 ACAATAAAGAACCAGGTGTTAGACGTCAGCAAGTCAAGGCCTTGGTTCTCAGCAGACAGA
 AGACAGCCCTTCTCAATCCTCATCCCTTCCCTGAACAGACATGTCTTCTGCAAGCTTCTC
 CAAGTCAGTTGTTCACAGGAACATCATCAGAATAAATTTGAAATTATGATTAGTATCTGA
 TAAAGCA
 ```

Download Trinity fasta to Desktop if needed, too large to have stored there always

 ```
scp -r danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trinity_out_dir/trinity_out_dir.Trinity.fasta /Users/Danielle/Desktop/Putnam_Lab/Gametogenesis/bioinformatics

  ```

# 7) Assessing quality of the assessment

Use [Trinity toolkit utilities](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats) for a assembly quality assessment

a) Run TrinityStats.pl script for stats output, need to use the path of trinity module and then add the util folder before you can access the .pl script

This script will compute the contig Nx statistics (eg. the contig N50 value), in addition to a modification of the Nx statistic that takes into consideration transcript expression (read support) data, which we call the ExN50 statistic.

Based on the lengths of the assembled transcriptome contigs, we can compute the conventional Nx length statistic, such that at least x% of the assembled transcript nucleotides are found in contigs that are at least of Nx length. The traditional method is computing N50, such that at least half of all assembled bases are in transcript contigs of at least the N50 length value.

```
/opt/software/Trinity/2.15.1-foss-2022a/trinityrnaseq-v2.15.1/util/TrinityStats.pl trinity_out_dir.Trinity.fasta > trinity_assembly_stats

```

```
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':  1019918
Total trinity transcripts:      1476390
Percent GC: 44.78

########################################
Stats based on ALL transcript contigs:
########################################

        Contig N10: 4395
        Contig N20: 2790
        Contig N30: 1956
        Contig N40: 1417
        Contig N50: 1021

        Median contig length: 376
        Average contig: 675.99
        Total assembled bases: 998026397

```

The N10 through N50 values are shown computed based on all assembled contigs. In this example, 10% of the assembled bases are found in transcript contigs at least 4,395 bases in length (N10 value), and the N50 value indicates that at least half the assembled bases are found in contigs that are at least 1,021 bases in length.

The contig N50 values can often be exaggerated due to an assembly program generating too many transcript isoforms, especially for the longer transcripts. To mitigate this effect, the script will also compute the Nx values based on using only the single longest isoform per 'gene':

```
#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

        Contig N10: 3840
        Contig N20: 2221
        Contig N30: 1479
        Contig N40: 1026
        Contig N50: 731

        Median contig length: 335
        Average contig: 565.18
        Total assembled bases: 576438524
```

You can see that the Nx values based on the single longest isoform per gene are lower than the Nx stats based on all assembled contigs, as expected, and even though the Nx statistic is really not a reliable indicator of the quality of a transcriptome assembly, the Nx value based on using the longest isoform per gene is perhaps better for reasons described above.


b) Run BUSCO ( Benchmarking Universal Single-Copy Orthologs)

- *Commands and overview for running BUSCO here*: https://busco.ezlab.org/busco_userguide.html
- uses highly conserved single-copy orthologs; evolutionary informed expectations of gene content.
- appears that youu can focus a BUSCO analysis to orthologs related to your target taxa.

- image below shows a BUSCO analysis comparing the crayfish targetted for tde novo transcriptome assembly to 44 other arthropod species assemblies and a single vertebrate assembly:

*Theissinger et al. 2016* https://www.sciencedirect.com/science/article/abs/pii/S1874778716300137

Citation: Theissinger, K., Falckenhayn, C., Blande, D., Toljamo, A., Gutekunst, J., Makkonen, J., ... & Kokko, H. (2016). De Novo assembly and annotation of the freshwater crayfish Astacus astacus transcriptome. Marine Genomics, 28, 7-10.


Make run-busco-transcriptome.sh script:

```
nano /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/busco.sh

```

```
#!/bin/bash

#SBATCH --job-name="busco"
#SBATCH --time="100:00:00"
#SBATCH --nodes 1 --ntasks-per-node=20
#SBATCH --mem=250G
##SBATCH --output="busco-%u-%x-%j"
##SBATCH --account=putnamlab
##SBATCH --export=NONE

echo "START" $(date)

labbase=/data/putnamlab
busco_shared="${labbase}/shared/busco"
[ -z "$query" ] && query="${labbase}/dbecks/DeNovo_transcriptome/2023_A.pul/data/trinity_out_dir/trinity_out_dir.Trinity.fasta" # set this to the query (genome/transcriptome) you are running
[ -z "$db_to_compare" ] && db_to_compare="${busco_shared}/downloads/lineages/metazoa_odb10"

source "${busco_shared}/scripts/busco_init.sh"  # sets up the modules required for this in the right order

# This will generate output under your $HOME/busco_output
cd "${labbase}/${dbecks}"
busco --config "${busco_shared}/scripts/busco-config.ini"  -f -c 20 --long -i "${query}" -l "${db_to_compare}" -o busco_output -m transcriptome

echo "STOP" $(date)

```

Run BUSCO script

```

sbatch /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/busco.sh


```

My first run of BUSCO, I was getting this repeated errors in my slurm output error file:

"Message: BatchFatalError(AttributeError("'NoneType' object has no attribute 'hmmer_results_lines'"))"

The issue was that the above script was referencing busco-config.ini which had old versions (based on the old modules) instead of using the installed copy $EBROOTBUSCO/config/config.ini which has the correct paths. So I changed the line in the script below to use the new correct path and it began to run!

There was also another error that the metazoa_odb10 file was an old version and to make sure BUSCO always updates these file versions, we need to use the command -l metazoa_odb10 instead of a set file path for BUSCO to automatically update the versions. Further explanation could be found in the [download and automatically update lineaages BUSCO guide section](https://busco.ezlab.org/busco_userguide.html#download-and-automated-update).


```
#!/bin/bash

#SBATCH --job-name="busco"
#SBATCH --time="100:00:00"
#SBATCH --nodes 1 --ntasks-per-node=20
#SBATCH --mem=250G
##SBATCH --output="busco-%u-%x-%j"
##SBATCH --account=putnamlab
##SBATCH --export=NONE

echo "START" $(date)

labbase=/data/putnamlab
busco_shared="${labbase}/shared/busco"
[ -z "$query" ] && query="${labbase}/dbecks/DeNovo_transcriptome/2023_A.pul/data/trinity_out_dir/trinity_out_dir.Trinity.fasta" # set this to the query (genome/transcriptome) you are running
[ -z "$db_to_compare" ] && db_to_compare="${busco_shared}/downloads/lineages/metazoa_odb10"

source "${busco_shared}/scripts/busco_init.sh"  # sets up the modules required for this in the right order

# This will generate output under your $HOME/busco_output
cd "${labbase}/${dbecks}"
busco --config "$EBROOTBUSCO/config/config.ini"  -f -c 20 --long -i "${query}" -l metazoa_odb10 -o busco_output -m transcriptome

echo "STOP" $(date)

```

```

sbatch /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/busco.sh

Submitted batch job 288779

```



### BUSCO Results

```
cd /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/busco/busco_output
less short_summary.specific.metazoa_odb10.busco_output.txt

# BUSCO version is: 5.2.2
# The lineage dataset is: metazoa_odb10 (Creation date: 2021-02-17, number of genomes: 65, number of BUSCOs: 954)
# Summarized benchmarking in BUSCO notation for file /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trinity_out_dir/trinity_out_dir.Trinity.fasta
# BUSCO was run in mode: transcriptome

        ***** Results: *****

        C:99.9%[S:21.0%,D:78.9%],F:0.1%,M:0.0%,n:954       
        953     Complete BUSCOs (C)                        
        200     Complete and single-copy BUSCOs (S)        
        753     Complete and duplicated BUSCOs (D)         
        1       Fragmented BUSCOs (F)                      
        0       Missing BUSCOs (M)                         
        954     Total BUSCO groups searched               

```

BUSCO completeness looks great. Completeness looks for the presence or absence of highly conserved genes in an assembly. The aim is to have the highest percentage of genes identified in your assembly, with a BUSCO complete score above 95% considered good (we have 99.9% so yay!).

However, the complete and duplicated BUSCOs are high. Transcriptomes and protein sets that are not filtered for isoforms will lead to a high proportion of duplicates. So, next step suggestions:

1. Filter for isoforms
2. Map to closest genome (*Acropora millepora*)
3. Filter symbiont genes to check if it helps duplication

# 8) Map to *Acropora millepora* genome

#### a) Obtain Reference Genome and Create Folder Structure

I am using the Maldives *Acropora millepora* genome to compare alignment statistics.

[Ying, Hua, et al. "The whole-genome sequence of the coral Acropora millepora." Genome biology and evolution 11.5 (2019): 1374-1379.](https://academic.oup.com/gbe/article/11/5/1374/5480306)

Location on Andromeda, the HPC server for URI:
```
cd /data/putnamlab/REFS/Amil_ref
```

All of the genome files (i.e., genome scaffolds, CDS, proteins, GFF, and functional annotations) were downloaded from the [NCBI RefBase](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=45264) and added to the Andromeda HPC folder.

```
scp -r /Users/Danielle/Downloads/ncbi_dataset/ncbi_dataset/data/GCF_013753865.1/GCF_013753865.1_Amil_v2.1_genomic.fna danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/REFS/Amil_ref/
```


#### b) Generate genome build

##### HiSat2 Align reads to refernece genome
[HiSat2](https://daehwankimlab.github.io/hisat2/main/)
[HiSat2 Github](https://github.com/DaehwanKimLab/hisat2)


```
nano /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/Hisat2_genome_build.sh
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
#SBATCH -D /data/putnamlab/REFS/Amil_ref/
#SBATCH --error="script_error"
#SBATCH --output="output_script"

module load HISAT2/2.2.1-gompi-2022a

# index the reference genome for Amil output index to working directory
hisat2-build -f /data/putnamlab/REFS/Amil_ref/GCF_013753865.1_Amil_v2.1_genomic.fna ./Amil_ref # called the reference genome (scaffolds)
echo "Reference genome indexed. Starting alignment" $(date)
```

```
sbatch /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/Hisat2_genome_build.sh

Submitted batch job 289152
```

#### c) Align reads to genome

```
mkdir data/mapped
```

- Trinity assembly in FASTA format and you want to align it to a reference genome, you would generally use a different tool like BLAT, GMAP, or minimap2, as they are designed for aligning longer sequences like transcripts or assembled contigs.
- Using minimap2
  - In this example, -ax splice specifies that the input is spliced transcripts, -uf specifies the reference genome, and the output is redirected to a SAM file (Trinity_aligned.sam).

```
nano /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/minimap_align2.sh
```

```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/ref_genome_Amil/mapped
#SBATCH --cpus-per-task=3
#SBATCH --error="script_error"
#SBATCH --output="output_script"


module load minimap2/2.24-GCCcore-11.3.0

#Aligning to Trinity output file

minimap2 -d /data/putnamlab/REFS/Amil_ref/reference_index.idx /data/putnamlab/REFS/Amil_ref/GCF_013753865.1_Amil_v2.1_genomic.fna

minimap2 -ax splice -uf /data/putnamlab/REFS/Amil_ref/reference_index.idx /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/trinity_out_dir/trinity_out_dir.Trinity.fasta > trinity_aligned.sam

```

```
sbatch /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/minimap_align2.sh

Submitted batch job 290621

```

#### d) Sort and convert sam to bam and check number of mapped reads and mapping percentages

- Explanation:
  - samtools sort -o Trinity_aligned.sorted.bam Trinity_aligned.sam: This command sorts the SAM file and creates a sorted BAM file (Trinity_aligned.sorted.bam).
  - samtools index Trinity_aligned.sorted.bam: This command creates an index for the sorted BAM file. The index file (Trinity_aligned.sorted.bam.bai) is necessary for certain operations and viewers.
  - samtools flagstat Trinity_aligned.sorted.bam: This command generates statistics about the alignment, including the number of mapped reads and mapping percentages.

```
nano /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/SAMtoBAM.sh

#There will be lots of .tmp file versions in your folder, this is normal while this script runs and they should delete at the end to make one sorted.bam file
```
```
#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --export=NONE
#SBATCH --mem=500GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/ref_genome_Amil/mapped
#SBATCH --cpus-per-task=3
#SBATCH --error="script_error"
#SBATCH --output="output_script"

module load SAMtools/1.16.1-GCC-11.3.0 #Preparation of alignment for assembly: SAMtools

samtools view -bS trinity_aligned.sam > trinity_aligned.bam
samtools sort trinity_aligned.bam -o trinity_aligned_sorted.bam
samtools index trinity_aligned_sorted.bam
samtools flagstat trinity_aligned_sorted.bam > alignment_stats.txt

```

```
sbatch /data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/scripts/SAMtoBAM.sh

Submitted batch job 290625

```

**Alignment statistics**

```

2337097 + 0 in total (QC-passed reads + QC-failed reads)
1476390 + 0 primary
805750 + 0 secondary
54957 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
1714368 + 0 mapped (73.35% : N/A)
853661 + 0 primary mapped (57.82% : N/A)

Explanation:

1. Total Reads:
  - 2337097 + 0 in total (QC-passed reads + QC-failed reads): Indicates the total number of reads, including both QC-passed and QC-failed reads. In this case, there are 2,337,097 reads in total.

2. Primary and Secondary Alignments:
  - 1476390 + 0 primary: The number of primary alignments. These are the primary alignment records for each read. In this case, there are 1,476,390 primary alignments.
  - 805750 + 0 secondary: The number of secondary alignments.
  - Secondary alignments can occur for reads that map equally well to multiple locations in the reference genome.

3. Supplementary Alignments:
  - 54957 + 0 supplementary: The number of supplementary alignments. - Supplementary alignments are used to represent chimeric or novel splice junctions.

4. Duplicates:
  - 0 + 0 duplicates: The number of duplicate reads. Duplicate reads can result from PCR artifacts and are often removed in quality control.

5. Primary Duplicates:
  - 0 + 0 primary duplicates: The number of duplicate primary reads. This specifically refers to duplicate primary alignment records.

6. Mapped Reads:
  - 1714368 + 0 mapped (73.35% : N/A): The total number and percentage of mapped reads. In this case, 1,714,368 reads are mapped, and they constitute 73.35% of the total reads.

7. Primary Mapped Reads:
  - 853661 + 0 primary mapped (57.82% : N/A)
  - The number and percentage of primary mapped reads.
  - Primary mapped reads refer to those reads where the primary alignment is reported.
  - In this case, 853,661 reads are primary mapped, constituting 57.82% of the total reads.

```


#### e) Download mapping percentages and statistics to desktop

```

scp -r danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/DeNovo_transcriptome/2023_A.pul/data/ref_genome_Amil/mapped/alignment_stats.txt /Users/Danielle/Desktop/Putnam_Lab/Gametogenesis/bioinformatics/transcriptome

```

# 9) Identify and analyze isoforms in **de novo** traanscriptome

#### a) Extract isoform sequences and gene-isoform mapping to visualize isoform positions per gene and identify "stacked" or mutually exclusive isoform expression patterns

```
nano /data/putnamlab/dbecks/Heatwave_A.pul_2022Project/scripts/gene_isoform_map.sh

```
```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=32GB
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.ed
#SBATCH --account=putnamlab
#SBATCH --output=slurm-%A.out
#SBATCH -D /data/putnamlab/dbecks/Heatwave_A.pul_2022Project/data/trinity/trinity_mapped/all_samples/

module load Trinity/2.15.1-foss-2022a
module load SAMtools/1.16.1-GCC-11.3.0
module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2


## Generate isoform expression estimates
$EBROOTTRINITY/trinityrnaseq-v2.15.1/util/abundance_estimates_to_matrix.pl \
    --est_method RSEM \
    --out_prefix isoform_expr \
    --gene_trans_map /data/putnamlab/dbecks/Heatwave_A.pul_2022Project/data/trinity/trinity_mapped/all_samples/gene_to_isoform.map \
    --name_sample_file /data/putnamlab/dbecks/Heatwave_A.pul_2022Project/data/trinity/trinity_mapped/rsem_sample_names.txt \
    --name_sample_by_basedir \
    --quant_files /data/putnamlab/dbecks/Heatwave_A.pul_2022Project/data/trinity/trinity_mapped/rsem_isoform_files.txt

```

```
sbatch /data/putnamlab/dbecks/Heatwave_A.pul_2022Project/scripts/gene_isoform_map.sh

Submitted batch job 312338
```

#### b) Download isoform_expr.isoform.TPM.not_cross_norm and gene_to_isoform.map files for processing in RStudio

- The isoform_expr.isoform.TPM.not_cross_norm file contains Transcripts Per Million (TPM) values, which are length-normalized expression values. TPM values are more interpretable and easier to compare across isoforms within a gene, as they account for isoform length differences.

```

scp -r danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Heatwave_A.pul_2022Project/data/trinity/trinity_mapped/all_samples/isoform_expr.isoform.TPM.not_cross_norm /Users/Danielle/Desktop/

scp -r danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Heatwave_A.pul_2022Project/data/trinity/trinity_mapped/all_samples/gene_to_isoform.map /Users/Danielle/Desktop/

```
