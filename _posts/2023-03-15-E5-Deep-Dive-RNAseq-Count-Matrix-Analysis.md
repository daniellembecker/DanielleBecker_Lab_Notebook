---
layout: post
title: E5 Deep Dive RNAseq Count Matrix Analysis
date: '2023-03-15'
categories: E5
tags: Bioinformatics GeneExpression
---
This post details analysis and bioinformatic steps to generate a gene count matrix from RNAseq data for the E5 Deep Dive Project. 

More information on this project can be found on the [GitHub repo](https://github.com/urol-e5/deep-dive). 

# 1. Obtain RNAseq files from NCBI SRA 

## Obtain SRR numbers  

First, we needed to obtain RNAseq files from NCBI from the project [PRJNA744403](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA744403). I first obtained a list of SRR numbers that identify the specific sequence files for RNAseq. I did this by going to [all SRA links for the project](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=744403) and selecting all the RNAseq files (32 total). In the upper right hand corner, I selected "Send to" and "File". This outputs a .txt file with a list of all SRR numbers. 

## Download RNAseq files  

I then logged into Andromeda and created a folder for these files at `/data/putnamlab/ashuffmyer/e5-deepdive` with `raw`, `scripts`, and `refs` folders. Jill Ashey is also working on this project. 

I then copied the SRR identifiers from the .txt file produced above into a file in the `raw` folder using `nano SraAccList.txt` and pasting the identifiers. The file looks like this:  

```
SRR15101688
SRR15101689
SRR15101690
SRR15101691
SRR15101692
SRR15101693
SRR15101694
SRR15101695
SRR15101696
SRR15101697
SRR15101699
SRR15101700
SRR15101701
SRR15101702
SRR15101703
SRR15101704
SRR15101705
SRR15101706
SRR15101707
SRR15101708
SRR15101710
SRR15101711
SRR15101712
SRR15101713
SRR15101715
SRR15101718
SRR15101719
SRR15101720
SRR15101721
SRR15101722
SRR15101723
SRR15101724
```

Next, I downloaded the files using the [NCBI SRA Toolkit](https://github.com/ncbi/sra-tools) which is installed on Andromeda as module `SRA-Toolkit/2.10.9-gompi-2020b`.  

With the help of [Sam and Steven in the Roberts Lab](https://github.com/RobertsLab/resources/issues/1569), I wrote a script to download the files. 

`nano download_sra.sh` 

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="download_sra_error" #if your job fails, the error report will be put in this file
#SBATCH --output="download_sra_output" #once your job is completed, any final job report comments will be put in this file

module load SRA-Toolkit/2.10.9-gompi-2020b

prefetch --option-file ../raw/SraAccList.txt -O ../raw #this creates a folder for each SRR in the .txt list and outputs in the raw data folder  

# Enable recursive globbing
shopt -s globstar

# Run fasterq-dump on any SRA file in any directory and split into read 1 and 2 files and put in raw folder 
for file in ../raw/**/*.sra
do
  fasterq-dump --outdir ../raw --split-files "${file}"
done

#Remove the SRR directories that are no longer needed
rm -r ../raw/SRR*/

#zip all files 
gzip ../raw/*.fastq

#Generate checksums 
md5sum  ../raw/*.fastq.gz > ../raw/md5.original.download.20230315
```

`sbatch download_sra.sh`  

This script uses `prefetch` to retrieve the data files for each SRR number and stores them in the `raw` data folder with a directory for each file. 

The `fasterq-dump` command then takes the `.sra` file from each directory and converts them to read 1 and read 2 fastq files. 

Finally, the SRR directories are removed that are no longer needed and we generate md5 checksums.   

This script can be used for any SRA download with a custom list of SRR numbers. 

## Download reference files 

The genome is stored on [NCBI](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCA_014529365.1/). 

We downloaded the reference files from the [Reef Genomics Database](http://pver.reefgenomics.org/). 

```
cd refs

#genome scaffolds
wget http://pver.reefgenomics.org/download/Pver_genome_assembly_v1.0.fasta.gz

#gene models CDS
wget http://pver.reefgenomics.org/download/Pver_genes_names_v1.0.fna.gz

#gene models proteins
wget http://pver.reefgenomics.org/download/Pver_proteins_names_v1.0.faa.gz

#gene models GFF
wget http://pver.reefgenomics.org/download/Pver_genome_assembly_v1.0.gff3.gz

#Generate checksums 
md5sum  *.gz > md5.original.refs.download.20230315
```

Finally, I updated permissions after all files were downloaded so Jill can collaborate on this project.  

`chmod u=rwx,g=rwx,o=rwx,a=rwx -R e5-deepdive`  
 
Additional steps will be added to this post as we move forward. 

# 2. QC using FastQC and MultiQC 

Make directories for fastqc data 

`mkdir fastqc_raw`

Write script for raw QC

`nano fastqc_raw.sh`

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/scripts              
#SBATCH --error="fastqc_raw_error" #if your job fails, the error report will be put in this file
#SBATCH --output="fastqc_raw_output" #once your job is completed, any final job report comments will be put in this file

# Load modules needed 
module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.9-intel-2020a-Python-3.8.2

echo "Start raw QC" $(date)

cd /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq

for file in /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/*fastq.gz
do 
fastqc $file --outdir /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/fastqc_raw
done

echo "End raw QC" $(date)

# Compile MultiQC report from fastQC files 
cd /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/fastqc_raw
multiqc --interactive ./

echo "Cleaned MultiQC report generated." $(date)

```

Submitted batch job 243914. Script ran in ~4.5 hours.  


## Raw sequence MultiQC summary: 

The raw MultiQC report [can be found on the E5 Deep Dive GitHub repo here](https://github.com/urol-e5/deep-dive/blob/main/A-Pver/data/rna-seq/raw_multiqc_report.html). 

There was high adapter content in these raw sequences. Overall, the quality score was high and there were no red flags in other metrics.

![fastqc_sequence_counts_plot](https://user-images.githubusercontent.com/32178010/227280396-e9fac75c-8357-451d-af84-6112455a3618.png)  
![fastqc_per_base_sequence_quality_plot](https://user-images.githubusercontent.com/32178010/227280538-5fb9bccf-4551-46d3-92f3-dbf5078c9dfa.png)  
![fastqc_per_sequence_quality_scores_plot](https://user-images.githubusercontent.com/32178010/227280588-f238d863-086b-4dcd-9df6-fdc6db386f82.png)  
![fastqc_per_sequence_gc_content_plot](https://user-images.githubusercontent.com/32178010/227280634-384b948a-64c6-441d-9ca2-acbadd60e9ec.png)  
![fastqc_per_base_n_content_plot](https://user-images.githubusercontent.com/32178010/227280664-8a297ea2-4310-4a5b-b349-598517cb6e00.png)  
![fastqc_sequence_duplication_levels_plot](https://user-images.githubusercontent.com/32178010/227280709-f36ec7c2-d63c-46b5-b754-993e744006d5.png)  
![fastqc_overrepresented_sequencesi_plot](https://user-images.githubusercontent.com/32178010/227280742-53d5d8ca-c216-4598-af5d-38fdfc38a243.png)  
![fastqc_adapter_content_plot](https://user-images.githubusercontent.com/32178010/227280785-11c27c3f-d4f8-451e-a9a5-85e3b5021945.png)  

# 3. Trimming sequences with FastP

Make a script to trim sequences with FastP. We will be removing adapters, trimming poly-g, and trimming by quality phred score with a cut off of 30. 

`nano fastp.sh`

```
#!/bin/bash

#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/scripts              
#SBATCH --error="fastp_error" #if your job fails, the error report will be put in this file
#SBATCH --output="fastp_output" #once your job is completed, any final job report comments will be put in this file

# Load modules needed 
module load fastp/0.19.7-foss-2018b

echo "Start trimming with fastp" $(date)

# Make array of sequences to trim 
cd /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw
array1=($(ls *1.fastq.gz))

# fastq and fastqc loop
for i in ${array1[@]}; do
    fastp --in1 ${i} \
        --in2 $(echo ${i}|sed s/_1/_2/)\
        --out1 /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/trim/trim.${i} \
        --out2 /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/trim/trim.$(echo ${i}|sed s/_1/_2/) \        
        --detect_adapter_for_pe \
        --qualified_quality_phred 30 \
        --trim_poly_g 
done

echo "Read trimming of adapters complete." $(date)
```

Submitted batch job 244006. Took about 8 hrs to run 


# 4. QC trimmed sequences 

Now we will QC the trimmed sequences. Make a script. 

`nano fastqc_trim.sh`

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/scripts              
#SBATCH --error="fastqc_trim_error" #if your job fails, the error report will be put in this file
#SBATCH --output="fastqc_trim_output" #once your job is completed, any final job report comments will be put in this file

# Load modules needed 
module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.9-intel-2020a-Python-3.8.2

echo "Start trimmed QC" $(date)

cd /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq

for file in /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/trim/*fastq.gz
do 
fastqc $file --outdir /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/fastqc_trim
done

echo "End trimmed QC" $(date)

# Compile MultiQC report from fastQC files 
cd /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/fastqc_trim
multiqc --interactive ./

echo "Trimmed MultiQC report generated." $(date)
```

Submitted batch job 244018.

## Trimmed sequence MultiQC summary: 

The trimmed MultiQC report [can be found on the E5 Deep Dive GitHub repo here](https://github.com/urol-e5/deep-dive/blob/main/A-Pver/data/rna-seq/trim_multiqc_report.html). 

Trimming removed the adapter content in these sequences. The quality scores are very high across the sequence length, so we are not going to trim by length in this iteration. We can return to the trimming step in the future if we need to make changes.  

![fastqc_adapter_content_plot](https://user-images.githubusercontent.com/32178010/227281183-0d57f77f-dbd0-4004-84ed-27eae50f53f3.png)
![fastqc_sequence_duplication_levels_plot](https://user-images.githubusercontent.com/32178010/227281201-1625d008-0423-4eb2-b531-ba34c321914f.png)
![fastqc_sequence_length_distribution_plot](https://user-images.githubusercontent.com/32178010/227281214-c92857af-effc-4207-887e-f28059b178f8.png)
![fastqc_per_base_n_content_plot](https://user-images.githubusercontent.com/32178010/227281227-fa12fce3-c81f-4a62-888f-4b290601263c.png)
![fastqc_per_sequence_gc_content_plot](https://user-images.githubusercontent.com/32178010/227281242-b3c6dc36-4d8d-45a3-bef5-83203b077c30.png)
![fastqc_per_sequence_quality_scores_plot](https://user-images.githubusercontent.com/32178010/227281249-d675d5d9-c973-4882-bcc1-52953b2d030c.png)
![fastqc_per_base_sequence_quality_plot](https://user-images.githubusercontent.com/32178010/227281259-08c420fe-2da8-4293-bb8c-b3ff48580798.png)
![fastqc_sequence_counts_plot](https://user-images.githubusercontent.com/32178010/227281270-f3f72a60-e487-4c4c-b0be-8f994fecf15a.png)

## Calculate sequence length 

We also ran a script to calculate sequence length of raw and trimmed sequences.

Make a script. 

`nano length.sh`

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/scripts              
#SBATCH --error="seq_length_error" #if your job fails, the error report will be put in this file
#SBATCH --output="seq_length_output" #once your job is completed, any final job report comments will be put in this file

echo "Start sequence length counts" $(date)

cd /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq

zgrep -c "@SRR" /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/*.gz > /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/raw_seq_counts.txt

zgrep -c "@SRR" /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/trim/*.gz > /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/trim/trim_seq_counts.txt

echo "End sequence length counts" $(date)
```
Submitted batch job 244020.

The results are below. There are about 20-25 million sequences in each file.    

```
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101688_1.fastq.gz:23396439
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101688_2.fastq.gz:23396439
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101689_1.fastq.gz:24990687
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101689_2.fastq.gz:24990687
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101690_1.fastq.gz:23267238
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101690_2.fastq.gz:23267238
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101691_1.fastq.gz:29417106
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101691_2.fastq.gz:29417106
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101692_1.fastq.gz:22578178
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101692_2.fastq.gz:22578178
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101693_1.fastq.gz:20905338
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101693_2.fastq.gz:20905338
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101694_1.fastq.gz:22990550
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101694_2.fastq.gz:22990550
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101695_1.fastq.gz:16484060
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101695_2.fastq.gz:16484060
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101696_1.fastq.gz:24030431
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101696_2.fastq.gz:24030431
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101697_1.fastq.gz:24846067
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101697_2.fastq.gz:24846067
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101699_1.fastq.gz:22733345
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101699_2.fastq.gz:22733345
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101700_1.fastq.gz:25158606
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101700_2.fastq.gz:25158606
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101701_1.fastq.gz:29523059
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101701_2.fastq.gz:29523059
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101702_1.fastq.gz:23796710
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101702_2.fastq.gz:23796710
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101703_1.fastq.gz:22630044
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101703_2.fastq.gz:22630044
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101704_1.fastq.gz:24455094
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101704_2.fastq.gz:24455094
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101705_1.fastq.gz:27805087
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101705_2.fastq.gz:27805087
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101706_1.fastq.gz:18568153
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101706_2.fastq.gz:18568153
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101707_1.fastq.gz:29131006
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101707_2.fastq.gz:29131006
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101708_1.fastq.gz:25068848
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101708_2.fastq.gz:25068848
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101710_1.fastq.gz:23675842
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101710_2.fastq.gz:23675842
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101711_1.fastq.gz:24020166
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101711_2.fastq.gz:24020166
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101712_1.fastq.gz:24937261
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101712_2.fastq.gz:24937261
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101713_1.fastq.gz:16381619
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101713_2.fastq.gz:16381619
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101715_1.fastq.gz:21751368
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101715_2.fastq.gz:21751368
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101718_1.fastq.gz:27076800
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101718_2.fastq.gz:27076800
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101719_1.fastq.gz:24642131
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101719_2.fastq.gz:24642131
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101720_1.fastq.gz:26361873
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101720_2.fastq.gz:26361873
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101721_1.fastq.gz:17900262
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101721_2.fastq.gz:17900262
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101722_1.fastq.gz:25641209
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101722_2.fastq.gz:25641209
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101723_1.fastq.gz:23770610
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101723_2.fastq.gz:23770610
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101724_1.fastq.gz:25368402
/data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/raw/SRR15101724_2.fastq.gz:25368402
```


# 5. Align sequences to the reference genome 

We are using `HISAT2` and `samtools` to build the reference genome and align our sequences to this reference. See Step 1 above for reference genome information. 

## Build reference genome  

Create a script. 

`nano build.sh`

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=jillashey@uri.edu #your email to send notifications
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/scripts              
#SBATCH --error="build_error" #if your job fails, the error report will be put in this file
#SBATCH --output="build_output" #once your job is completed, any final job report comments will be put in this file

# load modules needed
module load HISAT2/2.2.1-gompi-2021b #Alignment to reference genome: HISAT2

# Unzip reference genome 
#gunzip /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/refs/Pver_genome_assembly_v1.0.fasta.gz

# Index reference genome 
hisat2-build -f /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/refs/Pver_genome_assembly_v1.0.fasta Pver_ref
echo "Reference genome indexed." $(date)
```

This built the reference genome we need to do the alignment. 

## Align sequences  

Create script. 

`nano align.sh`

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab  
#SBATCH -D /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/trim/
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications        
#SBATCH --error="align_error" #if your job fails, the error report will be put in this file
#SBATCH --output="align_output" #once your job is completed, any final job report comments will be put in this file

# load modules needed
module load HISAT2/2.2.1-foss-2019b #Alignment to reference genome: HISAT2
module load SAMtools/1.9-foss-2018b #Preparation of alignment for assembly: SAMtools

## Genome already referenced 

echo "Start alignment" $(date)

# Alignment of clean reads to the reference genome
# Make array of sequences to trim 
array1=($(ls /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/trim/*_1.fastq.gz | xargs -n 1 basename))

for i in ${array1[@]}; do
    hisat2 -p 8 --rna-strandness RF --dta -q -x /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/refs/Pver_ref -1 /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/trim/${i} -2 /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/trim/$(echo ${i}|sed s/_1/_2/) -S /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/mapped/${i}.sam
    samtools sort -@ 8 -o /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/mapped/${i}.bam /data/putnamlab/ashuffmyer/e5-deepdive/rna-seq/mapped/${i}.sam
    echo "${i} bam-ified!"
    rm ${i}.sam
done

echo "Alignment completed" $(date)

```

This script will align the sequences to the references, generate `.bam` and `.sam` files, and then remove the unneeded and large `.sam` files. 

This was started on 20230323 at 09:45 Pacific Time job 244260. 
