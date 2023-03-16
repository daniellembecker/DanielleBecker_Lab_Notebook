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