---
layout: post
title: Mcapitata Development 16S Analysis Part 3
date: '2021-12-30'
categories: Mcapitata_EarlyLifeHistory_2020
tags: 16S Mcapitata Molecular Protocol
---
This post details QC and QIIME analysis for the 16S analysis adapted from the pipeline [developed by Emma Strand](https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2021-06-21-16s-Analysis-Pipeline.md). 

# 16S QC and Analysis in QIIME2   

[Steps #6-7 in Part 2 post](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/Mcapitata-Development-16S-Analysis-Part-2/)  

In previous post, I created metadata, prepared the environment, and loaded data into a QIIME artifact.  

Next we can proceed with cleaning and QC in QIIME following [E Strand analysis post](https://emmastrand.github.io/EmmaStrand_Notebook/E5-16S-Analysis/).  

### 8. Denoising with DADA2  

```
ssh -l ashuffmyer ssh3.hac.uri.edu
cd /data/putnamlab/ashuffmyer/AH_MCAP_16S

cd scripts
nano denoise.sh

```

We are using the following parameters based on [MultiQC report of 16S data](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/Mcapitata-Development-16S-Analysis-Part-1/).   

P-trim values are based on primer length:  

- p-trim-left-r 20 (reverse is 20 bp long)
- p-trim-left-f 19 (forward is 19 bp long)

Primer information for our sequencing run [is here](https://emmastrand.github.io/EmmaStrand_Notebook/16s-Sequencing-HoloInt/).  

P-truncate values are based on where the mean quality scores of R1 and R2 files start to decrease seen in the MultiQC report:    

- p-trunc-len-r 192
- p-trunc-len-f 245

Forward sequences: 
![forward](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/16S/forwardquality.png)  

Reverse sequences:  
![reverse](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/16S/reversequality.png)  

Create a script to run denoising step with these parameters.  

`mkdir processed_data`  

Move data produced in step 7 to the raw_data folder.  

`mv AH-MCAP-16S-paired-end-sequences1.qza raw_data/AH-MCAP-16S-paired-end-sequences1.qza`

Other commands used in the script are as follows:  

- `--i-demultiplexed-seqs` followed by the sequences artifact to be denoised
- `--p-trunc-len-f INTEGER`: position to be truncated due to decreased quality. This truncates the 3' end of sequences which are the bases that were sequenced in the last cycles. On the forward read.
- `--p-trunc-len-r INTEGER`: same as above but on the reverse read.
- `p-trim-left-f INTEGER`: Position at which forward read sequences should be trimmed due to low quality. This trims the 5' end of the input sequences, which will be the bases that were sequenced in the first cycles.
- `p-trim-left-r INTEGER`: Position at which reverse read sequences should be trimmed due to low quality. This trims the 5' end of the input sequences, which will be the bases that were sequenced in the first cycles.
- `o-table`: The resulting feature table.
- `o-representative-sequences`: The resulting feature sequences. Each feature in the feature table will be represented by exactly one sequence, and these sequences will be the joined paired-end sequences.
- `o-denoising-stats`: SampleData[DADA2Stats]
- `p-n-threads`: The number of threads to use for multithreaded processing. If 0 is provided, all available cores will be used.

This script will also generate qiime metadata table, feature-table summarized and feature-table tabulated. This will be output in the `processed_data` directory and the script lives in the `script` directory.   

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                 
#SBATCH --error="dada_script_error" #if your job fails, the error report will be put in this file
#SBATCH --output="dada_output_script" #once your job is completed, any final job report comments will be put in this file

source /usr/share/Modules/init/sh # load the module function
module load QIIME2/2021.8

cd ../processed_data

# Metadata path
METADATA="../metadata/sample_metadata.csv"

qiime dada2 denoise-paired --verbose --i-demultiplexed-seqs ../raw_data/AH-MCAP-16S-paired-end-sequences1.qza \
  --p-trunc-len-r 192 --p-trunc-len-f 245 \
  --p-trim-left-r 20 --p-trim-left-f 19 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza \
  --p-n-threads 20

# Summarize feature table and sequences
  qiime metadata tabulate \
    --m-input-file denoising-stats.qza \
    --o-visualization denoising-stats.qzv
  qiime feature-table summarize \
    --i-table table.qza \
    --o-visualization table.qzv \
    --m-sample-metadata-file $METADATA
  qiime feature-table tabulate-seqs \
    --i-data rep-seqs.qza \
    --o-visualization rep-seqs.qzv

```

