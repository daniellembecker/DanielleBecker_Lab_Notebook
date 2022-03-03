---
layout: post
title: 16S Pipeline in Mothur
date: '2022-01-12'
categories: Analysis Mcapitata_EarlyLifeHistory_2020
tags: 16S Mcapitata Molecular Protocol R
---
This post details 16S data analysis using the Mothur pipeline for the 2020 *Montipora capitata* developmental time series.  

# **Pipeline: 16S analysis in Mothur using HPC**  

For more information on Mothur, visit the [Mothur website](https://mothur.org/).  

I previously did preliminary analysis on our 16S data using the QIIME pipeline. See posts [Part 1 here](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/Mcapitata-Development-16S-Analysis-Part-1/), [Part 2 here](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/Mcapitata-Development-16S-Analysis-Part-2/), and [Part 3 here](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/Mcapitata-Development-16S-Analysis-Part-3/).  

Mothur is written in C++. Mothur takes in sequence data, cleans the data, and outputs files that can be used for taxonomy and diversity analyses. Mothur is independent of operating system and requires no dependencies.  

Information on this pipeline and workflow can be viewed in the [Mothur MiSeq SOP wiki](https://mothur.org/wiki/miseq_sop/).

The associated R project for this data can be found on [GitHub here](https://github.com/AHuffmyer/EarlyLifeHistory_Energetics).  

This pipeline details 16S analysis of MiSeq sequencing data of the bacterial V4 16S region.   

General Workflow:  
1. [Prepare directory](#Directory)   
2. [Start mothur](#Start)    
3. [Preparing sequences](#Prepare)   
4. [QC sequences](#QC)    
5. [Unique sequences](#Unique)    
6. [Aligning](#Align)    
7. [Preclustering](#Precluster)    
8. [Identify chimeras](#Chimera)  
9. [Classify sequences](#Classify)     
10. [Cluster OTU's](#Cluster)    
11. [Subsampling](#Subsample)  
12. [Calculate ecological statistics](#Statistics)  
13. [Output data for R analysis](#Output)   

## <a name="Directory"></a> **1. Prepare Directory**  

First, prepare a directory that will contain all mothur analyses. This is located within the `AH_MCAP_16S` directory in the URI Andromeda HPCC.  

Copy files into the `mothur` directory.  

```
ssh -l ashuffmyer ssh3.hac.uri.edu
cd /data/putnamlab/ashuffmyer/AH_MCAP_16S
mkdir mothur

cp /data/putnamlab/ashuffmyer/AH_MCAP_16S/raw_data/*.gz /data/putnamlab/ashuffmyer/AH_MCAP_16S/mothur
```

The current mothur module available is: Mothur/1.46.1-foss-2020b

As of 21 January 2022, this is the current version of Mothur available.  

## <a name="Start"></a> **2. Start Mothur**  

Start in interactive mode to see how we can get mothur to run and test that the module is working.  

```
interactive
module load Mothur/1.46.1-foss-2020b
cd mothur/
mothur

```

This successfully activated mothur. This should see a display of citation, program information, and `mothur >`.  

At any time, use `quit` to leave mothur.  

At the top of the header when starting mothur, you will see version and last updated as well as citation and other information.  

`mothur >` is the "prompt". 

You can use interactive mode to run computationally small commands. 

*In this file, I show each command that will be used and then in each section I show how these commands are run in one bash script. If you want, you can test run each step and then combine all steps into a single bash file to run the job one at a time. However, I recommend running each step separately so that you can view and use information from each step to inform the next analysis.*      

## <a name="Prepare"></a> **3. Preparing Sequences: make.file, make.contig, and summary.seq**  

First, we need to prepare files and sequences to be analyzed.  

#### make.file()  

make.file() tells mothur to look for fastq files in the current directory (`mothur`) and identify forward and reverse reads. Put type=gz to look for gz files. If there are .fasta or .fastq files you can use type=fasta or type=fastq. The prefix gives the output file a name - this can be a project name. Here we will use "mcap". 

*For running a different project, all you have to do is change the prefix for all the files in all the commands in this document.*    

The code will look like this:    
 
```
make.file(inputdir=., type=gz, prefix=mcap)
```

This creates a file called `mcap.files` with sample name, R1, R2 listed in columns.  

Then we will align sequences together for R1 and R2 to generate contigs. 

#### make.contigs() 

We will write a bash script to run the `make.contigs()` command and make contigs of forward and reverse reads for each sample. This will assemble contigs for each pair of files and use R1 and R2 quality information to correct any sequencing errors in locations where there is overlap and the quality of the complementary sequence is sufficient. 

We also need to create an "oligos.txt" file in a text editor that has the primer information so that the primer sequences can be removed. Our primers are 515F and 806RB for V4 region. The file for this data has one line and looks like the below. If you have barcodes, you can look at the MiSeq SOP (link above) to view formatting requirements for the oligos file.    

Primer information for these samples is [here](https://emmastrand.github.io/EmmaStrand_Notebook/16s-Sequencing-HoloInt/).  

```
nano oligos.oligos 

primer GTGCCAGCMGCCGCGGTAA GGACTACNVGGGTWTCTAAT
```

In the make.contigs step, we need to use trimoverlap=T because we used 2x300 bp sequences that is longer than our amplicon length. 

The code will look like this:  

```
make.contigs(inputdir=., outputdir=., file=mcap.files, oligos=oligos.oligos, trimoverlap=T)
```

#### summary.seqs()  

We will also run the `summary.seqs()` command with the `make.file()` and `make.contig()` steps. This generates summary information about the sequences from the files we made above. The code will look like this: 

```
summary.seqs(fasta=mcap.trim.contigs.fasta)
```

#### Run all of these commands in a script       

Within the `mothur` directory, create a script. Everything will live within this directory - I did not make any subfolders. 

```
nano contigs.sh
```

Here is the script with the commands all together.  

```
#!/bin/bash
#SBATCH --job-name="contigs"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications                
#SBATCH --error="contigs_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="contigs_output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -q putnamlab

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#make.file(inputdir=., type=gz, prefix=mcap)"

mothur "#make.contigs(inputdir=., outputdir=., file=mcap.files, trimoverlap=T)"

mothur "#summary.seqs(fasta=mcap.trim.contigs.fasta)"

```

Note the syntax for running mothur commands in a slurm script. The syntax is to activate mothur with `mothur` followed by the command in quotes and a hashtag (`"#make.files()"`).  

Run script and check status.    

```
sbatch contigs.sh

squeue -u ashuffmyer -l
```

To view the output of each script, you can view with the following:  

```
nano contigs_output_script
```

*For all of these steps in mothur, you will be looking at these script output files to see any output and lists of any files that were created. I recommend copying these outputs so that you can track the output files.*  

For 39 samples, this took about 30 minutes.  

The following files were output from the `make.files()` step: 
 
```
mcap.files
```

The following files were output from the `make.contigs()`step:  

```
mcap.trim.contigs.fasta
mcap.scrap.contigs.fasta
mcap.contigs.report
mcap.contigs.groups
```

Descriptions of contig files:  
*Trim file* = sequences that were "good".  
*Scrap file* = sequences that were "bad".  
*Groups file* = what group each sequence belongs to map sequence to each sample from the trimmed sequence file.  
*Contigs report file* = information on sequences that were aligned and paired together. 

Count the number of sequences that were removed and the number that were kept by counting sequences in each fasta file.  

```
grep -c "^>" mcap.trim.contigs.fasta
grep -c "^>" mcap.scrap.contigs.fasta
```

I repeated the repeat contigs step with a few settings to try to troubleshoot the loss of sequences at this step. I suspect we are not treating our primers appropriately. Here are the results: 

317,884 sequences were kept and 1,194,315 were removed.  *with trimoverlap=T*    
317,884 sequences were kept and 1,194,315 were removed.  *with trimoverlap=T+checkorient=T*   
1,131,991 sequences were kept and 337,769 were removed.  *with trimoverlap=T+checkorient=T+pdiffs=2*   

```
#make.contigs(inputdir=., outputdir=., file=mcap.files, oligos=oligos.oligos, trimoverlap=T, checkorient=T, pdiffs=2)"
```

1,501,977 sequences were kept and 10,222 were removed.  *with trimoverlap=T+no oligo/primer file* 

```
#make.contigs(inputdir=., outputdir=., file=mcap.files, trimoverlap=T)"
```

View the scrap file to see what error codes are included for the sequences that were removed.  

```
head -100 mcap.scrap.contigs.fasta
```

*We get very different amounts of sequences kept depending on if we use an oligos file and what our pdiff value is. I keep way more sequences if I have pdiffs=2 and even more without an oligos file. But I think the oligos file is for when we need to demultiplex - which we dont.*  

View the kept sequences and search (command+F) for the primer sequences: 

```
head -1000 mcap.trim.contigs.fasta
```

The primers we are looking for are: 

`F GTGCCAGCMGCCGCGGTAA R GGACTACNVGGGTWTCTAAT`

When I search for the primer sequences in this file they do not show up. Grep to count for any total number of primer sequences.  

```
grep -c "GTGCCAGCMGCCGCGGTAA" mcap.trim.contigs.fasta
grep -c "GGACTACNVGGGTWTCTAAT" mcap.trim.contigs.fasta

grep -c "GTGCCAGCMGCCGCGGTAA" mcap.scrap.contigs.fasta
grep -c "GGACTACNVGGGTWTCTAAT" mcap.scrap.contigs.fasta
```
The primers show up 0 times in the kept contigs and scrap sequences. The primers are not present. Therefore, we do not need an oligos file.  

Results below are from running `make.contigs` with `trimoverlap=T` and no oligos file.  

View the contigs report.  

```
less mcap.contigs.report
```

This shows a long and detailed file of contig assembly information. 

If you have barcodes, you can use the `oligos` command to direct to a file that has the primer and barcode information during the contig step to demultiplex your sequences. [See more information here](https://mothur.org/wiki/oligos_file/#).  

When you run summary.seqs() you will get an output table like the one below and it will output a file with this summary information.  

```
nano contigs_output_script

```

Scroll to the bottom to find the summary information.  


``` 
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	1	1	0	1	1
2.5%-tile:	1	35	35	0	3	37550
25%-tile:	1	247     247     0	5	375495
Median:         1	247     247     0	7	750989
75%-tile:	1	247     247     1	7	1126483
97.5%-tile:     1	292     292     35	35	1464428
Maximum:        1	305     305     137     35	1501977
Mean:   1	214     214     6	10
# of Seqs:	1501977
``` 


This table shows quantile values about the distribution of sequences for a few things:  

- *Start position*: All at 1 now, will start at different point after some QC.   
- *End position*: We see that there are some sequences that are very short and we may need to remove those later.   
- *Number of bases*: length (we see most are in expected range here, but one is super long! This might tell us there is no overlap so they are butted up against each other. We will remove things like this.  
- *Ambigs*: Number of ambiguous calls in sequences. Here there are a few that have ambiguous base calls. We will remove any sequence with an ambiguous call or any longer than we would expect for V4 region.   
- *Polymer*: Length of polymer repeats.    
- *NumSeqs*: Number of sequences.  

### <a name="QC"></a> **4. QC'ing sequences with screen.seqs**    

Now, based on the output above, we need to remove "bad" sequences from the dataset. In the `screen.seqs()` function, we will specify the fasta file of the contigs generated in the previous step and remove any sequence with an ambiguous call ("N"). We will also remove sequences >350 nt. We will also set a minimum size amount (200). These parameters could be adjusted based on specific experiment and variable region.  

The code will look like this:  

```
screen.seqs(fasta=mcap.trim.contigs.fasta, group=mcap.contigs.groups, maxambig=0, maxlength=300, minlength=200)
```

Note that when making output files, Mothur keeps adding tags to indicate different pipeline steps in the file name. 

Generate a script to run `screen.seqs()` and re run `summary.seqs()` to get information about the sequences.    

```
nano screen.sh
```

```
#!/bin/bash
#SBATCH --job-name="screen"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="screen_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="screen_output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -q putnamlab

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#screen.seqs(inputdir=., outputdir=., fasta=mcap.trim.contigs.fasta, group=mcap.contigs.groups, maxambig=0, maxlength=350, minlength=200)"

mothur "#summary.seqs(fasta=mcap.trim.contigs.good.fasta)"

```

```
sbatch screen.sh
```

This generates the following output files: 
`mcap.trim.contigs.good.fasta`
`mcap.trim.contigs.bad.accnos`
`mcap.contigs.good.groups`

"good" seqs satisfied criteria and "bad" seqs did not meet criteria. 

Count the number of "good" and "bad" sequences.  

```
grep -c "^>" mcap.trim.contigs.good.fasta
grep -c ".*" mcap.trim.contigs.bad.accnos
```

950,129 sequences were kept and 551,848 were removed. 

This removed ~35% of sequences due to length and ambiguous bases.  

You can also view the `bad.accnos` file to see why sequences were removed.  

```
head -1000 mcap.trim.contigs.bad.accnos
``` 

The summary output as viewed by `nano screen_output_script` file now reads: 

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	200     200     0	3	1
2.5%-tile:	1	247     247     0	4	23754
25%-tile:	1	247     247     0	5	237533
Median:         1	247     247     0	7	475065
75%-tile:	1	292     292     0	7	712597
97.5%-tile:     1	292     292     0	7	926376
Maximum:        1	300     300     0	19	950129
Mean:   1	258     258     0	6
# of Seqs:	950129

```

We now see that we have removed all sequences with ambigous calls and the max sequence length is <300 with min >200. In steps below we will align to a reference V4 and filter again.    

### <a name="Unique"></a> **5. Determining and counting unique sequences**  

Next, determine the number of unique sequences. The code will look like this: 

```
unique.seqs(fasta=mcap.trim.contigs.good.fasta)
```

The unique.seqs() step will output the following files: 
`mcap.trim.contigs.good.names`  
`mcap.trim.contigs.good.unique.fasta`  

After determining the unique sequences, we can use these unique sequences to count how many times each sequence (which will later be classified to OTU) shows up in each sample.  

```
count.seqs(name=mcap.trim.contigs.good.names, group=mcap.contigs.good.groups) 
```

`count.seqs()` outputs a file named: 
`mcap.trim.contigs.good.count_table`  

Finally, we can run `summary.seq()` again to view the "clean" data.  

Now run summary seqs again.  

```
summary.seqs(fasta=mcap.trim.contigs.good.unique.fasta, count=mcap.trim.contigs.good.count_table)
```  

During the run (visualize at the end of the `unique_output_script` file), the first column is number of sequences looked at, and the second column is the number of unique sequences. 

Generate and run a script.  

```
nano unique.sh
```

```
#!/bin/bash
#SBATCH --job-name="unique"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="unique_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="unique_output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -q putnamlab

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#unique.seqs(fasta=mcap.trim.contigs.good.fasta)"

mothur "#count.seqs(name=mcap.trim.contigs.good.names, group=mcap.contigs.good.groups)"

mothur "#summary.seqs(fasta=mcap.trim.contigs.good.unique.fasta, count=mcap.trim.contigs.good.count_table)"

```


```
sbatch unique.sh
```

This script generate the following outputs: 

```
mcap.trim.contigs.good.names
mcap.trim.contigs.good.unique.fasta
mcap.trim.contigs.good.count_table
```

```
nano unique_output_script
``` 

In this run, there were 950,129 sequences and 109,985 were unique = ~11%.  

The output table from summary.seqs looks like this and shows the number of unique and the total: 

```

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	200     200     0	3	1
2.5%-tile:	1	247     247     0	4	23754
25%-tile:	1	247     247     0	5	237533
Median:         1	247     247     0	7	475065
75%-tile:	1	292     292     0	7	712597
97.5%-tile:     1	292     292     0	7	926376
Maximum:        1	300     300     0	19	950129
Mean:   1	258     258     0	6
# of unique seqs:	109985
total # of seqs:        950129

```

Now we can align just the unique sequences, which will be much faster than aligning the full data set and is an indicator of how polished and clean the data are.  

*From this, we have our unique sequences identified and can proceed with further cleaning and polishing of the data. Next we will look at alignment, error rate, chimeras, classification and further analysis.*  

### <a name="Align"></a> **6. Aligning to reference database**

#### Prepare the reference sequences  

First, download the Silva reference files from the [Mothur Wiki](https://mothur.org/wiki/silva_reference_files/) at the latest release.

The silva reference is used and recommended by the Mothur team. It is a manually curated data base with high diversity and high alignment quality.  

In your mothur directory, run `wget` to download the silva reference and training sets from the mothur website and unzip them and move into the right directory. 

```
wget https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.bacteria.zip

wget https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset9_032012.pds.zip

unzip silva.bacteria.zip
cp silva.bacteria.fasta ../silva.bacteria.fasta

unzip trainset9_032012.pds.zip
```

Now we have the reference files that we need in our directory.  

Now we are going to take our Silva database reference alignment and select the V4 region.  

We can do this with the `pcr.seqs` function. 

```
pcr.seqs(fasta=silva.bacteria.fasta, start=11894, end=25319, keepdots=F)
```

This region includes the primers, even though our sequences dont include the primers. If we did keepdots=T then the first hundreds of columns would be periods, which is not useful for us. These periods are placeholders in the silva database.   

Write a script to do this and generate a summary of the reference and rename the file to something more useful to us.   

The commands will be: 

```
summary.seqs(fasta=silva.bacteria.pcr.fasta)

rename.file(input=silva.bacteria.pcr.fasta, new=silva.v4.fasta)
```

```
nano silva_ref.sh
```

```
#!/bin/bash
#SBATCH --job-name="silva_ref"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="silva_ref_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="silva_ref_output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -q putnamlab

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#pcr.seqs(fasta=silva.bacteria.fasta, start=11894, end=25319, keepdots=F)"

mothur "#summary.seqs(fasta=silva.bacteria.pcr.fasta)"

mothur "#rename.file(input=silva.bacteria.pcr.fasta, new=silva.v4.fasta)"

```

```
sbatch silva_ref.sh
```

The script outputs the following files: 

```
Output File Names: 
silva.bacteria.pcr.fasta
```

The file was renamed to `silva.v4.fasta` as specified in the script. 

The summary of sequences looks like this in the script output: 

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	13424   270     0	3	1
2.5%-tile:	1	13425   292     0	4	374
25%-tile:	1	13425   293     0	4	3740
Median:         1	13425   293     0	4	7479
75%-tile:	1	13425   293     0	5	11218
97.5%-tile:     1	13425   294     1	6	14583
Maximum:        3	13425   351     5	9	14956
Mean:   1	13424   292     0	4
# of Seqs:	14956

```

We now have a reference to align to.  

#### Align sequences to the reference  

We next align sequences with the `align.seqs` command.   

```
align.seqs(fasta=mcap.trim.contigs.good.unique.fasta, reference=silva.v4.fasta)

summary.seqs(fasta= mcap.trim.contigs.good.unique.align)
```

Write a script to do this, generate a new summary of the output, and run.  

```
nano align.sh
```

```
#!/bin/bash
#SBATCH --job-name="align"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="align_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="align_output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -q putnamlab

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#align.seqs(fasta=mcap.trim.contigs.good.unique.fasta, reference=silva.v4.fasta)"

mothur "#summary.seqs(fasta=mcap.trim.contigs.good.unique.align)"

```

```
sbatch align.sh
``` 

The script generates these output files: 

```
Output File Names:
mcap.trim.contigs.good.unique.align
mcap.trim.contigs.good.unique.align.report
mcap.trim.contigs.good.unique.flip.accnos
```

View the report file.  

```
head mcap.trim.contigs.good.unique.align.report
```

View the accnos file.  

```
head mcap.trim.contigs.good.unique.flip.accnos
```

The summary now looks like this:  

```
nano align_output_script 
```


```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	3	1	0	1	1
2.5%-tile:	1	1258    8	0	2	2750
25%-tile:	1	13424   22	0	3	27497
Median:         1	13424   292     0	4	54993
75%-tile:	1	13424   292     0	5	82489
97.5%-tile:     13398   13425   293     0	6	107236
Maximum:        13425   13425   300     0	13	109985
Mean:   1769    11409   206     0	3
# of Seqs:	109985

```

From this, we see that 109,985 sequences aligned to the reference, which matches the number of sequences that we had after the unique.sh step (109,985). In the next steps we will filter out any sequences that don't meet alignment settings.  

#### QC sequences according to alignment to the reference  

Our sequences now align at the correct positions on the reference.

Now remove sequences that are outside the alignment window (1968-11550bp). This removes anything that starts after `start` and ends before `end`. Maxhomop=8 argument removes anything that has repeats greater than the threshold - e.g., 8 A's in a row = polymer 8. Here we will removes polymers >8 because we are confident these are likely not high quality sequences (see mothur MiSeq SOP for more information).  

This is the command we will use: 

```
screen.seqs(fasta=mcap.trim.contigs.good.unique.align, count=mcap.trim.contigs.good.count_table, start=1968, end=11550, maxhomop=8)
```

We will then run a summary.  

```
summary.seqs(fasta=mcap.trim.contigs.good.unique.good.align, count=mcap.trim.contigs.good.good.count_table)
```

Write a script and run this screening step.  

```
nano screen2.sh
```

```
#!/bin/bash
#SBATCH --job-name="screen2"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="screen2_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="screen2_output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -q putnamlab

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#screen.seqs(fasta=mcap.trim.contigs.good.unique.align, count=mcap.trim.contigs.good.count_table, start=1968, end=11550, maxhomop=8)"

mothur "#summary.seqs(fasta=mcap.trim.contigs.good.unique.good.align, count=mcap.trim.contigs.good.good.count_table)"

```

```
sbatch screen2.sh
```


This will output the following files: 
```
Output File Names:
mcap.trim.contigs.good.unique.good.align
mcap.trim.contigs.good.unique.bad.accnos
mcap.trim.contigs.good.good.count_table
```

View the accnos file to see why sequences will be removed and count the number of "bad" sequences.  

```
head mcap.trim.contigs.good.unique.bad.accnos

grep -c ".*" mcap.trim.contigs.good.unique.bad.accnos
``` 

34,704 uniques are tagged to be removed due to filtering at this step. This totals to 695,380 total sequences removed (this number is in the output file). 

*We therefore have a large portion (~73%) of our sequences that do not align to the bacterial V4 region and therefore may be from the host.*     

The summary looks like this: 

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	11552   274     0	3	1
2.5%-tile:	1	13424   292     0	4	6369
25%-tile:	1	13424   292     0	4	63688
Median:         1	13424   292     0	4	127375
75%-tile:	1	13424   292     0	5	191062
97.5%-tile:     1	13424   292     0	6	248381
Maximum:        1253    13425   300     0	8	254749
Mean:   2	13423   291     0	4
# of unique seqs:	75281
total # of seqs:        254749
```

#### Filter sequences  

Now we can filter out sequences that didn't meet our criteria above, which will generate a report and a new summary of our sequences.  

We will run the following code. We align vertically and use trump=. to align the sequences accounting for periods in the reference.   

```
filter.seqs(fasta=mcap.trim.contigs.good.unique.good.align, vertical=T, trump=.)
```

Write and run a script.  

```
nano filter.sh
```

```
#!/bin/bash
#SBATCH --job-name="filter"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="filter_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="filter_output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -q putnamlab

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#filter.seqs(fasta=mcap.trim.contigs.good.unique.good.align, vertical=T, trump=.)"

mothur "#summary.seqs(fasta=mcap.trim.contigs.good.unique.good.filter.fasta, count=mcap.trim.contigs.good.good.count_table)"

```

```
sbatch filter.sh
```

This script outputs the following files:  

```
Output File Names: 
mcap.filter
mcap.trim.contigs.good.unique.good.filter.fasta
```

We get a report on the filtering in the script output file that looks like this: 

```
Length of filtered alignment: 528
Number of columns removed: 12897
Length of the original alignment: 13425
Number of sequences used to construct filter: 75281
```

We also get a new summary that looks like this:  

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	526     250     0	3	1
2.5%-tile:	1	528     262     0	3	6369
25%-tile:	1	528     262     0	4	63688
Median:         1	528     262     0	4	127375
75%-tile:	1	528     262     0	4	191062
97.5%-tile:     1	528     262     0	6	248381
Maximum:        2	528     280     0	8	254749
Mean:   1	527     262     0	4
# of unique seqs:	75281
total # of seqs:        254749

```

From this summary we see that the alignment window spans ~500 bp and the length of our sequences is about 253 nt. We have a maximum polymer of 8 as specified in our settings above.  

### <a name="Precluster"></a> **7. Polish the data with pre clustering**     

Now we need to further polish and cluster the data with pre.cluster. The purpose of this step is to remove noise due to sequencing error. The rational behind this step assumes that the most abundant sequences are the most trustworthy and likely do not have sequencing errors. Pre-clustering then looks at the relationship between abundant and rare sequences - rare sequences that are "close" (e.g., 1 nt difference) to highly abundant sequences are likely due to sequencing error. This step will pool sequences and look at the maximum differences between sequences within this group to form ASV groupings. 

In this step, the number of sequences is not reduced, but they are grouped into amplicon sequence variants ASV's which reduces the error rate. V4 region has the lowest likelihood of errors, so the error rate is going to be lower than for other variable regions.  

Other programs that conduct this "denoising" are DADA2, UNOISE, and DEBLUR. However, these programs remove the rare sequences, which can distort the relative abundance of remaining sequences. DADA2 also removes all sigletons (sequences with single representation) which disproportionately affects the sequence relative abundance. Mothur avoids the removal of rare sequences for this reason. 

We will first add code to identify unique sequences after the filtering steps above.  

```
unique.seqs(fasta=mcap.trim.contigs.good.unique.good.filter.fasta, count=mcap.trim.contigs.good.good.count_table)
```

We will then perform the pre-clustering a default of 1 nt difference. Diffs can be changed according to your requirements.  

```
pre.cluster(fasta=mcap.trim.contigs.good.unique.good.filter.unique.fasta, count=mcap.trim.contigs.good.unique.good.filter.count_table, diffs=1) 
```

Finally, we will run another summary.  

```
summary.seqs(fasta=mcap.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.count_table)
```

Write and run the script.  

```
nano precluster.sh
```

```
#!/bin/bash
#SBATCH --job-name="precluster"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="precluster_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="precluster_output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -q putnamlab

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#unique.seqs(fasta=mcap.trim.contigs.good.unique.good.filter.fasta, count=mcap.trim.contigs.good.good.count_table)"

mothur "#pre.cluster(fasta=mcap.trim.contigs.good.unique.good.filter.unique.fasta, count=mcap.trim.contigs.good.unique.good.filter.count_table, diffs=1)"

mothur "#summary.seqs(fasta=mcap.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.count_table)"
```

```
sbatch precluster.sh
```

Mothur output the following files:   

```
Output File Names: 
mcap.trim.contigs.good.unique.good.filter.count_table
mcap.trim.contigs.good.unique.good.filter.unique.fasta
```

Then, the pre-clustering step outputs the following files: 

```
Output File Names:
mcap.trim.contigs.good.unique.good.filter.unique.precluster.fasta
mcap.trim.contigs.good.unique.good.filter.unique.precluster.count_table
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH174.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH175.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH176.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH177.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH178.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH179.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH180.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH181.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH182.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH182.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH183.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH184.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH185.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH186.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH187.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH188.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH189.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH190.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH191.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH192.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH193.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH194.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH195.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH196.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH201.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH202.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH203.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH204.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH205.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH206.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH207.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH208.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH209.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH210.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH211.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH212.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH213.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH214.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH215.map
mcap.trim.contigs.good.unique.good.filter.unique.precluster.WSH216.map
```  

There are a ton of files here, but the two most important files are: 

```
mcap.trim.contigs.good.unique.good.filter.unique.precluster.fasta
mcap.trim.contigs.good.unique.good.filter.unique.precluster.count_table
``` 

The other files have text for maps of sequence name, errors, abundance, differences, and the filtered sequence for each sample.   

Finally, we get the output from the summary:   

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	527     250     0	3	1
2.5%-tile:	1	528     262     0	3	6369
25%-tile:	1	528     262     0	4	63688
Median:         1	528     262     0	4	127375
75%-tile:       1 	528     262     0	4	191062
97.5%-tile:     1	528     262     0	6	248381
Maximum:        2       528     280     0	8	254749
Mean:   1	527     262     0	4
# of unique seqs:       10774
total # of seqs:        254749
```

Note that the number of unique sequences has decreased from 75,281 to 10,774 as expected since we are clustering sequences that are within 1 nt difference from each other. 

### <a name="Chimera"></a> **8. Identify chimeras**  

Now we will remove chimeras using the dereplicate method. In this method, we are again using the assumption that the highest abundance sequences are most trustworthy. Chimeras are sequences that did not extend during PCR and then served as templates for other PCR products, forming sequences that are partially from one PCR product and partially from another. This program looks for chimeras by comparing each sequences to the next highest abundance sequences to determine if a sequence is a chimera of the more abundance sequences.  

We will use the `chimera.vsearch` function to identify chimeras: 

```
chimera.vsearch(fasta=mcap.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=T)
``` 

We will then remove the identified chimeras with `remove.seqs`:  

```
remove.seqs(fasta=mcap.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)
```

Finally, we will run a new summary:  

```
summary.seqs(fasta=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)
```

This step requires an executable program called "vsearch". This is now available as a module on Andromeda. If you are working in Andromeda, we will load the module.  If the module is not available or you are working on a different system, you can install vsearch in your working directory with the following commands: 

```
wget https://github.com/torognes/vsearch/archive/v2.21.0.tar.gz
tar xzf v2.21.0.tar.gz
cd vsearch-2.21.0
./autogen.sh
./configure CFLAGS="-O3" CXXFLAGS="-O3"
make
make install  # as root or sudo make install
```

Write and run a script to do these steps.  

```
nano chimera.sh
```

```
#!/bin/bash
#SBATCH --job-name="chimera"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="chimera_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="chimera_output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -q putnamlab

module load Mothur/1.46.1-foss-2020b

module load VSEARCH/2.18.0-GCC-10.2.0

mothur

mothur "#chimera.vsearch(fasta=mcap.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=T)"

mothur "#remove.seqs(fasta=mcap.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)"

mothur "#summary.seqs(fasta=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)"

mothur "#count.groups(count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)"

```

```
sbatch chimera.sh
```


This script outputs the following files: 

```
Output File Names:
mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table
mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.chimeras
mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos
```

The new summary looks like this:

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	527     250     0	3	1
2.5%-tile:	1	528     262     0	3	6290
25%-tile:	1	528     262     0	4	62900
Median:         1	528     262     0	4	125799
75%-tile:	1	528     262     0	4	188698
97.5%-tile:     1	528     262     0	6	245308
Maximum:        2	528     280     0	8	251597
Mean:   1	527     262     0	4
# of unique seqs:	8818
total # of seqs:        251597
```

The program identified and removed ~3% chimeras.    

We can look at a count of the number of sequences per sample: 

```
WSH174 contains 11208.
WSH175 contains 1700.
WSH176 contains 7455.
WSH177 contains 28911.
WSH178 contains 6786.
WSH179 contains 46658.
WSH180 contains 922.
WSH181 contains 682.
WSH182 contains 279.
WSH183 contains 128.
WSH184 contains 318.
WSH185 contains 23303.
WSH186 contains 1534.
WSH187 contains 9333.
WSH188 contains 4071.
WSH189 contains 194.
WSH190 contains 630.
WSH191 contains 278.
WSH192 contains 192.
WSH193 contains 3254.
WSH194 contains 237.
WSH195 contains 3540.
WSH196 contains 286.
WSH201 contains 434.
WSH202 contains 457.
WSH203 contains 2426.
WSH204 contains 1931.
WSH205 contains 3171.
WSH206 contains 5005.
WSH207 contains 6691.
WSH208 contains 13514.
WSH209 contains 10770.
WSH210 contains 9071.
WSH211 contains 14187.
WSH212 contains 5830.
WSH213 contains 4578.
WSH214 contains 14119.
WSH215 contains 1786.
WSH216 contains 5728.

Size of smallest group: 128.
```

The smallest group has 128 sequences and we have 13 samples with <1,000 sequences. We will further look at this sampling depth.  

### <a name="Classify"></a> **9. Classifying sequences**  

Now our sequences are clean and ready for classification!  

We will use the training set downloaded above from the silva database through the [Mothur wiki](https://mothur.org/wiki/classify.seqs/). 

#### Classify sequences  

We will use the `classify.seqs` command: 

```
classify.seqs(fasta=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds.tax)
```

The training and taxonomy files were from downloads. The mothur [classify.seq wiki page](https://mothur.org/wiki/classify.seqs/) has the reference and training sets that you can download and put into the directory to use for analyzing data. These sets include adding chlorophyll, mitchondria, etc to identify and remove these. 

*Outside andromeda, download these files and then transfer to Andromeda.*     

The output file from `classify.seqs()` ending in .taxonomy has the name of sequence and the classification with % confidence in parentheses for each level. It will end at the level that is has confidence.  

The tax.summary file has the taxonimc level, the name of the taxonomic group, and the number of sequences in that group for each sample.  

We will also remove sequences that are classified to Chloroplast, Mitochondria, Unknown (not bacteria, archaea, or eukaryotes), Archaea, and Eukaryotes. 

```
remove.lineage(fasta=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon=Cyanobacteria_Chloroplast-Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
```

Write a script to classify and remove lineages. 

```
nano classify.sh
```

```
#!/bin/bash
#SBATCH --job-name="classify"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="classify_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="classify_output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -q putnamlab

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#classify.seqs(fasta=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds.tax)"

mothur "#remove.lineage(fasta=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon=Cyanobacteria_Chloroplast-Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)"

mothur "#summary.seqs(fasta=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)"
``` 

```
sbatch classify.sh
```

Classify.seqs will output the following files: 

```
Output File Names: 

mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy
mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary

```

The output file .taxonomy has name of sequence and the classification with % confidence in parentheses for each level. It will end at the level that is has confidence.  

The tax.summary file has the taxonimc level, the name of the taxonomic group, and the number of sequences in that group for each sample.  

```
head mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy
head mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary
```

The remove.lineage command will remove sequences and provide a report in the output file. The fasta file is removing unique sequences, the count file is removing individual sequences.   

I altered the script to only pull out each lineage individually and ran for each lineage to get the individual stats. I then ran the final script with removing all lineages together to proceed with analysis.  

*Removing only Chloroplast*     
Removed 1055 sequences from your fasta file.  
Removed 45530 sequences from your count file.  

*Removing only Mitochondria*     
Removed 2 sequences from your fasta file.  
Removed 4 sequences from your count file.   

*Removing only unknown domain*     
Removed 76 sequences from your fasta file.     
Removed 719 sequences from your count file.     

*Removing only Archaea*   
Removed 22 sequences from your fasta file.  
Removed 59 sequences from your count file.   

*Removing only Eukaryotes* 
No contaminants to remove. This is expected since we aligned to the bacterial database.  

*Removing Cyanobacteria_Chloroplast*  
Removed 1186 sequences from your fasta file.    
Removed 50976 sequences from your count file.    
Note that this is largely what is removed in the Chloroplast step, so the Chloroplast and Cyanobacteria chloroplast overlap in sequences.  

*Removing all lineages*   
Removed 1286 sequences from your fasta file.  
Removed 51758 sequences from your count file.    

We also can see the following summary:  

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	527     250     0	3	1
2.5%-tile:	1	528     262     0	3	4996
25%-tile:	1	528     262     0	4	49960
Median:         1	528     262     0	4	99920
75%-tile:	1	528     262     0	5	149880
97.5%-tile:     1	528     262     0	6	194844
Maximum:        2	528     280     0	8	199839
Mean:   1	527     262     0	4
# of unique seqs:	7532
total # of seqs:        199839
```

We now have 199,839 total sequences and 7,532 unique sequences.  

The script will then output the following files: 

```
Output File Names:

mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta
mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table
mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy
mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.accnos
mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table
mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta

```

Now we have the taxon removed that we do not want in our dataset.  

This is the end of the pipeline for curating our sequencing - we have corrected for pcr errors, removed chimeras, and removed sequences outside of taxon of interest. Now we can move onto OTU clustering!

View a sample of the taxonomy summary dataset. 
 
```
head mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy
```

```
M00763_26_000000000-K4TML_1_1115_18196_2703	Bacteria(100);"Proteobacteria"(82);"Proteobacteria"_unclassified(82);"Proteobacteria"_unclassified(82);"Proteobacteria"_unclassified(82);"Proteobacteria"_unclassified(82);
M00763_26_000000000-K4TML_1_2115_7713_19993	Bacteria(100);"Bacteroidetes"(100);Flavobacteria(100);"Flavobacteriales"(100);Flavobacteriaceae(100);Flavobacteriaceae_unclassified(100);
M00763_26_000000000-K4TML_1_1118_9895_8086	Bacteria(100);"Bacteroidetes"(100);Flavobacteria(100);"Flavobacteriales"(100);Flavobacteriaceae(100);Winogradskyella(87);
M00763_26_000000000-K4TML_1_1118_27290_8095	Bacteria(100);"Planctomycetes"(100);"Planctomycetacia"(100);Planctomycetales(100);Planctomycetaceae(100);Planctomycetaceae_unclassified(100);
M00763_26_000000000-K4TML_1_1106_19148_22444	Bacteria(100);"Planctomycetes"(100);"Planctomycetacia"(100);Planctomycetales(100);Planctomycetaceae(100);Planctomycetaceae_unclassified(100);
M00763_26_000000000-K4TML_1_2115_26839_20267	Bacteria(100);"Proteobacteria"(100);Gammaproteobacteria(100);Oceanospirillales(97);Oceanospirillaceae(97);Thalassolituus(86);
M00763_26_000000000-K4TML_1_1108_13024_9786	Bacteria(100);"Proteobacteria"(100);Gammaproteobacteria(100);Alteromonadales(88);Alteromonadales_unclassified(88);Alteromonadales_unclassified(88);
M00763_26_000000000-K4TML_1_1118_20226_8150	Bacteria(100);"Proteobacteria"(100);Alphaproteobacteria(100);Rhodobacterales(100);Rhodobacteraceae(100);Rhodobacteraceae_unclassified(100);
M00763_26_000000000-K4TML_1_2105_8942_6450	Bacteria(100);Bacteria_unclassified(100);Bacteria_unclassified(100);Bacteria_unclassified(100);Bacteria_unclassified(100);Bacteria_unclassified(100);
M00763_26_000000000-K4TML_1_1113_14691_17286	Bacteria(100);"Proteobacteria"(100);Gammaproteobacteria(100);Pseudomonadales(100);Pseudomonadaceae(100);Pseudomonas(98);
```

### <a name="Cluster"></a> **10. Cluster for OTUs**  

*In this analysis, we will cluster to OTU level (cutoff=0.03). For ASV clustering, you can move directly to the make.shared step, skipping the dist.seqs and cluster steps because mothur pre-clustering occurs as the ASV level.*  

First, we will calculate the pairwise distances between sequences.  

We will first use the `dist.seqs` command.  

```
dist.seqs(fasta=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta)
```

We will then `cluster` using the following command: 

```
cluster(column=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, cutoff=0.03)
```

This will run a line for each iteration of clustering. This is run until the Matthews correlation coefficient (MCC) value is maximized. A high MCC = high confidence in clustering. MCC is optimized by randomly aligning sequences to OTU's and calculating the correlation coefficient. Then sequences are moved between OTU's to see if the MCC is improved. This is repeated many times until the MCC is maximized. This method is fast and RAM efficient. AKA Opticlust.  

Alternatively, you could run `cluster.split` to calculate the distance and cluster within each taxonomic level (order in this case) in parallel, then bring the data back together. This will generate a file that indicates which sequence are in which OTU's. We do not need this in this analysis.   

```
cluster.split(fasta=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, taxlevel=4, cutoff=0.03, splitmethod=classify)
```

We would move directly to the make.shared step (below) if you want to use ASV since the distance and clustering steps do not need to be completed for ASVs.  

Next we will make a shared file. This .shared file has the label for the OTU, sample name, the number of OTU's and then the number of time each OTU appears in each sample. 

This file will be the basis of what we will do to measure richness of communities compared to each other.  

We want to keep the shared file and a consensus taxonomy file. 

```
make.shared(list=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
```

The classify.otu command will then output a concensus cons.taxonomy file that has the taxonomic information for each OTU. Use label=ASV to specify ASV in taxonomy names if starting from the make.shared step for ASV's.    

```
classify.otu(list=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy)
```

Then we can rename the files to something more useful.  

```
rename.file(taxonomy=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy, shared=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared)
```

Finally, view a count of the number of sequences in each sample.  

```
count.groups(shared=mcap.opti_mcc.shared)
``` 

Generate a script to run these commands.    

```
nano cluster.sh
```

```
#!/bin/bash
#SBATCH --job-name="cluster"
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="cluster_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="cluster_output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -q putnamlab

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#dist.seqs(fasta=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta)"

mothur "#cluster(column=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, cutoff=0.03)"

mothur "#make.shared(list=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)"

mothur "#classify.otu(list=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy)"

mothur "#rename.file(taxonomy=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy, shared=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared)"

mothur "#count.groups(shared=mcap.opti_mcc.shared)"
```

```
sbatch cluster.sh
```

Dist.seqs outputs the following files: 

```
Output File Names:
mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist

```

Cluster outputs the following files:  

```
Output File Names:
mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list
mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.steps
mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.sensspec

```

Make.shared outputs the following file: 

```
Output File Names:
mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared
```

Finally, classify.otu outputs the following file:  

```
Output File Names:
mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy
mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.tax.summary
```

The rename function at the end will rename our files to something more useful. We now have a "taxonomy" file and a "shared" file.  

The files we now care about are: 

```
Current files saved by mothur:
list=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list
shared=mcap.opti_mcc.shared
taxonomy=mcap.taxonomy

constaxonomy=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy
count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table

```

The count of sequences in each file are:  

```
WSH174 contains 9448.
WSH175 contains 1319.
WSH176 contains 5407.
WSH177 contains 26867.
WSH178 contains 4433.
WSH179 contains 24988.
WSH180 contains 835.
WSH181 contains 513.
WSH182 contains 266.
WSH183 contains 121.
WSH184 contains 309.
WSH185 contains 12908.
WSH186 contains 1111.
WSH187 contains 4007.
WSH188 contains 1678.
WSH189 contains 172.
WSH190 contains 621.
WSH191 contains 278.
WSH192 contains 186.
WSH193 contains 3236.
WSH194 contains 234.
WSH195 contains 3536.
WSH196 contains 286.
WSH201 contains 426.
WSH202 contains 449.
WSH203 contains 2381.
WSH204 contains 1927.
WSH205 contains 3156.
WSH206 contains 4987.
WSH207 contains 6619.
WSH208 contains 13451.
WSH209 contains 10739.
WSH210 contains 8992.
WSH211 contains 14019.
WSH212 contains 5805.
WSH213 contains 4498.
WSH214 contains 13313.
WSH215 contains 1656.
WSH216 contains 4662.

Size of smallest group: 121.
```

*Now we have classified taxonomic data and we are ready to move onto subsampling and calculating statistics*  


### <a name="Subsample"></a> **11. Subsampling for Sequencing Depth**   

**The remaining commands can all be run in interactive mode. In Andromeda, use the following commands:**   

```
interactive 
module load Mothur/1.46.1-foss-2020b
mothur
```

This will open the interactive mothur terminal.  In order to use any bash commands like `nano`, you have to first `quit` in mothur. Once you want to return to mothur, use the `mothur` command.  

#### Subsampling for sequencing depth 

As we can see in rarefaction curves, we know that OTU discovery is a function of sequencing depth. Therefore, we need to account for sequencing depth in our analysis.  

With the count table above, we see that the sample with the smallest number of groups is 95. To account for uneven sampling depth, we will need to rarify our samples. We can do this with the `sub.sample` command that will automatically use the smallest group size or we can manually set a cut off. If we set a threshold higher than the number of samples, those samples are removed from the group.  

```
sub.sample(shared=mcap.opti_mcc.shared, size=1000) 
```

```
Output files: 
mcap.opti_mcc.0.03.subsample.shared
```

It may be helpful to run this with multiple iterations. If you do not include a size=# argument, then mothur will automatically set to the lowest sequence. If I set this to 1000, for example the output shows that several samples are removed because there are <1000 samples: 

```
WSH180 contains 835. Eliminating.
WSH181 contains 513. Eliminating.
WSH182 contains 266. Eliminating.
WSH183 contains 121. Eliminating.
WSH184 contains 309. Eliminating.
WSH189 contains 172. Eliminating.
WSH190 contains 621. Eliminating.
WSH191 contains 278. Eliminating.
WSH192 contains 186. Eliminating.
WSH194 contains 234. Eliminating.
WSH196 contains 286. Eliminating.
WSH201 contains 426. Eliminating.
WSH202 contains 449. Eliminating.
Sampling 1000 from each group.
```  

Sample names are found in the [spreadsheet here](https://docs.google.com/spreadsheets/d/1lLvCp-RoRiBSGZ4NBPwi6cmZuozmfS20OJ7hBIueldU/edit#gid=0). QC information from extraction can be found in the [spreadsheet here](https://docs.google.com/spreadsheets/d/1Ew0AOs88n1i6wEyuqbv3tsRXwMjhqyYEeZX5sACFdDY/edit#gid=0).   

WSH180 = Attached Recruit (231 hpf)  
WSH181 = TP1  
WSH182 = TP1  
WSH183 = TP1  
WSH184 = TP1  
WSH189 = TP2  
WSH190 = TP2  
WSH191 = TP2  
WSH192 = TP2  
WSH194 = TP3  
WSH196 = TP3  
WSH201 = TP4  
WSH202 = TP4  

This sub-sampling to 1000 sequences removes all of TP1 (eggs) and TP2 (embryos 5 hpf). We remove 2/4 of TP3 and TP4 (embryos 38 hpf and embryos 65 hpf). We also remove  one attached recruit sample.  

We know that DNA concentrations were very low for the TP1-TP3 samples so we may need to remove them from the analysis or would have to consider that the microbial communities will be different due to sequencing depth and DNA extraction artifacts.    

If we need to standardized sampling depth for all of our samples, then we would need to remove these samples. 

#### Subsampling the sample set  

We are going to generate a rarefaction curve and subsample to 1000 sequences. Calculating the observed number of OTUs using the Sobs metric (observed number of taxa) with a freq=100 to output the data for every 100 sequences samples - if you did every single sample that would be way too much data to output.  

```
rarefaction.single(shared=mcap.opti_mcc.shared, calc=sobs, freq=100)
```

```
Output File Names: 
mcap.opti_mcc.groups.rarefaction
```

We will want to keep this rarefaction file for our final output.  

You can look at this file with `nano mcap.opti_mcc.groups.rarefaction`. It contains the number of OTUs detected for every 100 sequences for each sample. Remember OTU = 0.03 cut off here, so that is the OTU distinction. Series will be trucated after the number of sequences available in the sample.  

Next, rarefy the data to the minimum number of sequences, then the next command will rarefy and pull out that many sequences from each sample. If you want to sample a specific number, use subsample=#. You would want to do this if you have a super small sequence number in one group and you dont want to truncate the data for all groups by that much. If you leave subsample=T, mothur will calculate to the minimum size. Subsamples 1000 times and calculates the metrics we want with the number of otus. Sampling is random.   

All of the available metric calculations (calc) are [found here](https://mothur.org/wiki/calculators/).  

I am going to use subsample=1000 to remove the low sequence count stages.  

```
summary.single(shared=mcap.opti_mcc.shared, calc=nseqs-sobs-chao-shannon-invsimpson, subsample=1000)
```

In a new Andromeda window (outside mothur).  
```
less mcap.opti_mcc.groups.ave-std.summary
```

Open this file to view statistics for each sample. The method column has ave (average of each of the metrics across the 1000 iteractions) or std (this is the std deviation of the average calculation).  

Here, I will produce datasets both with and without subsampling. We will use the subsampled data for our analysis and we can explore/describe the data from the samples below the subsampling threshold.    

*I recommend running iterations of these analyses with different thresholds for the commands above to see how our data filtering is changed. Then we can make data-driven decisions.*    

### <a name="Statistics"></a> **12. Calculate Ecological Statistics**  

Alpha diversity = diversity within a sample 
Beta diversity = diversity between samples

We will calculate beta diversity metrics which will then be used for ecological analyses.  

#### Calculate Beta Diversity  

We will use the `dist.shared` function to calculate distances between samples. Here we are going to use Bray-Curtis distances. Many different metrics can be used (e.g., Theta YC). All of the available metric calculations (calc) are [found here](https://mothur.org/wiki/calculators/).  

We can use this function with or without `subsample=1000` or `subsample=T`. We will use the bray-curtis distances for this calculation. This step is important for generating files that we will use later in R. 

Run without subsampling: 

```
dist.shared(shared=mcap.opti_mcc.shared, calc=braycurtis)
```

```
Output File Names: 
mcap.opti_mcc.braycurtis.0.03.lt.dist
```

If we run a command with subsample=1000, then the following files are output: 

```
dist.shared(shared=mcap.opti_mcc.shared, calc=braycurtis, subsample=1000)
```

```
Output File Names: 
mcap.opti_mcc.braycurtis.0.03.lt.ave.dist
mcap.opti_mcc.braycurtis.0.03.lt.std.dist
``` 
We have average and std files because the function is run over 1000 iterations generating ave and std.    

We can use these different files to compare a non-subset to a subset dataset. These comparison files will be:    

- No subsampling `mcap.opti_mcc.braycurtis.0.03.lt.dist`   
- Subsampled `mcap.opti_mcc.braycurtis.0.03.lt.ave.dist`     

*After this point, analyses can be run in R. See the bottom of this document for transferring files to your computer. The following sections are examples of how to run the statistical analyses in mothur if desired.*    

#### Run a PCoA and NMDS   

Run a PCoA on the Bray Curtis distances without subsampling.    

```
pcoa(phylip=mcap.opti_mcc.braycurtis.0.03.lt.dist)
```

```
Processing...
Rsq 1 axis: 0.419531
Rsq 2 axis: 0.661573
Rsq 3 axis: 0.645606

Output File Names: 
mcap.opti_mcc.braycurtis.0.03.lt.pcoa.axes
mcap.opti_mcc.braycurtis.0.03.lt.pcoa.loadings
```

Run a PCoA on the Bray Curtis distances with subsampling.    

```
pcoa(phylip=mcap.opti_mcc.braycurtis.0.03.lt.ave.dist)
```

```
Processing...
Rsq 1 axis: 0.48495
Rsq 2 axis: 0.512541
Rsq 3 axis: 0.71141

Output File Names: 
mcap.opti_mcc.braycurtis.0.03.lt.ave.pcoa.axes
mcap.opti_mcc.braycurtis.0.03.lt.ave.pcoa.loadings
```

These outputs show the R squared value on each axis.  

We see that the R squared values are increased with the data set with subsampling.  

Now run an NMDS on the dataset without subsampling.  

```
nmds(phylip=mcap.opti_mcc.braycurtis.0.03.lt.dist)
``` 

```
Number of dimensions:	2
Lowest stress :	0.308031
R-squared for configuration:	0.71123

Output File Names: 
mcap.opti_mcc.braycurtis.0.03.lt.nmds.iters
mcap.opti_mcc.braycurtis.0.03.lt.nmds.stress
mcap.opti_mcc.braycurtis.0.03.lt.nmds.axes
``` 

Now run an NMDS on the dataset with subsampling.  

```
nmds(phylip=mcap.opti_mcc.braycurtis.0.03.lt.ave.dist)
``` 

```
Number of dimensions:	2
Lowest stress :	0.281742
R-squared for configuration:	0.703092

Output File Names: 
mcap.opti_mcc.braycurtis.0.03.lt.ave.nmds.iters
mcap.opti_mcc.braycurtis.0.03.lt.ave.nmds.stress
mcap.opti_mcc.braycurtis.0.03.lt.ave.nmds.axes
``` 

If you run NMDS each time, you will get the same clustering/points but they might be rotated on different axes.  

The iters file will also have information on stress values for each time it ran the tests.

#### Run AMOVA tests  

To run statistical tests, we first need to generate a metadata sheet.  

This will be a file named mcap.design and created in Excel as a tab delimited file. One column will be group (e.g., WSH192) and one column will be treatment (e.g., egg). 

Outside of Andromeda, copy the file.  

```
scp  ~/MyProjects/EarlyLifeHistory_Energetics/Mcap2020/Data/16S/mcap_design.txt ashuffmyer@ssh3.hac.uri.edu:/data/putnamlab/ashuffmyer/AH_MCAP_16S/mothur
```

An AMOVA tests whether the centroids are different between groupings.  

Test for differences in lifestage in data without subsampling:  

```
amova(phylip=mcap.opti_mcc.braycurtis.0.03.lt.dist, design=mcap_design.txt, iters=1000)
```

This outputs pairwise comparisons between stages. There is a significant effect of life stage.  

```
EggFertilized-Embryo1-Larvae1-Larvae2-Larvae3-Larvae4-Larvae6-Recruit1-Recruit1plug-Recruit2-Recruit2plug	Among	Within	Total
SS	8.24943	6.40084	14.6503
df	10	28	38
MS	0.824943	0.228602

Fs:	3.60865
p-value: <0.001*
```

```
Output File Names: 
mcap.opti_mcc.braycurtis.0.03.lt.amova
```

Test for differences in lifestage in data with subsampling: 

```
amova(phylip=mcap.opti_mcc.braycurtis.0.03.lt.ave.dist, design=mcap_design.txt, iters=1000)

```

```
Larvae1-Larvae2-Larvae3-Larvae4-Larvae6-Recruit1-Recruit1plug-Recruit2-Recruit2plug	Among	Within	Total
SS	5.26248	2.72787	7.99034
df	8	17	25
MS	0.65781	0.160463

Fs:	4.09946
p-value: <0.001*

```

There are significant differences by lifestage in the data that is subsampled as well. 

```
Output File Names: 
mcap.opti_mcc.braycurtis.0.03.lt.ave.amova
```

#### Run HOMOVA tests  

Run a HOMOVA test to test whether the variation between groups is different. In the last step we tested for differences in group centroids. Now we will test for differences in spread. This would be a good test for questions around variance and stability in microbiomes.  

Test without subsampling.  

```
homova(phylip=mcap.opti_mcc.braycurtis.0.03.lt.dist, design=mcap_design.txt)
```

```
HOMOVA	BValue	P-value	SSwithin/(Ni-1)_values
EggFertilized-Embryo1-Larvae1-Larvae2-Larvae3-Larvae4-Larvae6-Recruit1-Recruit1plug-Recruit2-Recruit2plug	-nan	<0.001*	0.212107	0.171277	0.298598	0.218722	0.127314	0.0710606	0.286777	-nan	0.324554	0.254493	0.380313


mcap.opti_mcc.braycurtis.0.03.lt.homova
```

There is a significant difference in variation in microbiomes between lifestages without subsampling. 


Test with subsampling.  

```
homova(phylip=mcap.opti_mcc.braycurtis.0.03.lt.ave.dist, design=mcap_design.txt)
```

```
HOMOVA	BValue	P-value	SSwithin/(Ni-1)_values
Larvae1-Larvae2-Larvae3-Larvae4-Larvae6-Recruit1-Recruit1plug-Recruit2-Recruit2plug	-nan	<0.001*	0.062777	0.0539804	0.0634763	0.0732464	0.25485	-nan	0.260393	0.219021	0.21717

mcap.opti_mcc.braycurtis.0.03.lt.ave.homova
```

Note there are NA's due to removal of egg and embryo groups.  

There is a significant difference in variation in microbiomes between lifestages with subsampling as well.  

### <a name="Output"></a> **13. Output data for analysis in R**  

We can now move some files into R for further analysis. These are the primary files we will be using - the bray curtis distance matrix (distances between samples), the taxonomy file (taxonomy of each otu), and the shared file (otu abundance in each sample). From these we can run PCoA, NMDS, and statistical testing in R (as done above in mothur, but I prefer working in R).   

Note that here I am taking files for both subsampling and file without subsampling since I am comparing these methods.  

Bray-Curtis matrices:     

```
mcap.opti_mcc.braycurtis.0.03.lt.dist
mcap.opti_mcc.braycurtis.0.03.lt.ave.dist
```

Taxonomy file:  

```
mcap.taxonomy
```

Shared files:  

```
mcap.opti_mcc.0.03.subsample.shared
mcap.opti_mcc.shared
```

Rarefaction file:  
```
mcap.opti_mcc.groups.rarefaction
```

Outside Andromeda do the following for each file:  

```
scp ashuffmyer@ssh3.hac.uri.edu:/data/putnamlab/ashuffmyer/AH_MCAP_16S/mothur/mcap.opti_mcc.braycurtis.0.03.lt.dist ~/MyProjects/EarlyLifeHistory_Energetics/Mcap2020/Output/16S/mothur

scp ashuffmyer@ssh3.hac.uri.edu:/data/putnamlab/ashuffmyer/AH_MCAP_16S/mothur/mcap.opti_mcc.braycurtis.0.03.lt.ave.dist ~/MyProjects/EarlyLifeHistory_Energetics/Mcap2020/Output/16S/mothur

scp ashuffmyer@ssh3.hac.uri.edu:/data/putnamlab/ashuffmyer/AH_MCAP_16S/mothur/mcap.taxonomy ~/MyProjects/EarlyLifeHistory_Energetics/Mcap2020/Output/16S/mothur

scp ashuffmyer@ssh3.hac.uri.edu:/data/putnamlab/ashuffmyer/AH_MCAP_16S/mothur/mcap.opti_mcc.0.03.subsample.shared ~/MyProjects/EarlyLifeHistory_Energetics/Mcap2020/Output/16S/mothur

scp ashuffmyer@ssh3.hac.uri.edu:/data/putnamlab/ashuffmyer/AH_MCAP_16S/mothur/mcap.opti_mcc.shared ~/MyProjects/EarlyLifeHistory_Energetics/Mcap2020/Output/16S/mothur

scp ashuffmyer@ssh3.hac.uri.edu:/data/putnamlab/ashuffmyer/AH_MCAP_16S/mothur/mcap.opti_mcc.groups.rarefaction ~/MyProjects/EarlyLifeHistory_Energetics/Mcap2020/Output/16S/mothur
```

After copying to the computer, I resave these files at tab delimited (.txt) before loading into R.  

Finally, move sequences to personal computer. 

```
scp ashuffmyer@ssh3.hac.uri.edu:/data/putnamlab/ashuffmyer/AH_MCAP_16S/raw_data/*.gz ~/MyProjects/EarlyLifeHistory_Energetics/Mcap2020/Data/16S/Sequences

scp ashuffmyer@ssh3.hac.uri.edu:/data/putnamlab/ashuffmyer/AH_MCAP_16S/raw_data/*.txt ~/MyProjects/EarlyLifeHistory_Energetics/Mcap2020/Data/16S/Sequences

```