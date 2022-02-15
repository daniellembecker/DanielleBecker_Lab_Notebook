---
layout: post
title: 16S Pipeline in Mothr
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
13. [Popluation analyses](#Population)  
14. [Output data for R analysis](#Output)   

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

The code will look like this:  

```
make.contigs(inputdir=., outputdir=., file=mcap.files, oligos=oligos.oligos)
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
#SBATCH --account=putnamlab
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="contigs_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="contigs_output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -q putnamlab

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#make.file(inputdir=., type=gz, prefix=mcap)"

mothur "#make.contigs(inputdir=., outputdir=., file=mcap.files, oligos=oligos.oligos)"

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

If you have barcodes, you can use the `oligos` command to direct to a file that has the primer and barcode information during the contig step. [See more information here](https://mothur.org/wiki/oligos_file/#).  

The following files were output from the make.contigs() step:  

```
mcap.trim.contigs.summary
```

When you run summary.seqs() you will get an output table like the one below and it will output a file with this summary information.  

```
nano contigs_output_script

```
Scroll to the bottom to find the summary information.  

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	23	23	0	2	1
2.5%-tile:	1	292     292     0	4	7948
25%-tile:	1	292     292     0	4	79472
Median:         1	292     292     0	4	158943
75%-tile:	1	301     301     0	5	238414
97.5%-tile:     1	357     357     5	6	309938
Maximum:        1	563     563     56	250     317885
Mean:   1	298     298     0	4
# of Seqs:	317885

```

This table shows quantile values about the distribution of sequences for a few things:  

- *Start position*: All at 1 now, will start at different point after some QC.   
- *End position*: We see that there are some sequences that are very short and we may need to remove those later.   
- *Number of bases*: length (we see most are in expected range here, but one is super long! This might tell us there is no overlap so they are butted up against each other. We will remove things like this.  
- *Ambigs*: Number of ambiguous calls in sequences. Here there are a few that have ambiguous base calls. We will remove any sequence with an ambiguous call or any longer than we would expect for V4 region.   
- *Polymer*: Length of polymer repeats.    
- *NumSeqs*: Number of sequences.  

#### Check that primers are gone  

We should now check that the primers were removed in the `make.contigs` step. 

```
head mcap.trim.contigs.fasta 
```

The primers we are looking for are:  

`F GTGCCAGCMGCCGCGGTAA R GGACTACNVGGGTWTCTAAT`

Here are the first 10 sequences in the file: 

```
>M00763_26_000000000-K4TML_1_1101_13007_1601	ee=3.43915	fbdiffs=0(match), rbdiffs=0(match) fpdiffs=0(match), rpdiffs=0(match) 
AAGATACGGGGTCCAGCCGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGTTTGTTAAGTGAGATGTGAAAGCCCAGGGCTCAACCTTGGAACTGCATCTCATACTGGCAGGCTAGAGTATGGTAGAGGGAGGTAGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAGGAATACCAGTGGCGAAGGCGGCCTCCTGGACTAATACTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTAGTAGTCC
>M00763_26_000000000-K4TML_1_1101_15107_1629	ee=1.2146	fbdiffs=0(match), rbdiffs=0(match) fpdiffs=0(match), rpdiffs=0(match) 
ATGAGACAGGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTATCCGGAATCACTGGGTTTAAAGGGTGCGTAGGCGGCGCTATAAGTCAGAGGTGAAAGGCCACCGCTTAACGGTGGGACTGCCTTTGATACTGTAGTGCTTGAATCAGGTTGAGGTAGGCGGAATGTGACATGTAGCGGTGAAATGCTTAGATATGTCATAGAACACCAATTGCGAAGGCAGCTTGCTAGACCTTGATTGACGCTGAGGCACGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTAGTAGTCC
>M00763_26_000000000-K4TML_1_1101_16047_1642	ee=2.46866	fbdiffs=0(match), rbdiffs=0(match) fpdiffs=0(match), rpdiffs=0(match) 
AAGGGACAGGTGCCAGCATCCGCGGTAATACGGAGGATGCAAGCGTTATCCGGAATCATTGGGTTTAAAGGGTCCGTAGGTGGACAATTAAGTCAGAGGTGAAATCCTGCAGCTCAACTGTAGAATTGCCTTTGATACTGGTTGTCTTGAATTATTGTGAAGTGGTTAGAATATGTAGTGTAGCGGTGAAATGCATAGATATTACATAGAATACCAATTGCGAAGGCAGATCACTAACAATATATTGACACTGATGGACGAAAGCGTGGGGAGCGAACGGGATTAGATACCCGGGTAGTCC
>M00763_26_000000000-K4TML_1_1101_23424_1609	ee=15.1987	fbdiffs=0(match), rbdiffs=0(match) fpdiffs=0(match), rpdiffs=0(match) 
AAGGGCAGGGTGACAGCGTCCGGGTAAATACGGAGGATGCAAGCGTTAATCGGAATTACTGGGCGTAAAGNGCTGTTAGGTGNTTTGCTAAGTCGANNTNTGAAAGCNCCGGGCTTAACCTAGNAAATGCATATGAACTGGCAAGCTTGAGTACAGTAGNGGGTGGCGGAATTTCCGGTGTAGNGGTGAAATGCGTAGAGATGGNAAGGAACATCAGTNGCGAAGGCGGCCACCTGGACTNATACTGACACTGAGGGACGAAAGCNAGGGTAGCGAANAGGATTAAGAAACCCGATAGTCCCTGTCCTTT
>M00763_26_000000000-K4TML_1_1101_8044_1728	ee=3.85617	fbdiffs=0(match), rbdiffs=0(match) fpdiffs=0(match), rpdiffs=0(match) 
GTGCCAACCGCCGCGGGAAAAAGGGAGGGGCAAGCGTTATNCGGCATAACTGGGCGTAAAGAGTNCGTAGACGGTAAAGTAAGTTTTTTGTTAAATTGTAAACCTTAANTTTAAAACNAGCATTAAATACTGCTTTACTTTGAGTTTAGTACAGAAAAGTAGAATTTTATATGGAAGGGTGAAATCTGCTAATATATAAAGGAATGTCATTTAGCGAAGGCGACTTTTTAGTATAAACTGACGTTGAGGGACGAAAGTGTGGGTATCGAACAGGATTAGATACCCCAGTAGTCC

```

It looks like the primers were removed successfully.  


### <a name="QC"></a> **4. QC'ing sequences with screen.seqs**    

Now, based on the output above, we need to remove "bad" sequences from the dataset. In the `screen.seqs()` function, we will specify the fasta file of the contigs generated in the previous step and remove any sequence with an ambiguous call ("N"). We will also remove sequences >350 nt. We will also set a minimum size amount (200). These parameters could be adjusted based on specific experiment and variable region.  

The code will look like this:  

```
screen.seqs(fasta=mcap.trim.contigs.fasta, group=mcap.contigs.groups, maxambig=0, maxlength=300, minlength=250)
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

This removed ~70,000 sequences out of 317,885 total.  

The summary output as viewed in the `screen_output_script` file now reads: 

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	204     204     0	3	1
2.5%-tile:	1	292     292     0	4	6178
25%-tile:	1	292     292     0	4	61778
Median:         1	292     292     0	4	123556
75%-tile:	1	301     301     0	5	185334
97.5%-tile:     1	310     310     0	6	240934
Maximum:        1	348     348     0	26	247111
Mean:   1	295     295     0	4
# of Seqs:	247111

```

We now see that we have removed all sequences with ambigous calls and the max sequence length is <300. Note that in this dataset we have sequences longer than expected for V4. In steps below we will align to a reference V4.  

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

In this run, there were 247,111 sequences and 176,977 were unique = ~70% (seen at the end of the output).  

The output table from summary.seqs looks like this and shows the number of unique and the total: 

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	204     204     0	3	1
2.5%-tile:	1	292     292     0	4	6178
25%-tile:	1	292     292     0	4	61778
Median:         1	292     292     0	4	123556
75%-tile:	1	301     301     0	5	185334
97.5%-tile:     1	310     310     0	6	240934
Maximum:        1	348     348     0	26	247111
Mean:   1	295     295     0	4
# of unique seqs:	176977
total # of seqs:        247111

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

The summary now looks like this:  

```
nano align_output_script 
```

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        0	0	0	0	1	1
2.5%-tile:	1	13424   292     0	3	4425
25%-tile:	1	13424   292     0	4	44245
Median:         1	13424   292     0	4	88489
75%-tile:	1	13424   292     0	5	132733
97.5%-tile:     1	13425   293     0	6	172553
Maximum:        13425   13425   301     0	18	176977
Mean:   53	13373   289     0	4
# of Seqs:	176977

```

From this, we see that 176,977 sequences aligned to the reference, which matches the number of sequences that we had after the unique.sh step (176,977). In the next steps we will filter out any sequences that don't meeting alignment settings.  

#### QC sequences according to alignment to the reference  

Our sequences now align at the correct positions on the reference.

Now remove sequences that are outside the alignment window (1968-11550bp). This removes anything that starts after `start` and ends before `end`. Maxhomop=8 argument removes anything that has repeats greater than the threshold - e.g., 8 A's in a row = polymer 8. Here we will removes polymers >8 because we are confident these are likely not "real".  

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

The summary looks like this: 

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	11552   271     0	3	1
2.5%-tile:	1	13424   292     0	4	6130
25%-tile:	1	13424   292     0	4	61296
Median:         1	13424   292     0	4	122592
75%-tile:	1	13424   292     0	5	183888
97.5%-tile:     1	13425   293     0	6	239054
Maximum:        1968    13425   300     0	8	245183
Mean:   1	13424   292     0	4
# of unique seqs:	175479
total # of seqs:        245183
```


We have now "tagged" sequences outside of the window of interest in our alignment to the Silva reference 16S V4 region. In the next step we will actually remove them.     

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
Length of filtered alignment: 513
Number of columns removed: 12912
Length of the original alignment: 13425
Number of sequences used to construct filter: 175479

```

We also get a new summary that looks like this:  

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       511     242     0	3	1
2.5%-tile:      1       513     254     0	3	6130
25%-tile:	1	513     254     0	4	61296
Median:         1       513     254     0	4	122592
75%-tile:       1       513     254     0       4       183888
97.5%-tile:     1	513     254     0	6	239054
Maximum:        1       513     272     0       8	245183
Mean:   1	512     254     0	4
# of unique seqs:       175479
total # of seqs:        245183

```

From this summary we see that the alignment window spans ~500 bp and the length of our sequences is about 254 nt. We have a maximum polymer of 8 as specified in our settings above.  

### <a name="Precluster"></a> **7. Polish the data with pre clustering**     

Now we need to further polish and cluster the data with pre.cluster. In Mothur, pre-clustering is done to obtain ASV, or amplicon sequence variants. This will serve as the input for identifying OTUs in future steps. 

The purpose of this step is to remove noise due to sequencing error. The rational behind this step assumes that the most abundant sequences are the most trustworthy and likely do not have sequencing errors. Pre-clustering then looks at the relationship between abundant and rare sequences - rare sequences that are "close" (e.g., 1 nt difference) to highly abundant sequences are likely due to sequencing error. This step will pool sequences and look at the maximum differences between sequences within this group to form ASV groupings. 

In this step, the number of sequences is not reduced, but they are grouped into ASV's which reduces the error rate. V4 region has the lowest likelihood of errors, so the error rate is going to be lower than for other variable regions.  

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

The program identified 25,035 unique sequences after looking at 175,479 sequences. After identifying unique sequences, Mothur output the following files:   

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

Finally, we get the output from summary:   

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	511     242     0	3	1
2.5%-tile:	1	513     254     0	3	6130
25%-tile:	1	513     254     0	4	61296
Median:         1	513     254     0	4	122592
75%-tile:	1	513     254     0	4	183888
97.5%-tile:     1	513     254     0	6	239054
Maximum:        1	513     272     0	8	245183
Mean:   1	512     254     0	4
# of unique seqs:	10918
total # of seqs:        245183
```

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
Minimum:        1	511     242     0	3	1
2.5%-tile:	1	513     254     0	3	6057
25%-tile:	1	513     254     0	4	60561
Median:         1	513     254     0	4	121121
75%-tile:	1	513     254     0	4	181681
97.5%-tile:     1	513     254     0	6	236185
Maximum:        1	513     272     0	8	242241
Mean:   1	512     254     0	4
# of unique seqs:	9045
total # of seqs:        242241
```


The program identified and removed 2,953 chimeras out of 245,183 sequences = 1.2% chimeras.    

We can look at a count of the number of sequences per sample: 

```
WSH174 contains 10909.
WSH175 contains 1598.
WSH176 contains 7156.
WSH177 contains 27867.
WSH178 contains 6468.
WSH179 contains 45074.
WSH180 contains 829.
WSH181 contains 633.
WSH182 contains 205.
WSH183 contains 106.
WSH184 contains 292.
WSH185 contains 22387.
WSH186 contains 1255.
WSH187 contains 8981.
WSH188 contains 3825.
WSH189 contains 158.
WSH190 contains 567.
WSH191 contains 250.
WSH192 contains 140.
WSH193 contains 3127.
WSH194 contains 224.
WSH195 contains 3440.
WSH196 contains 268.
WSH201 contains 418.
WSH202 contains 436.
WSH203 contains 2358.
WSH204 contains 1858.
WSH205 contains 3072.
WSH206 contains 4849.
WSH207 contains 6582.
WSH208 contains 13055.
WSH209 contains 10465.
WSH210 contains 8781.
WSH211 contains 13718.
WSH212 contains 5636.
WSH213 contains 4384.
WSH214 contains 13592.
WSH215 contains 1744.
WSH216 contains 5534.

Size of smallest group: 106.

```

The smallest group has 106 sequences. We will further look at this sampling depth.  

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
remove.lineage(fasta=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
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

mothur "#remove.lineage(fasta=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)"

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

*Removing only Chloroplast* 
Removed 1088 sequences from your fasta file.
Removed 43819 sequences from your count file.

*Removing only Mitochondria* 
Removed 2 sequences from your fasta file.
Removed 4 sequences from your count file.

*Removing only unknown domain* 
Removed 91 sequences from your fasta file.
Removed 694 sequences from your count file.

*Removing only Archaea* 
Removed 16 sequences from your fasta file.
Removed 39 sequences from your count file.

*Removing only Eukaryotes* 
No contaminants to remove.  

*Removing all lineages*  
Removed 1197 sequences from your fasta file.
Removed 44556 sequences from your count file.

We also can see the following summary:  

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	511     242     0	3	1
2.5%-tile:	1	513     254     0	3	4943
25%-tile:	1	513     254     0	4	49422
Median:         1	513     254     0	4	98843
75%-tile:	1	513     254     0	4	148264
97.5%-tile:     1	513     254     0	6	192743
Maximum:        1	513     272     0	8	197685
Mean:   1	512     254     0	4
# of unique seqs:	7848
total # of seqs:        197685
```

We now have 197685 total sequences and 7848 unique sequences.  

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

This is the end of the pipeline for curating our sequencing - we have corrected for pcr errors, removed chimeras, and removed sequences outside of taxon of interest. These are ASV's. Now we can move onto OTU clustering!

View a sample of the taxonomy summary dataset. 
 
```
head mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy
```

```
M00763_26_000000000-K4TML_1_1104_16138_15614	Bacteria(100);"Proteobacteria"(92);"Proteobacteria"_unclassified(92);"Proteobacteria"_unclassified(92);"Proteobacteria"_unclassified(92);"Proteobacteria"_unclassified(92);
M00763_26_000000000-K4TML_1_1107_28903_13620	Bacteria(100);"Proteobacteria"(100);Gammaproteobacteria(100);Alteromonadales(100);Alteromonadaceae(100);Alteromonadaceae_unclassified(100);
M00763_26_000000000-K4TML_1_1105_21427_2529	Bacteria(100);"Proteobacteria"(100);Gammaproteobacteria(100);Oceanospirillales(100);Oceanospirillaceae(100);Neptuniibacter(100);
M00763_26_000000000-K4TML_1_1101_20650_23111	Bacteria(100);"Proteobacteria"(99);Gammaproteobacteria(97);Oceanospirillales(82);Oceanospirillales_unclassified(82);Oceanospirillales_unclassified(82);
M00763_26_000000000-K4TML_1_1101_12790_23130	Bacteria(100);"Proteobacteria"(92);"Proteobacteria"_unclassified(92);"Proteobacteria"_unclassified(92);"Proteobacteria"_unclassified(92);"Proteobacteria"_unclassified(92);
M00763_26_000000000-K4TML_1_1102_3812_14988	Bacteria(100);"Bacteroidetes"(100);"Sphingobacteria"(98);"Sphingobacteriales"(98);"Flammeovirgaceae"(87);"Flammeovirgaceae"_unclassified(87);
M00763_26_000000000-K4TML_1_1119_8635_14800	Bacteria(100);"Proteobacteria"(89);Gammaproteobacteria(81);Gammaproteobacteria_unclassified(81);Gammaproteobacteria_unclassified(81);Gammaproteobacteria_unclassified(81);
M00763_26_000000000-K4TML_1_1101_16613_7668	Bacteria(100);"Proteobacteria"(89);"Proteobacteria"_unclassified(89);"Proteobacteria"_unclassified(89);"Proteobacteria"_unclassified(89);"Proteobacteria"_unclassified(89);
M00763_26_000000000-K4TML_1_2112_18154_7641	Bacteria(100);"Proteobacteria"(87);Gammaproteobacteria(86);Gammaproteobacteria_unclassified(86);Gammaproteobacteria_unclassified(86);Gammaproteobacteria_unclassified(86);
M00763_26_000000000-K4TML_1_2112_5217_7648	Bacteria(100);Bacteria_unclassified(100);Bacteria_unclassified(100);Bacteria_unclassified(100);Bacteria_unclassified(100);Bacteria_unclassified(100);
```






DO FOR ASV'S?????? 


### <a name="Cluster"></a> **10. Cluster for OTUs**  

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

Third, we will run `cluster.split` to calculate the distance and cluster within each taxonomic level (order in this case) in parallel, then bring the data back together. This will generate a file that indicates which sequence are in which OTU's. 

cluster.split uses a cutoff of 0.03 which is 3% difference.   

```
cluster.split(fasta=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, taxlevel=4, cutoff=0.03, splitmethod=classify)
```

Next we will make a shared file. This .shared file has the label for 3% OTU, sample name, the number of OTU's and then the number of time each OTU appears in each sample. 

This file will be the basis of what we will do to measure richness of communities compared to each other.  

We want to keep the shared file and a consensus taxonomy file. 

```
make.shared(list=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
```

The classify.otu command will then output a concensus cons.taxonomy file that has the taxonomic information for each OTU.   

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

mothur "#cluster.split(fasta=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=mcap.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, taxlevel=4, cutoff=0.03, splitmethod=classify)"

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

Cluster.split outputs the following files: 

```
Output File Names:
mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist
mcap.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list
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
WSH174 contains 9311.
WSH175 contains 1246.
WSH176 contains 5397.
WSH177 contains 26221.
WSH178 contains 4287.
WSH179 contains 28273.
WSH180 contains 750.
WSH181 contains 601.
WSH182 contains 193.
WSH183 contains 95.
WSH184 contains 284.
WSH185 contains 12547.
WSH186 contains 880.
WSH187 contains 3893.
WSH188 contains 1556.
WSH189 contains 128.
WSH190 contains 554.
WSH191 contains 243.
WSH192 contains 136.
WSH193 contains 3119.
WSH194 contains 223.
WSH195 contains 3433.
WSH196 contains 267.
WSH201 contains 416.
WSH202 contains 435.
WSH203 contains 2319.
WSH204 contains 1854.
WSH205 contains 3059.
WSH206 contains 4838.
WSH207 contains 6520.
WSH208 contains 12997.
WSH209 contains 10435.
WSH210 contains 8710.
WSH211 contains 13555.
WSH212 contains 5614.
WSH213 contains 4315.
WSH214 contains 12831.
WSH215 contains 1627.
WSH216 contains 4523.

Size of smallest group: 95.

Total seqs: 197685

```

*Now we have classified taxonomic data and we are ready to look at ecological statistics!*  


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

It may be helpful to run this with multiple iterations. If you do not include a size=# argument, then mothur will automatically set to the lowest sequence. If I set this to 300, for example the output shows that several samples are removed because there are <300 samples: 

```
WSH180 contains 750. Eliminating.
WSH181 contains 601. Eliminating.
WSH182 contains 193. Eliminating.
WSH183 contains 95. Eliminating.
WSH184 contains 284. Eliminating.
WSH186 contains 880. Eliminating.
WSH189 contains 128. Eliminating.
WSH190 contains 554. Eliminating.
WSH191 contains 243. Eliminating.
WSH192 contains 136. Eliminating.
WSH194 contains 223. Eliminating.
WSH196 contains 267. Eliminating.
WSH201 contains 416. Eliminating.
WSH202 contains 435. Eliminating.
Sampling 1000 from each group.
```  

Sample names are found in the [spreadsheet here](https://docs.google.com/spreadsheets/d/1lLvCp-RoRiBSGZ4NBPwi6cmZuozmfS20OJ7hBIueldU/edit#gid=0). QC information from extraction can be found in the [spreadsheet here](https://docs.google.com/spreadsheets/d/1Ew0AOs88n1i6wEyuqbv3tsRXwMjhqyYEeZX5sACFdDY/edit#gid=0).   

WSH180 = Attached Recruit (231 hpf)
WSH181 = TP1
WSH182 = TP1
WSH183 = TP1
WSH184 = TP1
WSH186 = TP10 (metamorphosed recruit 231 hpf)
WSH189 = TP2
WSH190 = TP2
WSH191 = TP2
WSH192 = TP2
WSH194 = TP3
WSH196 = TP3
WSH201 = TP4
WSH202 = TP4

This sub-sampling to 1000 sequences removes all of TP1 (eggs) and TP2 (embryos 5 hpf). We remove 2/4 of TP3 and TP4 (embryos 38 hpf and embryos 65 hpf). We also remove one metamorphosed recruit and one attached recruit.  

We know that DNA concentrations were very low for the TP1-TP3 samples so we may need to remove them from the analysis or would have to consider that the microbial communities will be different due to sequencing depth and DNA extraction artifacts.    

We did not have too low of DNA for the other lifestages, however, so we need to think about how to deal with these.  

If we need to standardized sampling depth for all of our samples, then we would need to remove these samples. 

#### Subsampling the sample set  

We are going to generate a rarefaction. Calculating the observed number of OTUs using the Sobs metric (observed number of taxa) with a freq=100 to output the data for every 100 sequences samples - if you did every single sample that would be way too much data to output.  

```
rarefaction.single(shared=mcap.opti_mcc.shared, calc=sobs, freq=100)
```

```
Output File Names: 
mcap.opti_mcc.groups.rarefaction
```

You can look at this file with `nano mcap.opti_mcc.groups.rarefaction`. It contains the number of OTUs detected for every 100 sequences for each sample. Remember OTU = 0.03 cut off here, so that is the OTU distinction. Series will be trucated after the number of sequences available in the sample.  

Next, rarefy the data to the minimum number of sequences, then the next command will rarefy and pull out that many sequences from each sample. If you want to sample a specific number, use subsample=#. You would want to do this if you have a super small sequence number in one group and you dont want to truncate the data for all groups by that much. If you leave subsample=T, mothur will calculate to the minimum size. Subsamples 1000 times and calculates the metrics we want with the number of otus. Sampling is random.   

All of the available metric calculations (calc) are [found here](https://mothur.org/wiki/calculators/).  

I am going to use subsample=1000 to remove the low sequence count stages.  

```
summary.single(shared=mcap.opti_mcc.shared, calc=nseqs-sobs-shannon-invsimpson, subsample=1000)
```

```
Output File Names: 
mcap.opti_mcc.groups.ave-std.summary
```

Open this file to view statistics for each sample. The method column has ave (average of each of the metrics across the 1000 iteractions) or std (this is the std deviation of the average calculation).  

We will clearly get different answers depending on how we subsample. 

For now, we can proceed with and without subsampling and then re run analyses with different settings to see how our answer changes.  

*I recommend running iterations of these analyses with different thresholds for the commands above to see how our data filtering is changed. Then we can make data-driven decisions.*    

### <a name="Statistics"></a> **12. Calculate Ecological Statistics**  

Alpha diversity = diversity within a sample 
Beta diversity = diversity between samples

We will calculate beta diversity metrics which will then be used for ecological analyses.  

#### Calculate Beta Diversity  

We will use the `dist.shared` function to calculate distances between samples. Here we are going to use Bray-Curtis distances. Many different metrics can be used (e.g., Theta YC). All of the available metric calculations (calc) are [found here](https://mothur.org/wiki/calculators/).  

We can use this function with or without `subsample=300` or `subsample=T`. We will use the bray-curtis distances for this calculation. This step is important for generating files that we will use later in R. 

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

```
mcap.opti_mcc.braycurtis.0.03.lt.dist
mcap.opti_mcc.braycurtis.0.03.lt.ave.dist
``` 

*After this point, analyses can be run in R as well. See the bottom of this document for transferring files to your computer.*   

#### Run a PCoA and NMDS   

Run a PCoA on the Bray Curtis distances without subsampling.    

```
pcoa(phylip=mcap.opti_mcc.braycurtis.0.03.lt.dist)
```

```
Processing...
Rsq 1 axis: 0.3997
Rsq 2 axis: 0.615637
Rsq 3 axis: 0.597504

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
Rsq 1 axis: 0.507001
Rsq 2 axis: 0.544563
Rsq 3 axis: 0.72678

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
Processing Dimension: 2
1
2
3
4
5
6
7
8
9
10

Number of dimensions:	2
Lowest stress :	0.314127
R-squared for configuration:	0.686885

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
Processing Dimension: 2
1
2
3
4
5
6
7
8
9
10

Number of dimensions:	2
Lowest stress :	0.278464
R-squared for configuration:	0.705154

Output File Names: 
mcap.opti_mcc.braycurtis.0.03.lt.ave.nmds.iters
mcap.opti_mcc.braycurtis.0.03.lt.ave.nmds.stress
mcap.opti_mcc.braycurtis.0.03.lt.ave.nmds.axes
``` 

If you run NMDS each time, you will get the same clustering/points but they might be rotated on different axes.  

Stress value was <1 so that is good. 

The iters file will also have information on stress values for each time it ran the tests.

R squared values for NMDS are slightly lower with subsampling.  

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
SS	8.11866	6.70995	14.8286
df	10	28	38
MS	0.811866	0.239641

Fs:	3.38784
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
SS	5.11832	2.51575	7.63407
df	8	16	24
MS	0.63979	0.157234

Fs:	4.06903
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
EggFertilized-Embryo1-Larvae1-Larvae2-Larvae3-Larvae4-Larvae6-Recruit1-Recruit1plug-Recruit2-Recruit2plug	-na<0.001*	0.255767	0.203303	0.300343	0.221489	0.128419	0.0715921	0.287117	-nan	0.328122	0.262186	0.39856

mcap.opti_mcc.braycurtis.0.03.lt.homova
```

There is a significant difference in variation in microbiomes between lifestages without subsampling. 


Test with subsampling.  

```
homova(phylip=mcap.opti_mcc.braycurtis.0.03.lt.ave.dist, design=mcap_design.txt)
```

```
HOMOVA	BValue	P-value	SSwithin/(Ni-1)_values
Larvae1-Larvae2-Larvae3-Larvae4-Larvae6-Recruit1-Recruit1plug-Recruit2-Recruit2plug	-nan	<0.001*	0.0636349	0.0548249	0.0643005	0.0728028	0.254806	-nan	0.261264	0.19983	0.23794

mcap.opti_mcc.braycurtis.0.03.lt.ave.homova
```

Note there are NA's due to removal of egg and embryo groups.  

There is a significant difference in variation in microbiomes between lifestages with subsampling as well.  

### <a name="Population"></a> **13. Population Level Analyses**    

There are several tools in mothur we can use to run further analyses on OTU's. Note that we can also do these steps in R, which is what I do with this data.  

#### metastats  

We can now use a couple approaches to look for differences in OTU's by group. If we only have two groups, we can use the `metastats` function. For example:  

```
metastats(shared=mcap.opti_mcc.0.03.subsample.shared, design=mcap_design.txt)
``` 

However, because I have >2 groups for this dataset, I need to use something else - this generates a file for each pairwise comparison. We will instead use the `lefse` function.  

#### lefse analysis  

Without subsampling: mcap.opti_mcc.shared

```
lefse(shared=mcap.opti_mcc.shared, design=mcap_design.txt)
```


With subsampling: mcap.opti_mcc.0.03.subsample.shared

```
lefse(shared=mcap.opti_mcc.0.03.subsample.shared, design=mcap_design.txt)
```

I had errors in this analysis, skipping for now for analysis in R.  


### <a name="Output"></a> **14. Output data for analysis in R**  

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

Outside Andromeda do the following for each file:  

```
scp ashuffmyer@ssh3.hac.uri.edu:/data/putnamlab/ashuffmyer/AH_MCAP_16S/mothur/mcap.opti_mcc.shared ~/MyProjects/EarlyLifeHistory_Energetics/Mcap2020/Output/16S/mothur

```