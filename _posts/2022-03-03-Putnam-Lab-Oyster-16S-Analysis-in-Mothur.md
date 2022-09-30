---
layout: post
title: Putnam Lab Oyster 16S Analysis in Mothur
date: '2022-03-03'
categories: PutnamLab_Oysters
tags: 16S Molecular Bioinformatics
---
This post details 16S analysis of the Putnam Lab oyster project using the Mothur pipeline based on [A. Huffmyer's pipeline here](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/16S-Analysis-in-Mothr-Part-1/).   

# **16S analysis workflow in mothur**  

This workflow details the analysis of the V6 region data set incorporating all sample types and subsetting the target samples (gut) after the pipeline.  

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

First, we need to prepare a directory in Andromeda to run our analyses.  

I tried to create a shared folder for our analysis but did not have permission. So I am creating a directory in my own files for now.  

```
ssh -l ashuffmyer ssh3.hac.uri.edu
cd /data/putnamlab/ashuffmyer/
mkdir oyster_16S
mkdir v6
cd v6
```  

Give permissions for all to read, write, and execute.   

`chmod u=rwx,g=rwx,o=rwx,a=rwx -R v6`  

Next, create symlinks to the files that we will need for analysis.  

```
ln -s /data/putnamlab/shared/PointJudithData_Rebecca/amplicons16s/allsamples_V6/*.txt /data/putnamlab/ashuffmyer/oyster_16S/v6

ln -s /data/putnamlab/shared/PointJudithData_Rebecca/amplicons16s/allsamples_V6/*.csv /data/putnamlab/ashuffmyer/oyster_16S/v6

ln -s /data/putnamlab/shared/PointJudithData_Rebecca/amplicons16s/allsamples_V6/00_RAW_gz/*.fastq.gz /data/putnamlab/ashuffmyer/oyster_16S/v6
```

This creates links to the original location of the data without the need to copy/download the data.  

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

This creates a file called `oyster.files` with sample name, R1, R2 listed in columns.  

Then we will align sequences together for R1 and R2 to generate contigs. 

We will write a bash script to run the `make.contigs()` command and make contigs of forward and reverse reads for each sample. This will assemble contigs for each pair of files and use R1 and R2 quality information to correct any sequencing errors in locations where there is overlap and the quality of the complementary sequence is sufficient. 

We use `trimoverlap=T` in the `make.contigs` step. The region is ~55bp in length and we used 2x75bp sequencing to obtain overlap.   

We will also use and oligos file to remove primers `oligos=oligos.oligos` specifying `pdiffs=2` and `checkorient=T`.  

We will also run the `summary.seqs()` command with the `make.file()` and `make.contig()` steps. This generates summary information about the sequences from the files we made above.  

Within the `v6` directory, create a script. Everything will live within this directory - I did not make any subfolders.  

Before we run the script we need to make an `oligos` file. This file contains the primers for our sequences. Mothur will remove these primers while building the contigs.  

The primers that we have are:  

```
Huber et al. 2007

967F: TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG CTAACCGANGAACCTYACC
      TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG CNACGCGAAGAACCTTANC
      TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG CAACGCGMARAACCTTACC
      TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG ATACGCGARGAACCTTACC
      
1046R:GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG CGACRRCCATGCANCACCT
```

The illumina adapter forward sequence is `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`. The last 19 nt's of the forward primers are the primers we need to remove. The reverse adapter sequence is `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG`. The last 19 nt's is the reverse primer.   

Our primers are:  
```
FORWARD
CTAACCGANGAACCTYACC
CNACGCGAAGAACCTTANC
CAACGCGMARAACCTTACC
ATACGCGARGAACCTTACC

REVERSE
CGACRRCCATGCANCACCT
```

For example, search for primers in a forward sequence file using the following:     
`zcat RS181_S11_L001_R1_001.fastq.gz | grep -c "NNNN"`    

*Need to replace degenerate bases with complement*  
* Y = C or T  
* N = Remove if on the end - if in the middle, replace with ACTG   
* R = A or G 
* M = A or C  

`zcat RS181_S11_L001_R1_001.fastq.gz | grep -c "CTAACCGAAGAACCTTACC"` 
Found 1991 times.  

`zcat RS181_S11_L001_R1_001.fastq.gz | grep -c "CAACGCGAAGAACCTTACC"`
Found 4638 times.  

We will check for these primers at the end of our contigs step.  

```
nano oligos.oligos
```

```
primer CTAACCGANGAACCTYACC CGACRRCCATGCANCACCT
primer CNACGCGAAGAACCTTANC CGACRRCCATGCANCACCT
primer CAACGCGMARAACCTTACC CGACRRCCATGCANCACCT
primer ATACGCGARGAACCTTACC CGACRRCCATGCANCACCT

```


```
nano contigs.sh
``` 

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
#SBATCH --account=putnamlab  

echo "Job started" $(date)  

module load Mothur/1.46.1-foss-2020b

module list

mothur

mothur "#make.file(inputdir=., type=gz, prefix=oyster)"

mothur "#make.contigs(inputdir=., outputdir=., file=oyster.files, trimoverlap=T, oligos=oligos.oligos, pdiffs=2, checkorient=T)"

mothur "#summary.seqs(fasta=oyster.trim.contigs.fasta)"

echo "Job ended" $(date)  

```

Note the syntax for running mothur commands in a slurm script. The syntax is to activate mothur with `mothur` followed by the command in quotes and a hashtag (`"#make.files()"`).  

Run script and check status (takes ~2hrs).      

```
sbatch contigs.sh

squeue -u ashuffmyer -l
```

To view the output of each script, you can view with the following:  

```
nano contigs_output_script
```

*For all of these steps in mothur, you will be looking at these script output files to see any output and lists of any files that were created. I recommend copying these outputs so that you can track the output files.*   

The following files were output from the `make.files()` step: 
 
```
oyster.files
```

The following files were output from the `make.contigs()`step:  

```
oyster.trim.contigs.fasta
oyster.scrap.contigs.fasta
oyster.contigs.report
oyster.contigs.groups
```

Descriptions of contig files:  
*Trim file* = sequences that were "good".  
*Scrap file* = sequences that were "bad".  
*Groups file* = what group each sequence belongs to map sequence to each sample from the trimmed sequence file.  
*Contigs report file* = information on sequences that were aligned and paired together. 

Count the number of sequences that were removed and the number that were kept by counting sequences in the output fasta files.  

```
grep -c "^>" oyster.trim.contigs.fasta
grep -c "^>" oyster.scrap.contigs.fasta
```
5,179,429 sequences were kept and 270,837 were removed.  

Next, check for the presence of primers in the output.  

For example, check for a couple combinations of primers: 

`grep -c "CTAACCGAAGAACCTTACC" oyster.trim.contigs.fasta` 
Found 0 times.  

`grep -c "CAACGCGAAGAACCTTACC" oyster.trim.contigs.fasta`
Found 0 times. 

Success! Our primers are removed.  

View the output.  

```
less contigs_output_script
```

You will get an output table like the one below and it will output a file with this summary information.  

Scroll to the bottom to find the summary information.  


``` 
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	1	1	0	1	1
2.5%-tile:	1	48	48	0	2	129486
25%-tile:	1	52	52	0	3	1294858
Median:         1	54	54	0	3	2589715
75%-tile:	1	54	54	0	4	3884572
97.5%-tile:     1	56	56	1	6	5049944
Maximum:        1	60	60	24	14	5179429
Mean:   1	52	52	0	3
# of Seqs:	5179429
``` 

This table shows quantile values about the distribution of sequences for a few things:  

- *Start position*: All at 1 now, will start at different point after some QC.   
- *End position*: We see that there are some sequences that are very short and we may need to remove those later.   
- *Number of bases*: length (we see most are in expected range here, but one is super long! This might tell us there is no overlap so they are butted up against each other. We will remove things like this.  
- *Ambigs*: Number of ambiguous calls in sequences. Here there are a few that have ambiguous base calls. We will remove any sequence with an ambiguous call or any longer than we would expect for V4 region.   
- *Polymer*: Length of polymer repeats.    
- *NumSeqs*: Number of sequences.  


### <a name="QC"></a> **4. QC'ing sequences with screen.seqs**    

Now, based on the output above, we need to remove "bad" sequences from the dataset. In the `screen.seqs()` function, we will specify the fasta file of the contigs generated in the previous step and remove any sequence with an ambiguous call ("N"). We will also remove sequences outside the length expected. We will also set a minimum size amount. These parameters could be adjusted based on specific experiment and variable region.  

We are optimizing this for the V6 region.  

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

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#screen.seqs(inputdir=., outputdir=., fasta=oyster.trim.contigs.fasta, group=oyster.contigs.groups, maxambig=0, maxlength=60, minlength=40)"

mothur "#summary.seqs(fasta=oyster.trim.contigs.good.fasta)"

```

```
sbatch screen.sh
```

This generates the following output files: 
`oyster.trim.contigs.good.fasta`
`oyster.trim.contigs.bad.accnos`
`oyster.contigs.good.groups`

"good" seqs satisfied criteria and "bad" seqs did not meet criteria. 







Count the number of "good" and "bad" sequences.  

```
grep -c "^>" oyster.trim.contigs.good.fasta
grep -c ".*" oyster.trim.contigs.bad.accnos
```

4,977,936 sequences were kept and 201,493 were removed. 

This removed ~4% of sequences due to length and ambiguous bases.  

You can also view the `bad.accnos` file to see why sequences were removed.  

```
head -1000 oyster.trim.contigs.bad.accnos
``` 

The summary output as viewed by `less screen_output_script` file now reads: 

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	40	40	0	2	1
2.5%-tile:	1	48	48	0	2	124449
25%-tile:	1	52	52	0	3	1244485
Median:         1	54	54	0	3	2488969
75%-tile:	1	54	54	0	4	3733453
97.5%-tile:     1	56	56	0	6	4853488
Maximum:        1	59	59	0	14	4977936
Mean:   1	53	53	0	3
# of Seqs:	4977936

```

We now see that we have removed all sequences with ambigous calls and the max sequence length is 59 with min 40. 

### <a name="Unique"></a> **5. Determining and counting unique sequences**  

Next, determine the number of unique sequences. 

After determining the unique sequences, we can use these unique sequences to count how many times each sequence (which will later be classified to OTU) shows up in each sample.  

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

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#unique.seqs(fasta=oyster.trim.contigs.good.fasta)"

mothur "#count.seqs(name=oyster.trim.contigs.good.names, group=oyster.contigs.good.groups)"

mothur "#summary.seqs(fasta=oyster.trim.contigs.good.unique.fasta, count=oyster.trim.contigs.good.count_table)"

```

```
sbatch unique.sh
```

This script generate the following outputs: 

```
oyster.trim.contigs.good.names
oyster.trim.contigs.good.unique.fasta
oyster.trim.contigs.good.count_table
```

```
nano unique_output_script
``` 

In this run, there were 4,977,936 sequences and 133,806 were unique = ~3%.  

The output table from `summary.seqs` looks like this and shows the number of unique and the total sequences: 

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	40	40	0	2	1
2.5%-tile:      1	48	48	0	2	124449
25%-tile:       1	52	52	0	3	1244485
Median:         1	54	54	0	3	2488969
75%-tile:       1	54	54	0	4	3733453
97.5%-tile:     1	56	56	0	6	4853488
Maximum:        1	59	59	0	14	4977936
Mean:   1       53	53	0	3
# of unique seqs:	133806
total # of seqs:        4977936

```

Now we can align just the unique sequences, which will be much faster than aligning the full data set.    

*From this, we have our unique sequences identified and can proceed with further cleaning and polishing of the data. Next we will look at alignment, error rate, chimeras, classification and further analysis.*  

### <a name="Align"></a> **6. Align our sequences to a reference database** 

#### Optimize the reference for the V6 region  

We need to find the coordinates for our sequences in the Silva database. We can do this by amplifying the E. coli 16S gene using our specific primers. We will then align this region of the E. coli gene to the Silva reference database. This will show us the coordinates that we need to use to align our references to the Silva database.  

More information at:  
https://mothur.org/blog/2016/Customization-for-your-region/  

First download a file for the E. coli 16S gene. I downloaded the fasta file from here: https://www.ncbi.nlm.nih.gov/nuccore/174375?report=fasta.    

Outside Andromeda, copy the file to our folder.  

```
scp ~/Desktop/ecoli_sequence.fasta ashuffmyer@ssh3.hac.uri.edu:/data/putnamlab/ashuffmyer/oyster_16S/v6 

```

Create a file that has our primers.  

```
nano pcrTest.oligos
```

```
primer CTAACCGANGAACCTYACC CGACRRCCATGCANCACCT
primer CNACGCGAAGAACCTTANC CGACRRCCATGCANCACCT
primer CAACGCGMARAACCTTACC CGACRRCCATGCANCACCT
primer ATACGCGARGAACCTTACC CGACRRCCATGCANCACCT
``` 

First, download the Silva reference files from the [Mothur Wiki](https://mothur.org/wiki/silva_reference_files/) at the latest release.

The silva reference is used and recommended by the Mothur team. It is a manually curated data base with high diversity and high alignment quality.  

In your mothur directory, run `wget` to download the silva reference and training sets from the mothur website and unzip them and move into the right directory. 

```
wget https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.bacteria.zip

wget https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset9_032012.pds.zip

unzip silva.bacteria.zip
cd silva.bacteria
cp silva.bacteria.fasta ../silva.bacteria.fasta
cd ../

unzip trainset9_032012.pds.zip
```

Now we have the reference files that we need in our directory.  

Start an interactive mode to align the E. coli sequence. We only need to do this step once and once we have the coordinates we do not need to run it again.    

```
interactive
module load Mothur/1.46.1-foss-2020b
mothur
```

Next we need to run `pcr.seqs` to align the E. coli gene using our primers.  

```
pcr.seqs(fasta=ecoli_sequence.fasta, oligos=pcrTest.oligos)
```

Next, align the sequences that we created from amplifying with our primers to the silva reference.  

```
align.seqs(fasta=ecoli_sequence.pcr.fasta, reference=silva.bacteria.fasta)
```

Next, view the location of the V6 region.  

```
summary.seqs(fasta=ecoli_sequence.pcr.align)
```

```
		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	31189	33183	60	0	4	1
2.5%-tile:	31189	33183	60	0	4	1
25%-tile:	31189	33183	60	0	4	1
Median: 	31189	33183	60	0	4	1
75%-tile:	31189	33183	60	0	4	1
97.5%-tile:	31189	33183	60	0	4	1
Maximum:	31189	33183	60	0	4	1
Mean:	31189	33183	60	0	4
# of Seqs:	1
```

We can see our coordinates are start=31189 and end=33183 with the NBases=60. Note that the positions will be longer because there are "dots" or open spaces in the silva reference.    

Exit interactive mode.  
`quit`  
`exit`  

#### Build the reference using the V6 region  

Now write a script to build a reference using our new coordinates.  

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

mothur "#pcr.seqs(fasta=silva.bacteria.fasta, start=31189, end=33183, keepdots=F)"

mothur "#summary.seqs(fasta=silva.bacteria.pcr.fasta)"

mothur "#rename.file(input=silva.bacteria.pcr.fasta, new=silva.v6.fasta)"

```

```
sbatch silva_ref.sh
```

The script outputs the following files: 

```
Output File Names: 
silva.bacteria.pcr.fasta
```

The file was renamed to `silva.v6.fasta` as specified in the script. 

The summary of sequences looks like this in the script output `less silva_ref_output_script`: 

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	1642    51	0	2	1
2.5%-tile:	1	1994    56	0	3	374
25%-tile:	1	1994    58	0	3	3740
Median:         1	1994    59	0	3	7479
75%-tile:	1	1994    61	0	4	11218
97.5%-tile:     1	1994    65	0	6	14583
Maximum:        4	1994    120     4	9	14956
Mean:   1	1993    59	0	3
# of Seqs:	14956

```

We now have a reference to align to for the V6 region!  


#### Align sequences to the reference  

We next align sequences with the `align.seqs` command.   

```
align.seqs(fasta=oyster.trim.contigs.good.unique.fasta, reference=silva.v6.fasta)

summary.seqs(fasta=oyster.trim.contigs.good.unique.align)
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

mothur "#align.seqs(fasta=oyster.trim.contigs.good.unique.fasta, reference=silva.v6.fasta)"

mothur "#summary.seqs(fasta=oyster.trim.contigs.good.unique.align)"

```

```
sbatch align.sh
``` 

The script generates these output files: 

```
Output File Names:
oyster.trim.contigs.good.unique.align
oyster.trim.contigs.good.unique.align.report
oyster.trim.contigs.good.unique.flip.accnos
```

View the report file.  

```
head oyster.trim.contigs.good.unique.align.report
```

View the accnos file.  

```
head oyster.trim.contigs.good.unique.flip.accnos
```

The summary now looks like this:  

```
less align_output_script 
```

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        0       0       0       0       1       1
2.5%-tile:      1       1610    37      0       2       3346
25%-tile:       5       1629    51      0       3       33452
Median:         5       1636    53      0       3       66904
75%-tile:       1209    1639    54      0       4       100355
97.5%-tile:     1221    1994    56      0       6       130461
Maximum:        1994    1994    58      0       14      133806
Mean:   384     1671    51      0       3
# of Seqs:      133806

```

From this, we see that 133,806 sequences aligned to the reference, which matches the number of sequences that we had after the unique.sh step (133,806). In the next steps we will filter out any sequences that don't meet alignment settings.  

#### QC sequences according to alignment to the reference  

Our sequences now align at the correct positions on the reference.

Now remove sequences that are outside the alignment window. This removes anything that starts after `start` and ends before `end`. We will use a start=500 and we will use end=1610. Maxhomop=8 argument removes anything that has repeats greater than the threshold - e.g., 8 A's in a row = polymer 8. Here we will removes polymers >8 because we are confident these are likely not high quality sequences (see mothur MiSeq SOP for more information).   

I selected the start/stop by trying different settings to see where our outliers are and using the percentiles in the summary files. There was a proportion starting at >500 start position with the rest starting at <10. This is how I arrived at the 500 cut off.    

This is the command we will use: 

```
screen.seqs(fasta=oyster.trim.contigs.good.unique.align, count=oyster.trim.contigs.good.count_table, start=500, end=1610, maxhomop=8)
```

We will then run a summary.  

```
summary.seqs(fasta=oyster.trim.contigs.good.unique.good.align, count=oyster.trim.contigs.good.good.count_table)
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

mothur "#screen.seqs(fasta=oyster.trim.contigs.good.unique.align, count=oyster.trim.contigs.good.count_table, start=500, end=1610, maxhomop=8)"

mothur "#summary.seqs(fasta=oyster.trim.contigs.good.unique.good.align, count=oyster.trim.contigs.good.good.count_table)"

```

```
sbatch screen2.sh
```

This will output the following files: 
```
Output File Names:
oyster.trim.contigs.good.unique.good.align
oyster.trim.contigs.good.unique.bad.accnos
oyster.trim.contigs.good.good.count_table
```

View the accnos file to see why sequences will be removed and count the number of "bad" sequences.  

```
head oyster.trim.contigs.good.unique.bad.accnos

grep -c ".*" oyster.trim.contigs.good.unique.bad.accnos
``` 

43,505 uniques are tagged to be removed due to filtering at this step. This totals to 1,182,121 total sequences removed (this number is in the output file) out of the 4,977,936 total (23%).     

The summary looks like this: 

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1	1610    46	0	2	1
2.5%-tile:	1	1625    50	0	3	94896
25%-tile:	4	1636    53	0	3	948954
Median:         5	1636    54	0	3	1897908
75%-tile:	5	1642    55	0	4	2846862
97.5%-tile:     7	1994    56	0	6	3700920
Maximum:        7	1994    58	0	8	3795815
Mean:   4	1717    53	0	3
# of unique seqs:	90301
total # of seqs:        3795815
```

#### Filter sequences  

Now we can filter out sequences that didn't meet our criteria above, which will generate a report and a new summary of our sequences.  

We will run the following code. We align vertically and use trump=. to align the sequences accounting for periods in the reference.   

```
filter.seqs(fasta=oyster.trim.contigs.good.unique.good.align, vertical=T, trump=.)
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

mothur "#filter.seqs(fasta=oyster.trim.contigs.good.unique.good.align, vertical=T, trump=.)"

mothur "#summary.seqs(fasta=oyster.trim.contigs.good.unique.good.filter.fasta, count=oyster.trim.contigs.good.good.count_table)"

```

```
sbatch filter.sh
```

This script outputs the following files:  

```
Output File Names: 
oyster.filter
oyster.trim.contigs.good.unique.good.filter.fasta
```

We get a report on the filtering in the script output file that looks like this: 

```
Length of filtered alignment: 140
Number of columns removed: 1854
Length of the original alignment: 1994
Number of sequences used to construct filter: 90301
```

We also get a new summary that looks like this:  

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       125     38      0       2       1
2.5%-tile:	1	132     43	0	3	94896
25%-tile:       1       135     44      0	3	948954
Median:         1	140     46	0	3	1897908
75%-tile:       1	140     46	0	4	2846862
97.5%-tile:     1	140     47	0	6	3700920
Maximum:        5       140     54	0	8	3795815
Mean:   1       137     45	0	3
# of unique seqs:       90301
total # of seqs:        3795815

```

From this summary we see that the alignment window spans ~140 bp and the length of our sequences is about ~50 nt. We have a maximum polymer of 8 as specified in our settings above.  

### <a name="Precluster"></a> **7. Pre cluster the sequences**     

Now we need to further polish and cluster the data with pre.cluster. The purpose of this step is to remove noise due to sequencing error. The rational behind this step assumes that the most abundant sequences are the most trustworthy and likely do not have sequencing errors. Pre-clustering then looks at the relationship between abundant and rare sequences - rare sequences that are "close" (e.g., 0 nt difference) to highly abundant sequences are likely due to sequencing error. We chose to precluster to 0 nt difference because with a short sequence we want to keep as much information as possible before we cluster to similarity by OTU. This step will pool sequences and look at the maximum differences between sequences within this group to form ASV groupings. 

In this step, the number of sequences is not reduced, but they are grouped into amplicon sequence variants ASV's which reduces the error rate.  

Other programs that conduct this "denoising" are DADA2, UNOISE, and DEBLUR. However, these programs remove the rare sequences, which can distort the relative abundance of remaining sequences. DADA2 also removes all sigletons (sequences with single representation) which disproportionately affects the sequence relative abundance. Mothur avoids the removal of rare sequences for this reason. 

We will first add code to identify unique sequences after the filtering steps above.  

```
unique.seqs(fasta=oyster.trim.contigs.good.unique.good.filter.fasta, count=oyster.trim.contigs.good.good.count_table)
```

We will then perform the pre-clustering a default of 1 nt difference. Diffs can be changed according to your requirements.  

```
pre.cluster(fasta=oyster.trim.contigs.good.unique.good.filter.unique.fasta, count=oyster.trim.contigs.good.unique.good.filter.count_table, diffs=1) 
```

Finally, we will run another summary.  

```
summary.seqs(fasta=oyster.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.count_table)
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

mothur "#unique.seqs(fasta=oyster.trim.contigs.good.unique.good.filter.fasta, count=oyster.trim.contigs.good.good.count_table)"

mothur "#pre.cluster(fasta=oyster.trim.contigs.good.unique.good.filter.unique.fasta, count=oyster.trim.contigs.good.unique.good.filter.count_table, diffs=0)"

mothur "#summary.seqs(fasta=oyster.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.count_table)"
```

```
sbatch precluster.sh
```

Mothur output the following files:   

```
Output File Names: 
oyster.trim.contigs.good.unique.good.filter.count_table
oyster.trim.contigs.good.unique.good.filter.unique.fasta
```

Then, the pre-clustering step outputs a file for each sample. The two most important files are:  

```
oyster.trim.contigs.good.unique.good.filter.unique.precluster.fasta
oyster.trim.contigs.good.unique.good.filter.unique.precluster.count_table
``` 

The other files have text for maps of sequence name, errors, abundance, differences, and the filtered sequence for each sample.   

Finally, we get the output from the summary:   

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       125     38      0       2       1
2.5%-tile:      1       132     43      0       3       94896
25%-tile:       1       135     44      0       3       948954
Median:         1       140     46      0       3       1897908
75%-tile:       1       140     46      0       4       2846862
97.5%-tile:     1       140     47      0       6       3700920
Maximum:        5       140     54      0       8       3795815
Mean:   1       137     45      0       3
# of unique seqs:       62039
total # of seqs:        3795815

```

Note that the number of unique sequences has decreased from 90,301 to 62,039 as expected since we are clustering sequences that are within 0 nt difference from each other (similar to ASV).  

### <a name="Chimera"></a> **8. Identify chimeras**  

Now we will remove chimeras using the dereplicate method. In this method, we are again using the assumption that the highest abundance sequences are most trustworthy. Chimeras are sequences that did not extend during PCR and then served as templates for other PCR products, forming sequences that are partially from one PCR product and partially from another. This program looks for chimeras by comparing each sequences to the next highest abundance sequences to determine if a sequence is a chimera of the more abundance sequences.  

We will use the `chimera.vsearch` function to identify chimeras: 

```
chimera.vsearch(fasta=oyster.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=T)
``` 

We will then remove the identified chimeras with `remove.seqs`:  

```
remove.seqs(fasta=oyster.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)
```

Finally, we will run a new summary:  

```
summary.seqs(fasta=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)
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

mothur "#chimera.vsearch(fasta=oyster.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=T)"

mothur "#remove.seqs(fasta=oyster.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)"

mothur "#summary.seqs(fasta=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)"

mothur "#count.groups(count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)"

```

```
sbatch chimera.sh
```


This script outputs the following files: 

```
Output File Names:
oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table
oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.chimeras
oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos
```

The new summary looks like this:

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       125     38      0       2       1
2.5%-tile:      1       132     43      0       3       94517
25%-tile:       1       135     44      0       3       945169
Median:         1       140     46      0       3       1890337
75%-tile:       1       140     46      0       4       2835505
97.5%-tile:     1       140     47      0       6       3686157
Maximum:        5       140     54      0       8       3780673
Mean:   1       137     45      0       3
# of unique seqs:       61489
total # of seqs:        3780673
```

The program identified and removed <1% chimeras.    

We can look at a count of the number of sequences per sample: 

```
RS187 contains 19215.
RS188 contains 28278.
RS189 contains 35919.
RS190 contains 38396.
RS191 contains 29307.
RS192 contains 48023.
RS193 contains 35852.
RS194 contains 31030.
RS195 contains 25906.
RS196 contains 29962.
RS197 contains 46391.
RS198 contains 38682.
RS199 contains 37317.
RS200 contains 31466.
RS201 contains 45964.
RS202 contains 24889.
RS203 contains 62435.
RS204 contains 43186.
RS205 contains 36019.
RS206 contains 48068.
RS207 contains 31736.
RS208 contains 34733.
RS209 contains 27081.
RS210 contains 32014.
RS211 contains 41452.
RS212 contains 45103.
RS213 contains 30833.
RS214 contains 37673.
RS215 contains 38619.
RS216 contains 34170.
RS217 contains 43151.
RS218 contains 41897.
RS219 contains 44560.
RS220 contains 38910.
RS221 contains 31660.
RS222 contains 32908.
RS223 contains 20861.
RS224 contains 24804.
RS225 contains 9477.
RS226 contains 39385.
RS227 contains 39599.
RS228 contains 39530.
RS229 contains 37025.
RS230 contains 37139.
RS231 contains 43080.
RS232 contains 33468.
RS233 contains 42194.
RS234 contains 31141.
RS235 contains 36737.
RS246 contains 38049.
RS247 contains 24967.

Size of smallest group: 9477.

Total seqs: 3780673.
```

The smallest group has 9,477 sequences. We will further look at this sampling depth and use this number for subsampling.  


### <a name="Classify"></a> **9. Classifying sequences**  

Now our sequences are clean and ready for classification!  

We will use the training set downloaded above from the silva database through the [Mothur wiki](https://mothur.org/wiki/classify.seqs/). 

#### Classify sequences  

We will use the `classify.seqs` command: 

```
classify.seqs(fasta=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds.tax)
```

The training and taxonomy files were from downloads. The mothur [classify.seq wiki page](https://mothur.org/wiki/classify.seqs/) has the reference and training sets that you can download and put into the directory to use for analyzing data. These sets include adding chlorophyll, mitchondria, etc to identify and remove these.     

The output file from `classify.seqs()` ending in .taxonomy has the name of sequence and the classification with % confidence in parentheses for each level. It will end at the level that is has confidence.  

The tax.summary file has the taxonimc level, the name of the taxonomic group, and the number of sequences in that group for each sample.  

We will also remove sequences that are classified to Chloroplast, Mitochondria, Unknown (not bacteria, archaea, or eukaryotes), Archaea, and Eukaryotes. 

```
remove.lineage(fasta=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
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

mothur "#classify.seqs(fasta=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds.tax)"

mothur "#remove.lineage(fasta=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)"

mothur "#summary.seqs(fasta=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)"
``` 

```
sbatch classify.sh
```

Classify.seqs will output the following files: 

```
Output File Names: 

oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy
oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary

```

The output file .taxonomy has name of sequence and the classification with % confidence in parentheses for each level.  

The tax.summary file has the taxonimc level, the name of the taxonomic group, and the number of sequences in that group for each sample.  

```
head oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy

head oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary
```

The remove.lineage command will remove sequences and provide a report in the output file. The fasta file is removing unique sequences, the count file is removing individual sequences.   

I altered the script to only pull out each lineage individually and ran for each lineage to get the individual stats. I then ran the final script with removing all lineages together to proceed with analysis.  

*Removing only Chloroplast*     
Removed 352 sequences from your fasta file.  
Removed 173441 sequences from your count file.   

*Removing only Mitochondria*     
Removed 1 sequences from your fasta file.  
Removed 6 sequences from your count file.     

*Removing only unknown domain*     
Removed 1598 sequences from your fasta file.  
Removed 43456 sequences from your count file.      

*Removing only Archaea*   
No contaminants to remove.   

*Removing only Eukaryotes*  
Removed 27 sequences from your fasta file.  
Removed 129 sequences from your count file.   

*Removing all lineages*   
Removed 1978 sequences from your fasta file.  
Removed 217032 sequences from your count file.     

We also can see the following summary:  

```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       125     38      0       2       1
2.5%-tile:      1       132     43      0       3       89092
25%-tile:       1       135     44      0       3       890911
Median:         1       140     46      0       3       1781821
75%-tile:       1       140     46      0       4       2672731
97.5%-tile:     1       140     47      0       6       3474550
Maximum:        5       140     54      0       8       3563641
Mean:   1       137     45      0       3
# of unique seqs:       59511
total # of seqs:        3563641
```

We now have 3,563,641 total sequences and 59,511 unique sequences.  

The script will then output the following files: 

```
Output File Names:

oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta
oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table
oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy
oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.accnos
oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table
oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta

```

Now we have the taxon removed that we do not want in our dataset.    

This is the end of the pipeline for curating our sequencing - we have corrected for pcr errors, removed chimeras, and removed sequences outside of taxon of interest. Now we can move onto OTU clustering!  

View a sample of the taxonomy summary dataset.   
 
```
head oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy
```

```
M00763_298_000000000-CBVP7_1_2108_14435_14093	Bacteria(94);Bacteria_unclassified(94);Bacteria_unclassified(94);Bacteria_unclassified(94);Bacteria_unclassified(94);Bacteria_unclassified(94);
M00763_298_000000000-CBVP7_1_2108_26527_14323	Bacteria(97);Bacteria_unclassified(97);Bacteria_unclassified(97);Bacteria_unclassified(97);Bacteria_unclassified(97);Bacteria_unclassified(97);
M00763_298_000000000-CBVP7_1_2103_18671_12758	Bacteria(95);Bacteria_unclassified(95);Bacteria_unclassified(95);Bacteria_unclassified(95);Bacteria_unclassified(95);Bacteria_unclassified(95);
M00763_298_000000000-CBVP7_1_2103_15136_12835	Bacteria(100);"Proteobacteria"(100);Gammaproteobacteria(100);Xanthomonadales(87);Xanthomonadaceae(87);Xanthomonadaceae_unclassified(87);
M00763_298_000000000-CBVP7_1_2109_28108_15591	Bacteria(95);Bacteria_unclassified(95);Bacteria_unclassified(95);Bacteria_unclassified(95);Bacteria_unclassified(95);Bacteria_unclassified(95);
M00763_298_000000000-CBVP7_1_2109_22001_15763	Bacteria(100);Bacteria_unclassified(100);Bacteria_unclassified(100);Bacteria_unclassified(100);Bacteria_unclassified(100);Bacteria_unclassified(100);
M00763_298_000000000-CBVP7_1_2108_9895_14769	Bacteria(100);"Proteobacteria"(96);Gammaproteobacteria(94);Gammaproteobacteria_unclassified(94);Gammaproteobacteria_unclassified(94);Gammaproteobacteria_unclassified(94);
M00763_298_000000000-CBVP7_1_2109_9770_16354	Bacteria(97);Bacteria_unclassified(97);Bacteria_unclassified(97);Bacteria_unclassified(97);Bacteria_unclassified(97);Bacteria_unclassified(97);
M00763_298_000000000-CBVP7_1_2109_25533_16525	Bacteria(100);Bacteria_unclassified(100);Bacteria_unclassified(100);Bacteria_unclassified(100);Bacteria_unclassified(100);Bacteria_unclassified(100);
M00763_298_000000000-CBVP7_1_2117_14902_17205	Bacteria(99);Bacteria_unclassified(99);Bacteria_unclassified(99);Bacteria_unclassified(99);Bacteria_unclassified(99);Bacteria_unclassified(99);
```

### <a name="Cluster"></a> **10. Cluster for OTUs**  

*In this analysis, we will cluster to OTU level (cutoff=0.03). For ASV clustering, we will move directly to the make.shared step, skipping the dist.seqs and cluster steps because mothur pre-clustering occurs at the ASV level.*  

First, we will calculate the pairwise distances between sequences.  

We will first use the `dist.seqs` command.  

```
dist.seqs(fasta=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta)
```

We will then `cluster` using the following command: 

```
cluster(column=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, cutoff=0.03)
```

This will run a line for each iteration of clustering. This is run until the Matthews correlation coefficient (MCC) value is maximized. A high MCC = high confidence in clustering. MCC is optimized by randomly aligning sequences to OTU's and calculating the correlation coefficient. Then sequences are moved between OTU's to see if the MCC is improved. This is repeated many times until the MCC is maximized. This method is fast and RAM efficient. AKA Opticlust.  

We would move directly to the make.shared step (below) if you want to use ASV since the distance and clustering steps do not need to be completed for ASVs.  

Next we will make a shared file. This .shared file has the label for the OTU, sample name, the number of OTU's and then the number of time each OTU appears in each sample. 

This file will be the basis of what we will do to measure richness of communities compared to each other.  

We want to keep the shared file and a consensus taxonomy file. 

```
make.shared(list=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
```

The classify.otu command will then output a concensus cons.taxonomy file that has the taxonomic information for each OTU. Use label=ASV to specify ASV in taxonomy names if starting from the make.shared step for ASV's.    

```
classify.otu(list=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy)
```

Then we can rename the files to something more useful.  

```
rename.file(taxonomy=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy, shared=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared)
```

Finally, view a count of the number of sequences in each sample.  

```
count.groups(shared=oyster.opti_mcc.shared)
``` 

Generate a script to run these commands to cluster for OTUs.    

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

mothur "#dist.seqs(fasta=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta)"

mothur "#cluster(column=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, cutoff=0.03)"

mothur "#make.shared(list=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)"

mothur "#classify.otu(list=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy)"

mothur "#rename.file(taxonomy=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy, shared=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared)"

mothur "#count.groups(shared=oyster.opti_mcc.shared)"
```

```
sbatch cluster.sh
```

Dist.seqs outputs the following files: 

```
Output File Names:
oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist

```

Cluster outputs the following files:  

```
Output File Names:
oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list
oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.steps
oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.sensspec

```

Make.shared outputs the following file: 

```
Output File Names:
oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared
```

Finally, classify.otu outputs the following file:  

```
Output File Names:
oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy
oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.tax.summary
```

The rename function at the end will rename our files to something more useful. We now have a "taxonomy" file and a "shared" file.  

The files we now care about are (OTUs): 

```
Current files saved by mothur:
list=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list
shared=oyster.opti_mcc.shared
taxonomy=oyster.taxonomy

constaxonomy=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy
count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table

```

If you wish to instead cluster to ASV level, use the following script:  

```
nano cluster_asv.sh
``` 

```
#!/bin/bash
#SBATCH --job-name="cluster_asv"
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

mothur "#make.shared(count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)"

mothur "#classify.otu(list=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.asv.list, count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=oyster.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, label=ASV)"

mothur "#rename.file(taxonomy=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.asv.ASV.cons.taxonomy, shared=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.asv.shared)"

mothur "#count.groups(shared=oyster.shared)"
```

```
sbatch cluster_asv.sh
``` 

make.shared outputs the following files:  

```
Output File Names:
oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.asv.list
oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.asv.shared
```

classify.otu outputs the following files:  

```
oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.asv.ASV.cons.taxonomy
oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.asv.ASV.cons.tax.summary
```  

The files we now care about are (if analyzed using the ASV script): 

Rename file outputs the following files:   

```
Current files saved by mothur:
list=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.asv.list
shared=oyster.shared
taxonomy=oyster.taxonomy

constaxonomy=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.asv.ASV.cons.taxonomy
count=oyster.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table

```

After clustering (either by ASV or OTU), the count of sequences in each file are:  

```
RS126 contains 26537.
RS127 contains 31369.
RS128 contains 27263.
RS129 contains 35508.
RS130 contains 33024.
RS131 contains 25644.
RS132 contains 29424.
RS133 contains 31335.
RS134 contains 32259.
RS135 contains 28863.
RS136 contains 23133.
RS137 contains 30410.
RS138 contains 30809.
RS139 contains 34207.
RS140 contains 31989.
RS141 contains 42299.
RS142 contains 28037.
RS143 contains 36516.
RS144 contains 37232.
RS145 contains 34310.
RS146 contains 25838.
RS147 contains 21026.
RS148 contains 31533.
RS149 contains 30380.
RS150 contains 27386.
RS151 contains 36060.
RS152 contains 19564.
RS153 contains 38932.
RS154 contains 28522.
RS155 contains 34888.
RS156 contains 30309.
RS157 contains 27190.
RS158 contains 27565.
RS159 contains 28632.
RS160 contains 38897.
RS161 contains 24631.
RS162 contains 20522.
RS163 contains 31028.
RS164 contains 37334.
RS165 contains 41914.
RS166 contains 27834.
RS167 contains 22934.
RS168 contains 21993.
RS169 contains 38077.
RS170 contains 30631.
RS171 contains 28954.
RS172 contains 26484.
RS173 contains 26953.
RS174 contains 26065.
RS175 contains 28324.
RS176 contains 33790.
RS177 contains 33733.
RS178 contains 29406.
RS179 contains 30546.
RS180 contains 41806.
RS181 contains 29137.
RS182 contains 27364.
RS183 contains 28041.
RS184 contains 26525.
RS185 contains 27731.
RS186 contains 24950.
RS187 contains 18869.
RS188 contains 27609.
RS189 contains 35365.
RS190 contains 37965.
RS191 contains 28830.
RS192 contains 47045.
RS193 contains 35197.
RS194 contains 30503.
RS195 contains 25633.
RS196 contains 29556.
RS197 contains 45666.
RS198 contains 34972.
RS199 contains 35204.
RS200 contains 26415.
RS201 contains 44030.
RS202 contains 22521.
RS203 contains 57790.
RS204 contains 36876.
RS205 contains 33721.
RS206 contains 45922.
RS207 contains 29156.
RS208 contains 31886.
RS209 contains 15011.
RS210 contains 29306.
RS211 contains 37072.
RS212 contains 42098.
RS213 contains 29690.
RS214 contains 34229.
RS215 contains 37034.
RS216 contains 31320.
RS217 contains 41697.
RS218 contains 41084.
RS219 contains 41901.
RS220 contains 37663.
RS221 contains 31146.
RS222 contains 31936.
RS223 contains 19927.
RS224 contains 24245.
RS225 contains 9082.
RS226 contains 38262.
RS227 contains 38422.
RS228 contains 38152.
RS229 contains 34276.
RS230 contains 35053.
RS231 contains 38347.
RS232 contains 32695.
RS233 contains 41450.
RS234 contains 30628.
RS235 contains 36436.
RS246 contains 38049.
RS247 contains 23102.

Size of smallest group: 9082.

Total seqs: 3563641.
```

*Now we have classified taxonomic data and we are ready to move onto subsampling and calculating statistics.*    

### <a name="Subsample"></a> **11. Subsampling for Sequencing Depth**   

**The remaining commands can all be run in interactive mode. In Andromeda, use the following commands. These commands could also be combined into a script. I am using interactive mode for now.**  

*All analyses from this point forward are for ASV level analysis. If using OTU's, revise file names in these scripts accordingly.*   

```
interactive 
module load Mothur/1.46.1-foss-2020b
mothur
```

This will open the interactive mothur terminal.  In order to use any bash commands like `nano`, you have to first `quit` in mothur. Once you want to return to mothur, use the `mothur` command.  

#### Subsampling for sequencing depth 

As we can see in rarefaction curves, we know that OTU discovery is a function of sequencing depth. Therefore, we need to account for sequencing depth in our analysis.  

To account for uneven sampling depth, we will need to rarify our samples. We can do this with the `sub.sample` command that will automatically use the smallest group size or we can manually set a cut off. If we set a threshold higher than the number of samples, those samples are removed from the group.  

We will base this number on the smallest sequence number from our target samples - the gut samples.    

The samples we want are:  
```
RS126
RS127
RS128
RS129
RS130
RS131
RS132
RS133
RS134
RS135
RS136
RS137
RS138
RS139
RS140
RS141
RS142
RS143
RS144
RS145
RS146
RS147
RS148
RS149
RS150
RS151
RS152
RS153
RS154
RS155
RS156
RS157
RS158
RS159
RS160
RS161

RS247 (negative control)  
RS246 (mock community)  
```

From this sample list, the smallest sequence per sample was RS152 with 19,564 sequences.  Therefore, we will subset to this value. If other samples are dropped from the list because they are lower than this value that is ok, we are targetting the gut samples because we will only be analyzing this data set.  


*NOTE FOR AH: ADD IN MOCK COMMUNITY ANALYSIS*  

```
sub.sample(shared=oyster.shared, size=19564) 
```

```
Output files: 
oyster.ASV.subsample.shared
```

It may be helpful to run this with multiple iterations. If you do not include a size=# argument, then mothur will automatically set to the lowest sequence. If I set this to 19564, for example the output shows samples are removed because there are <19564 sequences: 

```
RS187 contains 18869. Eliminating.
RS209 contains 15011. Eliminating.
RS225 contains 9082. Eliminating.
Sampling 19564 from each group.
```  

THere are three samples below our subsampling threshold. But again these are not in our target dataset, so that is OK.  

#### Subsampling the sample set  

We are going to generate a rarefaction curve and subsample to 19564 sequences. Calculating the observed number of OTUs using the Sobs metric (observed number of taxa) with a freq=100 to output the data for every 100 sequences samples - if you did every single sample that would be way too much data to output.  

```
rarefaction.single(shared=oyster.shared, calc=sobs, freq=100)
```

```
Output File Names: 
oyster.groups.rarefaction
```

We will want to keep this rarefaction file for our final output.  

You can look at this file with `nano oyster.groups.rarefaction`. It contains the number of ASVs detected for every 100 sequences for each sample. Series will be trucated after the number of sequences available in the sample.  

Next, rarefy the data to the minimum number of sequences, then the next command will rarefy and pull out that many sequences from each sample. If you want to sample a specific number, use subsample=#. You would want to do this if you have a super small sequence number in one group and you dont want to truncate the data for all groups by that much. If you leave subsample=T, mothur will calculate to the minimum size by default (we want to manually set the level in this case). Subsamples 1000 times and calculates the metrics we want with the number of otus. Sampling is random.   

All of the available metric calculations (calc) are [found here](https://mothur.org/wiki/calculators/).  

I am going to use subsample=19564 to remove the low sequence count stages.  

```
summary.single(shared=oyster.shared, calc=nseqs-sobs-chao-shannon-invsimpson, subsample=19564)
```

In a new Andromeda window (outside mothur).  
```
less oyster.groups.ave-std.summary
```

Open this file to view statistics for each sample. The method column has ave (average of each of the metrics across the 1000 iteractions) or std (this is the std deviation of the average calculation).    

### <a name="Statistics"></a> **12. Calculate Ecological Statistics**  

Alpha diversity = diversity within a sample 
Beta diversity = diversity between samples

We will calculate beta diversity metrics which can then be used for ecological analyses.  

#### Calculate Beta Diversity  

We will use the `dist.shared` function to calculate distances between samples. Here we are going to use Bray-Curtis distances. Many different metrics can be used (e.g., Theta YC). All of the available metric calculations (calc) are [found here](https://mothur.org/wiki/calculators/).  

We can use this function with `subsample=19564`. We will use the bray-curtis distances for this calculation. This step is important for generating files that we will use later in R. 

If we run a command with subsample=19564, then the following files are output: 

```
dist.shared(shared=oyster.shared, calc=braycurtis, subsample=19564)
```

```
Output File Names: 
oyster.braycurtis.lt.ave.dist
oyster.braycurtis.lt.std.dist
``` 

We have average and std files because the function is run over 1000 iterations generating ave and std.         

*After this point, analyses can be run in R. See the bottom of this document for transferring files to your computer.*    

### <a name="Output"></a> **13. Output data for analysis in R**  

We can now move some files into R for further analysis. These are the primary files we will be using - the bray curtis distance matrix (distances between samples), the taxonomy file (taxonomy of each otu), and the shared file (otu abundance in each sample). From these we can run PCoA, NMDS, and statistical testing in R.    

Bray-Curtis matrices:     

```
oyster.braycurtis.lt.ave.dist
```

Taxonomy file:  

```
oyster.taxonomy
```

Shared file:  

```
oyster.ASV.subsample.shared
```

Rarefaction file:  
```
oyster.groups.rarefaction
```

Outside Andromeda do the following for each file:  

```
scp ashuffmyer@ssh3.hac.uri.edu:/data/putnamlab/ashuffmyer/oyster_16S/v6/oyster.braycurtis.lt.ave.dist ~/MyProjects/Cvir_Nut_Int/output/16S_allv6/mothur

scp ashuffmyer@ssh3.hac.uri.edu:/data/putnamlab/ashuffmyer/oyster_16S/v6/oyster.taxonomy ~/MyProjects/Cvir_Nut_Int/output/16S_allv6/mothur

scp ashuffmyer@ssh3.hac.uri.edu:/data/putnamlab/ashuffmyer/oyster_16S/v6/oyster.ASV.subsample.shared ~/MyProjects/Cvir_Nut_Int/output/16S_allv6/mothur

scp ashuffmyer@ssh3.hac.uri.edu:/data/putnamlab/ashuffmyer/oyster_16S/v6/oyster.groups.rarefaction ~/MyProjects/Cvir_Nut_Int/output/16S_allv6/mothur

scp ashuffmyer@ssh3.hac.uri.edu:/data/putnamlab/ashuffmyer/oyster_16S/v6/*.txt ~/MyProjects/Cvir_Nut_Int/output/16S_allv6/mothur

scp ashuffmyer@ssh3.hac.uri.edu:/data/putnamlab/ashuffmyer/oyster_16S/v6/*.csv ~/MyProjects/Cvir_Nut_Int/output/16S_allv6/mothur
```

After copying to the computer, I resave these files at tab delimited (.txt) before loading into R.  




 




