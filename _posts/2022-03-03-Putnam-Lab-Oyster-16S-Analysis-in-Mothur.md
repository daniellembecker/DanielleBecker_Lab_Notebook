---
layout: post
title: Putnam Lab Oyster 16S Analysis in Mothur
date: '2022-03-03'
categories: PutnamLab_Oysters
tags: 16S Molecular
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
6. [Create a V6 database](#V6)   
7. [Aligning](#Align)    
8. [Preclustering](#Precluster)    
9. [Identify chimeras](#Chimera)  
10. [Classify sequences](#Classify)     
11. [Cluster OTU's](#Cluster)    
12. [Subsampling](#Subsample)  
13. [Calculate ecological statistics](#Statistics)  
14. [Output data for R analysis](#Output)   

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

### <a name="V6"></a> **6. Generate a reference for the V6 region** 

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


  




