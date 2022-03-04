---
layout: post
title: Putnam Lab Oyster 16S Analysis in Mothur
date: '2022-03-03'
categories: PutnamLab_Oysters
tags: 16S Molecular
---
This post details 16S analysis of the Putnam Lab oyster project using the Mothur pipeline based on [A. Huffmyer's pipeline here](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/16S-Analysis-in-Mothr-Part-1/).   

# **16S analysis workflow in mothur**  

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
cd oyster_16S
```  

Give permissions for all to read, write, and execute.   

`chmod u=rwx,g=rwx,o=r -R oyster_16S`  

Next, create symlinks to the files that we will need for analysis.  

```
ln -s /data/putnamlab/shared/PointJudithData_Rebecca/amplicons16s/gut_v4v5/*.txt /data/putnamlab/ashuffmyer/oyster_16S/

ln -s /data/putnamlab/shared/PointJudithData_Rebecca/amplicons16s/gut_v4v5/*.csv /data/putnamlab/ashuffmyer/oyster_16S/

ln -s /data/putnamlab/shared/PointJudithData_Rebecca/amplicons16s/gut_v4v5/00_RAW_gz/*.fastq.gz /data/putnamlab/ashuffmyer/oyster_16S/
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

We will also run the `summary.seqs()` command with the `make.file()` and `make.contig()` steps. This generates summary information about the sequences from the files we made above.  

Within the `oyster_16S` directory, create a script. Everything will live within this directory - I did not make any subfolders.  

```
nano contigs.sh
```

*CURRENTLY RE RUNNING WITH TRIM OVERLAP = F*  

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

module load Mothur/1.46.1-foss-2020b

mothur

mothur "#make.file(inputdir=., type=gz, prefix=oyster)"

mothur "#make.contigs(inputdir=., outputdir=., file=oyster.files, trimoverlap=F)"

mothur "#summary.seqs(fasta=oyster.trim.contigs.fasta)"

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

Count the number of sequences that were removed and the number that were kept by counting sequences in each fasta file.  

```
grep -c "^>" oyster.trim.contigs.fasta
grep -c "^>" oyster.scrap.contigs.fasta
```
4,800,906 sequences were kept and 8,128 were removed.  

Next, check whether we need to remove primers.  

The primers we are looking for are: 

Quince et al. 2011

518F: 5’ TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCAGCAGCYGCGGTAAN
926R: 5’ GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGCCGTCAATTCNTTTRAGT 
      5’ GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGCCGTCAATTTCTTTGAGT
      5’ GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGCCGTCTATTCCTTTGANT

When I search for the primer sequences in this file they do not show up. Grep to count for any total number of primer sequences.  

```
grep -c "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCAGCAGCYGCGGTAAN" oyster.trim.contigs.fasta
grep -c "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGCCGTCAATTCNTTTRAGT" oyster.trim.contigs.fasta
grep -c "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGCCGTCAATTTCTTTGAGT" oyster.trim.contigs.fasta
grep -c "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGCCGTCTATTCCTTTGANT" oyster.trim.contigs.fasta
```
The primers show up 0 times in the sequences.  

You will get an output table like the one below and it will output a file with this summary information.  

```
less contigs_output_script

```

Scroll to the bottom to find the summary information.  


``` 
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       1       1       0       1       1
2.5%-tile:      1       4       4       0       2       120023
25%-tile:       1       9       9       1       2       1200227
Median:         1       10      10      1       2       2400454
75%-tile:       1       10      10      3       2       3600680
97.5%-tile:     1       94      94      8       5       4680884
Maximum:        1       251     251     59      35      4800906
Mean:   1       15      15      2       2
# of Seqs:      4800906
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

mothur "#screen.seqs(inputdir=., outputdir=., fasta=oyster.trim.contigs.fasta, group=oyster.contigs.groups, maxambig=0, maxlength=XXX, minlength=XXX)"

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

X sequences were kept and X were removed. 

This removed ~X% of sequences due to length and ambiguous bases.  

You can also view the `bad.accnos` file to see why sequences were removed.  

```
head -1000 oyster.trim.contigs.bad.accnos
``` 

The summary output as viewed by `less screen_output_script` file now reads: 

```
INSERT HERE

```

We now see that we have removed all sequences with ambigous calls and the max sequence length is X with min X. Next we will align to the reference database.  

### <a name="Unique"></a> **5. Determining and counting unique sequences**  

Next, determine the number of unique sequences. The code will look like this:  

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

In this run, there were XX sequences and XX were unique = ~XX%.  

The output table from `summary.seqs` looks like this and shows the number of unique and the total: 

```
INSERT HERE

```

Now we can align just the unique sequences, which will be much faster than aligning the full data set and is an indicator of how polished and clean the data are.  

*From this, we have our unique sequences identified and can proceed with further cleaning and polishing of the data. Next we will look at alignment, error rate, chimeras, classification and further analysis.*  

https://mothur.org/blog/2016/Customization-for-your-region/  

