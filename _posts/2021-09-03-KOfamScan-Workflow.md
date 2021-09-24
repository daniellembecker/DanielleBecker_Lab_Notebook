---
layout: post
title: P. verrucosa KofamScan pipeline
Author: Danielle Becker-Polinski
Last Updated: 2021/09/03
tags: [ Protocol, KEGG, KO terms, KOfamscan ]
---

## Overview

Follow step 3: map KEGG terms to a genome in @echille's lab notebook [2020-10-08-M-capitata-functional-annotation-pipeline.md ](https://github.com/echille/E.-Chille-Open-Lab-Notebook/blob/master/_posts/2020-10-08-M-capitata-functional-annotation-pipeline.md) post for general steps and explanations on the beginning of this process.

Also, see @echille GitHub for further [KEGG ontology steps in R](https://github.com/echille/Mcapitata_OA_Developmental_Gene_Expression_Timeseries/tree/main/5-Planula-GO-Enrichment-Analysis/a-Kegg-ontology).

### Step 3: Map Kegg terms to genome  
KofamScan is the command-line version of the popular KofamKOALA web-based tool, used to map Kegg terms (containing pathway information) to a genes. KofamScan and KofamKoala work by using HMMER/HMMSEARCH to search against KOfam (a customized HMM database of KEGG Orthologs (KOs). Mappings are considered robust because each Kegg term has an individual pre-defined threshold that a score has to exceed in order to map to a gene. While all mappings are outputted, high scoring (significant) assignments are highlighted with an asterisk.

The commands that I used are below. In order to run KofamScan, you will need a fasta file of predicted protein sequences (preferably the same one used to run InterProScan).

General Protocol:

#### i) Download and inflate the Kofam database.

To get the most up-to-date Kofam database, download it just before running KofamScan. You will also need to download the profiles associated with the Kofam database containing threshold information.

```
curl -O ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz | gunzip > ko_list #download and unzip KO database
curl -O ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz | tar xf > profiles #download and inflate profiles
```

#### ii) Run KofamScan

**exec_annotation: executes Kofamscan**

To execute Kofamscan, use the exec_annotation script provided with program download. To run, use: ```exec_annotation {options} path_to_query_file```.

*Options (See KofamScan [repo](https://github.com/takaram/kofam_scan) for more options):*

- **-o** - set name of output file
- **-k** - path to ko database (downloaded above)
- **-p** - path to profile database (downloaded above)
- **-E** - set minimum expect value for significant mappings
- **-f** - output format
- **--report-unannotated** - returns names of sequences with no mapped KO terms
- **-pa** - enables Kegg term mapping

```
/opt/software/kofam_scan/1.3.0-foss-2019b/exec_annotation -o Pver_KO_annot.tsv -k ./ko_list -p ./profiles/eukaryote.hal -E 0.00001 -f detail-tsv --report-unannotated /data/putnamlab/REFS/Pverr/Pver_proteins_names_v1.0.faa
```

Bluewaves Script:

```
#!/bin/bash
#SBATCH --job-name="KofamScan"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=100GB
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/KEGG/

echo "Loading modules" $(date)
module load kofam_scan/1.3.0-foss-2019b
module load libyaml/0.1.5
module unload HMMER/3.3.1-foss-2019b
module load HMMER/3.3.2-gompi-2019b
module list

#echo "Starting analysis... downloading KO database" $(date)
#wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz #download KO database
#wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
#gunzip ko_list.gz
#tar xf profiles.tar.gz

echo "Beginning mapping" $(date)
/opt/software/kofam_scan/1.3.0-foss-2019b/exec_annotation -o Pver_KO_annot.tsv -k ./ko_list -p ./profiles/eukaryote.hal -E 0.00001 -f detail-$

echo "Analysis complete!" $(date)
```

### Step 4: Compilation of the output of different methods

Done in RStudio. See RMarkdown [script](https://github.com/hputnam/Becker_E5/blob/master/RAnalysis/Scripts/RNA-seq/KEGG_Pathway_Analysis.Rmd).

Output files can be found [here](https://github.com/hputnam/Becker_E5/tree/master/RAnalysis/Genome/KEGG-ontology)

Issues: [GitHub Issue #21](https://github.com/Putnam-Lab/Lab_Management/issues/21)

I ran into an issue when downloading and using this part of @echille code to get the up-to-date Kofam database:

`curl -O ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz | tar xf > profiles #download and inflate profiles`

It seemed to not download/extract on of the HMM files for some reason which I realized after running the KofamScan script, adapted from @echille original found [here](https://github.com/echille/Mcapitata_OA_Developmental_Gene_Expression_Timeseries/blob/main/5-Planula-GO-Enrichment-Analysis/a-Kegg-ontology/Mcap_KofamScan.sh).

In the slurm output, this error was written out:

```
hmmsearch was not run successfully

Error: File existence/permissions problem in trying to open HMM file
/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/KEGG/profiles/K00637.hmm.
HMM file /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/KEGG/profiles/K00637.hmm not found (nor an .h3m binary of it)
```

Some troubleshooting was done to check this error, you can search the profiles folder to see if a .hmm file was extracted correctly:

#how to check for file extraction in profiles folder:

`$ ls profiles/K00637*`

ls: cannot access 'profiles/K00637*': No such file or directory

#how to extract .hmm file from profiles.tar.gz if they do not appear initially

`$ tar tf profiles.tar.gz | grep K00637`

profiles/K00637.hmm


We were able to figure out that another way to download and extract all of the HMM files was to use wget:

`wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz`

Make sure to unzip the profiles.tar.gz after:

`tar -xf profiles.tar.gz > profiles`

After I solved this issue, I got another error in the slurm output that read:

`Error: Unknown KO: K00960`

This error means that for some reason when extracting/downloading the Kofam database, some of the values (K00960) and a few others that had not been updated were still being found in the tmp > tabular output folder that is created during the download/extraction.

I used the code below to see the older versions of the values that are no longer used and then had to delete them before running the job again.

```
ls -ltr tmp/tabular/  | grep -v Aug.19

total 232797
```

```
drwxr-xr-x 3 danielle_becker putnamlab   4096 Jul 16 12:03 ..
-rw-r--r-- 1 danielle_becker putnamlab   2238 Jul 19 13:34 K00960
-rw-r--r-- 1 danielle_becker putnamlab  26732 Jul 19 13:36 K02471
-rw-r--r-- 1 danielle_becker putnamlab  18298 Jul 19 13:52 K09311
-rw-r--r-- 1 danielle_becker putnamlab  21632 Jul 19 13:53 K09376
-rw-r--r-- 1 danielle_becker putnamlab   1654 Jul 19 13:53 K09394
-rw-r--r-- 1 danielle_becker putnamlab   3422 Jul 19 14:12 K18122
```

These values did not appear to be in the new profiles.tar.gz file, so they had to have come from a previous version.

I removed them and re-ran the job.


If you are still getting an error, make sure to re-download the ko_list.gz as well before re-running

`curl -O ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz | gunzip > ko_list #download and unzip KO database`

Sometimes the gunzip on the end does not work, make sure to check that the ko_list has contents in it. If not, use:

`gunzip ko_list.gz`  in a separate line to properly unzip the gz file.

After all of these steps, I re-ran the script and it worked! You should have a  Pver_KO_annot.tsv file that looks like this:

<img width="832" alt="Screen Shot 2021-08-23 at 2 29 15 PM" src="https://user-images.githubusercontent.com/40968353/130498319-0e27daa0-0446-465a-a655-94535c78e27f.png">
