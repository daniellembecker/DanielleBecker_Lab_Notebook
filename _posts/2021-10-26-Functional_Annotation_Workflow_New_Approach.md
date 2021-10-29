---
layout: post
title: Functional Annotation Pipeline
Author: Danielle Becker-Polinski
Last Updated: 2021/10/26
tags: [ Protocol, annotation, RNASeq, GO, KEGG ]
---


## Overview

Testing out new approaches for functional annotation of  *Pocillopora verrucosa*. Previous analyses (Erin and Jill) used the program [DIAMOND BLAST](http://www.diamondsearch.org/index.php) against the NCBI nr database. Like regular BLAST, DIAMOND is a sequence aligner for nucleotide and protein sequences; unlike BLAST, it is optimized for a higher performance capability of large datasets at 100x-20,000x speed of BLAST.  The output .xml file was then funneled to [BLAST2GO](https://www.blast2go.com) (B2G) which is a bioinformatics tools for functional annotation and analysis of gene or protein sequences. It was originally developed to provide a user-friendly interface for GO annotation and now hosts many different functional annotation tools. The B2G software makes it possible to generate annotation without requiring writing any code. While B2G can do an extensive array of functions, this analysis primarily utilizes the GO mapping and annotation functions.

Through a dense literature search, multiple databases should be used for BLAST to expand the hits possible for your sequences (Buitrago-López et al 2020; Baumgarten et al. 2015; Cunning et al. 2018). The main concensus from this literature search is that the protein sequences should be searched against the SwissProt, TrEMBL, NCBI nr databases using BLASTp (Basic Local Alignment Search Tool, e-value cut-off = 1e-05) and retaining annotations from databases in this order. Then, BLAST2GO should be used to provide GO annotations, and KEGG, Pfam, InterProScan, should be searched to further annotated gene sets.


### Step 1: Obtain sequences of interest.

In order to conduct functional annotation steps, protein and transcript sequences are needed. There are two main hubs where coral genomic information is stored: [Reef Genomics](http://reefgenomics.org) and [NCBI](https://www.ncbi.nlm.nih.gov). Other researchers store their genomic infomation on their own personal webpages. Genomic information must be downloaded from one of these databases in order to proceed.

#### i) Identify species to work with.

For this project, the coral species of interest is *Pocillopora verrucosa*.

#### ii) Download genomic files for species of interest.

##### [*Pocillopora verrucosa* Full Transcripts ](http://pver.reefgenomics.org/download/)

##### [*Pocillopora verrucosa* Gene models (protein)](http://pver.reefgenomics.org/download/)

Ready to start annotating!

### Step 2: Identify homologous sequences

Homology refers to the similarity of structure or genes in different taxa due to shared ancestry. {example}

Sequence homology is the homology between DNA, RNA, and protein sequences in terms of shared ancestry. Sequence homology is usually inferred by the similarity of nucleotide or amino acid sequences. Strong sequence similarity (or percent homology) provides evidence that two or more sequences are related through shared ancestry.

Several software programs have the ability to compare sequences to assess homology, the most popular one being NCBI's [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) (Basic Local Alignment Search Tool - what a great name). BLAST compares nucleotide or protein sequences of interest (called a query) to their sequence databases to find regions of similarity. If a nucleotide or protein sequence of interest significantly matches a sequence/sequences in the databases, BLAST will tag the sequence of interest with the information about the known sequence(s). For example, if a new transcript sequence is identified in a mouse, BLAST could be used to see if any other animals carry a similar sequence, and if they do, what the biological functions of the sequence are.

Other software programs include [SWISS-PROT](https://iop.vast.ac.vn/theor/conferences/smp/1st/kaminuma/SWISSPROT/access.html) which is a curated protein sequence database that provides a high level of annotation (such as the description of the function of a protein, its domain structure, post-translational modifications, variants, etc), a minimal level of redundancy and a high level of integration with other databases. Recent developments of the database include: an increase in the number and scope of model organisms; cross-references to seven additional databases; a variety of new documentation files; the creation of [TREMBL](http://www3.cmbi.umcn.nl/wiki/index.php/TrEMBL), an unannotated supplement to SWISS-PROT.

### Step 3: Download your databases

#### i) On the Andromeda server, download Swiss-Prot database from [UniProt/Swiss-Prot](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz).

I created a script on Andromeda to download the updated Swiss-Prot database, unzip it, and make it into a blast database. I followed the general instructions on how to create code to download database [here](https://www.biostars.org/p/354449/).  Current databases found [here](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/).

```
pwd /data/putnamlab/shared/sbatch_executables/download_swissprot_database.sh
```

Full script:

```
#!/bin/bash
#SBATCH --job-name="ref"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu # CHANGE EMAIL
#SBATCH -D /data/putnamlab/shared/databases/swiss_db

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b

cd databases/swiss_db

echo "Making swissprot database" $date
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
makeblastdb -in uniprot_sprot.fasta -parse_seqids -dbtype prot -out swissprot_20211022
echo "STOP" $(date)

```

**It is normal to get multiple files per blast database. That is how makeblastdb is supposed to work. Just make sure files for a database stay together in the same directory and you use the "basename" for the database (a suggestion: name your database some thing other than your input file name) when you run your searches.**

**Options for makeblastdb which is an application that produces BLAST databases from FASTA files, other options [here](http://nebc.nerc.ac.uk/nebc_website_frozen/nebc.nerc.ac.uk/bioinformatics/documentation/blast+/user_manual.pdf)**

- ```makeblastdb``` - application produces BLAST databases from FASTA files
- ```in``` - database that will be used to search against
- ```parse_sequids``` - parse the seq-id(s) in the FASTA input provided
- ```dbtype``` - molecule type of target db ("nucl" or "prot")
- ```out``` - base output name for database file


#### ii) On the Andromeda server, download Trembl database from [UniProt/Trembl](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz).

I created a script on Andromeda to download the updated Trembl database, unzip it, and make it into a blast database. I followed the general instructions on how to create code to download database [here](https://www.biostars.org/p/354449/).  Current databases found [here](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/).

```
pwd /data/putnamlab/shared/sbatch_executables/download_trembl_database.sh
```

Full script:

```
#!/bin/bash
#SBATCH --job-name="ref"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu # CHANGE EMAIL
#SBATCH -D /data/putnamlab/shared/databases/trembl_db

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b

cd databases/trembl_db

echo "Making trembl database" $date
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
gunzip uniprot_trembl.fasta.gz

makeblastdb -in uniprot_trembl.fasta -parse_seqids -dbtype prot -out trembl_20211022
echo "STOP" $(date)

```

**It is normal to get multiple files per blast database. That is how makeblastdb is supposed to work. Just make sure files for a database stay together in the same directory and you use the "basename" for the database (a suggestion: name your database some thing other than your input file name) when you run your searches.**

**Options for makeblastdb which is an application that produces BLAST databases from FASTA files, other options [here](http://nebc.nerc.ac.uk/nebc_website_frozen/nebc.nerc.ac.uk/bioinformatics/documentation/blast+/user_manual.pdf)**

- ```makeblastdb``` - application produces BLAST databases from FASTA files
- ```in``` - database that will be used to search against
- ```parse_sequids``` - parse the seq-id(s) in the FASTA input provided
- ```dbtype``` - molecule type of target db ("nucl" or "prot")
- ```out``` - base output name for database file


#### iii) On the Andromeda server, download updated [NCBI](ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz) nr.gz database and make it into a binary format for DIAMOND BLAST.

The program [DIAMOND BLAST](http://www.diamondsearch.org/index.php) will be used. Like regular BLAST, DIAMOND is a sequence aligner for nucleotide and protein sequences; unlike BLAST, it is optimized for a higher performance capability of large datasets at 100x-20,000x speed of BLAST.

**If you are in the Putnam Lab and on Andromeda**  

**I used this step for this functional annotation because using the below commands in bash took way too long**

This script, created by Erin Chille on August 6, 2020, downloads the most recent nr database in FASTA format from [NCBI](ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz) and uses it to make a Diamond-formatted nr database.  This step was updated by Danielle Becker-Polinski on September 24th, 2021 because the scripts were not including the full CPUs to download and a couple other formatting errors. Go to the *sbatch_executables* subdirectory in the Putnam Lab *shared* folder and run the scripts, ```make_diamond_nr_db.sh```  and  ```make_diamond_nr_db.sh``` in this order:

```
$ sbatch download_nr_database.sh
Submitted batch job NNN
$ sbatch -d afterok:NNN make_diamond_nr_db.sh
```

**If you are not in the Putnam Lab**  
Download the nr database from [NCBI](ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz). Then, use Diamond's ```makedb``` command to format the database in a Diamond-friendly format. You can also use the command ```dbinfo``` to find version information for the database.

The nr (non-redundant) database is a collection of non-identical protein sequences compiled by NCBI. It is updated on a daily basis.

```
# On Andromeda

# Download db from NCBI
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz

# Load Diamond module
module load DIAMOND/2.0.0-GCC-8.3.0

diamond makedb --in nr.gz -d nr
diamond dbinfo -d nr.dmnd
```

**Options**

- ```makedb``` - create a Diamond binary database file
- ```in``` - database that will be used to search against
- ```d``` - base output name for database file
- ```dbinfo``` - print information about database file

The file nr.dmnd is now Diamond-readable and can be used as a reference.

### Step 3: Align query protein sequences against databases

Now that the reference databases have been properly generated, the sequences of interest can be aligned against them.

#### i) BLAST the protein sequences against Swiss-Prot

```
pwd /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/Functional_Annotation/Swissprot/swissprot_blast.sh

```

Full script:

```
#!/bin/bash
#SBATCH --job-name="swissprot-blastp-protein"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="swissprot_blastp_out_error"
#SBATCH --output="swissprot_blastp_out"
#SBATCH --exclusive

echo "START" $(date)
module load BLAST+/2.11.0-gompi-2020b #load blast module

echo "Blast against swissprot database" $(date)

blastp -max_target_seqs 5 -num_threads 20 -db /data/putnamlab/shared/databases/swiss_db/swissprot_20211022 -query /data/putnamlab/REFS/Pverr/Pver_proteins_names_v1.0.faa -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -out PverGeneModels_vs_sprot_1e-5_max5.out

echo "STOP" $(date)

```

# Get the best hit for each Gene Model (protein) Swiss-Prot

```
#Sort by 1. query name, 2. bitscore, 3. evalue, 4. protein identity, and extract the best line for each query (bitscore more important than evalue, evalue more important than nucleotide identity).

cat PverGeneModels_vs_sprot_1e-5_max5.out | sort -k1,1 -k2,2 -k3,3r -k4,4r -k11,11 | awk '!seen[$1]++' > PverGeneModels_vs_sprot_1e-5_besthit.out

wc -l PverGeneModels_vs_sprot_1e-5_besthit.out #19,540

```

# Select the gene model proteins without hits in the swiss prot

```
#first use awk to print a list of all the Gene Model names from besthits.out

awk '{print $1}' PverGeneModels_vs_sprot_1e-5_besthit.out > list_of_Pvergenemodelproteins_sprot.txt

#then exclude these Gene Model names from your original fasta/.faa/protein file

-exclude Pver_proteins_names_v1.0.faa list_of_Pvergenemodelproteins_sprot.txt Pver_proteins_names_v1.0.faa.prot4trembl

#exclude command not found in bash, emailed Kevin Bryan for follow up
```


#### ii) BLAST the remaining protein sequences against Trembl

```
pwd /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/Functional_Annotation/Trembl/trembl_blastp.sh

```


#### iii) BLAST the remaining protein sequences against nr

```
# On Andromeda

# Load Diamond module
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

#Run sequence alignment against the nr database
diamond blastx -d /data/putnamlab/shared/databases/nr.dmnd -q /data/putnamlab/REFS/Pverr/Pver_transcriptome_v1.0.fasta -o Pver_annot -f 100 -b20 --more-sensitive -e 0.00001 -k1 --unal 1 --threads $SLURM_CPUS_ON_NODE --tmpdir /tmp/

```

Blastx: align translated DNA query sequences against protein reference database

**Options**

- ```d``` - path to nr database file
- ```q``` - path to query fasta file
- ```o``` - base output name
- ```f``` - output format
    - 100 = Diamond format
- ```b``` - Block size in billions of sequence letters to be processed at a time. Larger block sizes increase the use of memory and temporary disk space, but also improve performance. Set at 20. 20 is the highest recommended value.
- ```more-sensitive```
- ```e``` - maximum expected value to report an alignment
    - 1e-05 is typical cutoff for sequence alignments
- ```k``` maximum top sequences
    - Set at 1 to only report top sequence for each gene

#### iii) Generate readable output files

```
# On Andromeda

#Converting format to XML format for BLAST2GO
diamond view -a Pver.annot.20210924.daa -o Pver.annot.20210924.xml -f 5
diamond view -a Pver.annot.20210924.daa -o Pver.annot.20210924.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

```

View: generate formatted output files

**Options**

- ```a``` - path to input DAA file
- ```o``` - base output name
- ```f``` - output format
    - 5 = XML output; 6 = TAB output

The output files (XML and TAB) will both be used downstream in this workflow.

#### iv) View output file

This is an example of part of the DIAMOND output .tab file:

qseqid | sseqid | pident | length | mismatch | gapopen | qstart | qend | sstart | send | evalue | bitscore | qlen | slen
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
Pver_evm.model.Segkk0_pilon.12 | XP_015753513.1 | 86.7 | 1103 | 15 | 1 | 1 | 2913 | 23 | 1125 | 0.0e+00 | 1807.7 | 2916 | 1125

**Column names**

- qseqid - query seq-id
- sseqid - subject seq-id (from BLAST databases)
- pident - % identical matches
- length - alignment length
- mismatch - number of mismatches
- gapopen - number of gap openings
- qstart - start of alignment in query
- qend - end of alignment in query
- sstart - start of alignment in subject
- send - end of alignment in subject
- evalue - number of expected hits of similar quality that could be found by chance (similar to a pvalue). The smaller the evalue, the better the match!
- bitscore - required size of a sequence database in which the current match could be found just by chance. The higher the bitscore, the better the sequence similarity!
- qlen - query sequence length
- slen - subject sequence length

#### v) Secure-copy output files to local computer

```
# From a new terminal window (ie not Andromeda or remote server)

scp danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/BLAST-GO-KO/Diamond/Pver_annot.xml /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/Functional_Annotation/Diamond

scp danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/BLAST-GO-KO/Diamond/Pver_annot.tab /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/Functional_Annotation/Diamond
```

DIAMOND BLAST results can now be used in further analyses.

**Full Andromeda Script:**
**This step will take over 20 hours, you can run this and the InterProScan script at the same time on Andromeda**

```
Pver_annot_diamond.sh:

#!/bin/bash
#SBATCH --job-name="diamond-blastx"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="diamond_blastx_out_error"
#SBATCH --output="diamond_blastx_out"
#SBATCH --exclusive

echo "START" $(date)
module load DIAMOND/2.0.0-GCC-8.3.0 #Load DIAMOND

echo "Updating Pver annotation" $(date)
diamond blastx -d /data/putnamlab/shared/databases/nr.dmnd -q /data/putnamlab/REFS/Pverr/Pver_transcriptome_v1.0.fasta -o Pver_annot -f 100 -b20 --more-sensitive -e 0.00001 -k1 --unal 1 --threads $SLURM_CPUS_ON_NODE --tmpdir /tmp/


echo "Search complete... converting format to XML and tab"

diamond view -a Pver_annot.daa -o Pver_annot.xml -f 5
diamond view -a Pver_annot.daa -o Pver_annot.tab -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

echo "STOP" $(date)
```
```
Submitted batch job 1931742
```

### Step 3: Assign gene ontology terms to sequences

After DIAMOND BLAST is completed, analysis can move to assigning gene ontology (GO) terms to sequences.

The [Gene Ontology](http://geneontology.org) is an extensive consortium that aims to understand gene function and provide functional annotation for organisms across the tree of life. It also maintains a controlled vocabulary of gene and gene attributes across species.

The Gene Ontology has a system to classify genes into terms. The terms are grouped into 3 categories:

1. Molecular funciton - molecular-level activities performed by gene products (ie RNA or proteins). Describes the activities rather than the entities (molecules or protein complexes)
    - Examples: Toll receptor binding, transcription regulator activity
2. Cellular component - cellular anatomy and/or locations where gene products perform function
    - Examples: mitochondrion, ribosome
3. Biological process - larger biological processes accomplished by interaction of multiple molecular activities
    - Examples: glucose transmembrane transport, DNA repair

With this information, genes can be described/annotated with multiple terms. These terms are called GO terms. Here is the basic structure of a GO term:

```
GOID: GO:0007165
Term: signal transduction
Ontology: BP
Definition: The cellular process in which a signal is conveyed to
    trigger a change in the activity or state of a cell. Signal
    transduction begins with reception of a signal (e.g. a ligand
    binding to a receptor or receptor activation by a stimulus such as
    light), or for signal transduction in the absence of ligand,
    signal-withdrawal or the activity of a constitutively active
    receptor. Signal transduction ends with regulation of a downstream
    cellular process, e.g. regulation of transcription or regulation of
    a metabolic process. Signal transduction covers signaling from
    receptors located on the surface of the cell and signaling via
    molecules located within the cell. For signaling between cells,
    signal transduction is restricted to events at and within the
    receiving cell.
Synonym: GO:0023033
Synonym: signaling pathway
Synonym: signalling pathway
Synonym: signaling cascade
Synonym: signalling cascade
Secondary: GO:0023033
```

**Elements**

- GOID - unique 7-digit identifier that is intended to be machine readable
- Term - human-readable name
- Ontology - molecular function (MF), cellular component (CC), biological process (BP)
- Definition - description of what the term means
- Synonym - alternate words closely related in meaning to term, relevant connections
- Secondary - ID created when 2 or more terms are identical in meaning and so are merged into a single term. Secondary ID preserves the excess GO terms

There is so much more information available with these terms, including relationships to other genes/gene products, graphical representation of related GO terms, and much more, but that is beyond the scope of this analysis.

There are many different methdods to assigning GO terms to genes/gene products of interest. This workflow will focus on InterProScan, BLAST2GO, and Uniprot. These tools were chosen because they all utilize different databases for annotation. Additionally, they will all be run in different places: InterProScan on HPC, BLAST2GO on local computer, and Uniprot online. Results can be compared across tools to assess how each performed.

#### i) Run InterProScan to obtain GO terms

[InterProScan](https://www.ebi.ac.uk/interpro/) is a website/software that provides functional analysis of proteins by using predictive models that look across protein databases to find homologies with query proteins. If a homology is identified in one of the databases, the information about that homology is used to assign the query protein a GO term.

Note: Many fasta files willl use an asterisk to denote a STOP codon. InterProScan does not accept special characters within the sequences, so I removed them prior to running the program using the code below:

```
cp Pver_proteins_names_v1.0.faa /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/BLAST-GO-KO/InterProScan/
sed -i 's/*//g' Pver_proteins_names_v1.0.faa
```

InterProScan utilizes several member databases to enhance the chance of obtaining good protein info:

- Structural domains - Gene3D, Superfamily
- Functional annotation of families/domains - PIRSF, TIGR, Panther, Pfam, SMART, Prints, Hamap, ProSite
- Protein features - ProSite
- Prediction of conserved domains - ProDom

##### a) Run InterProScan on Andromeda

```
# On Andromeda

# Load module
module load InterProScan/5.52-86.0-foss-2021a
module load Java/11.0.2
java -version

# Run InterProScan
interproscan.sh -version
interproscan.sh -f XML -i Pver_proteins_names_v1.0.faa -b ./Pver.interpro.20210927  -iprlookup -goterms -pa
interproscan.sh -mode convert -f GFF3 -i ./Pver.interpro.20210927.xml -b ././Pver.interpro.20210927

```

**Options**

- ```version``` - display version number of IPS and java
- ```f``` - output format
- ```i``` - path to input file
- ```b``` - base output name
- ```iprlookup``` - provides mapping from matched member databases to InterPro entries
- ```goterms``` - provides mapping to Gene Ontology
- ```pa``` - provides mapping from matches to pathway info, which is based on matched mantually curated InterPro enteries
- ```mode    ``` - convert - change file format

##### b) View output file

Example of part of IPS output file:

seqid | source | type | start | end | score | strand | phase | attributes
--- | --- | --- | --- | --- | --- | --- | --- | --- |
Pver_evm.model.Segkk4293_pilon.5     | Pver | protein_match | 158 | 257 | 1.1E-19 | + | . | date=26-09-2020;Target=Pver_evm.model.Segkk4293_pilon.5 158 257;Ontology_term="GO:0006811","GO:0016021";ID=match$3_158_257;signature_desc=Neurotransmitter-gated ion-channel transmembrane region;Name=PF02932;status=T;Dbxref="InterPro:IPR006029"

**Column names**

- seqid - query sequence  id
- source - algorithm or software that generated feature
- type - type of feature
- start - start coordinates of feature
- end - end coordinates of feature
- score - similar to evalue or pvalue
- strand - + for positive strand, - for negative strand
- phase - ?
- attributes - list of feature attributes in format tag=value. There can be multiple tag=value attributes
    - Date
    - Target - target sequence for analysis and start/stop coordinates
    - Ontology - link between feature and ontology databses, like Gene Ontology (ie GO terms)
    - ID - feature ID, required to be unique for features that have 'children' (gene or mRNA)
    - Signature desc
    - Name - feature name, does not have to be unique
    - Status
    - Dbxref - link between feature and other databases, like InterPro

I noticed some inconsistencies with my output GFF3 output table and Jill's GFF3 output table.

On my output table I noticed there were a lot of extra pathway annotations that occurred after the Dbxref="InterPro:IPR000504" section of the output. Compared to Jill's that did not have this. While their is nothing wrong with this, it takes up a lot of space and makes the file much longer.

My Output GFF3 Format:

![](https://raw.githubusercontent.com/daniellembecker/DanielleBecker_Lab_Notebook/master/images/interproscan.gff3.example.db.jpg)

Jill's Output GFF3 Format:

![](https://raw.githubusercontent.com/daniellembecker/DanielleBecker_Lab_Notebook/master/images/interproscan.gff3.output.ja.jpg)

Since this information was not needed for downstream analysis, I removed it using the code below:

```
sed -e 's#^\(.*Dbxref="InterPro:[A-Z0-9]*"\).*#\1#' Pver.interpro.20210927.gff3 > Pver.interpro.20210927-smaller.gff3

```

```
I then checked to make sure nothing was deleted by this change:

1. looked at number of gene names

zgrep -c "^>" Pver.interpro.20210927.gff3

455630

zgrep -c "^>" Pver.interpro.20210927-smaller.gff3

455630

2. looked at number of Dbxref= references

grep -c 'Dbxref=' Pver.interpro.20210927.gff3

257193

grep -c 'Dbxref=' Pver.interpro.20210927-smaller.gff3

257193

3. looked at number of Ontology_term= references

grep -c 'Ontology_term=' Pver.interpro.20210927.gff3

127451

grep -c 'Ontology_term=' Pver.interpro.20210927-smaller.gff3

127451

4. looked at number of status=T references

grep -c 'status=T' Pver.interpro.20210927.gff3

429951

grep -c 'status=T' Pver.interpro.20210927-smaller.gff3

429951

```


**Full Andromeda Script:**
**Took ~16 hours to complete**

```
Pver_InterProScan.sh:

#!/bin/bash
#SBATCH --job-name="InterProScan"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="interproscan_out_error"
#SBATCH --output="interproscan_out"
#SBATCH --exclusive

cd /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/BLAST-GO-KO/InterProScan/

echo "START $(date)"

# Load module
module load InterProScan/5.52-86.0-foss-2021a
module load Java/11.0.2
java -version

# Run InterProScan
interproscan.sh --cpu $SLURM_CPUS_ON_NODE ...
interproscan.sh -version
interproscan.sh -f XML -i Pver_proteins_names_v1.0.faa -b Pver.interpro.20210927  -iprlookup -goterms -pa
interproscan.sh -mode convert -f GFF3 -i Pver.interpro.20210927.xml -b Pver.interpro.20210927

# -i is the input data
# -b is the output file base
# -f is formats
# -iprlookup enables mapping
# -goterms is GO Term
# -pa is pathway mapping
# -version displays version number

echo "DONE $(date)"

```
```
Script error:
ModuleCmd_Load.c(213):ERROR:105: Unable to locate a modulefile for 'InterProScan/5.46-81.0-foss-2019b'
openjdk version "11.0.2" 2019-01-15
OpenJDK Runtime Environment 18.9 (build 11.0.2+9)
OpenJDK 64-Bit Server VM 18.9 (build 11.0.2+9, mixed mode)
/var/spool/slurmd/job1931747/slurm_script: line 21: interproscan.sh: command not found
/var/spool/slurmd/job1931747/slurm_script: line 22: interproscan.sh: command not found
/var/spool/slurmd/job1931747/slurm_script: line 23: interproscan.sh: command not found

Error searching for module:

ModuleCmd_Load.c(213):ERROR:105: Unable to locate a modulefile for 'InterProScan/5.52-86.0-foss-2019b'
```

```
#Kevin had to download the background packages on bluewaves and suggested to use andromeda for faster running times, ran on andromeda and it worked

Submitted batch job 88966

#noticed in my output script that it mentioned a new version/module of InterProtScan available, had Kevin Bryan update this to InterProScan/5.52-86.0-foss-2021a

#also Kevin suggested adding these two additions to the above code to speed up the process:

#SBATCH --exclusive
interproscan.sh --cpu $SLURM_CPUS_ON_NODE ...

#added them above and submitted

Submitted batch job 89016
```

##### c) Secure-copy output file to local computer

```
# From a new terminal window (ie not Andromeda or remote server)

scp danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/BLAST-GO-KO/InterProScan/Pver.interpro.20210927-smaller.gff3 /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/Functional_Annotation/InterProScan/

```

#### ii) Run BLAST2GO to obtain GO terms

[BLAST2GO](https://www.blast2go.com) (B2G) is a bioinformatics tools for functional annotation and analysis of gene or protein sequences. It was originally developed to provide a user-friendly interface for GO annotation and now hosts many different functional annotation tools. The B2G software makes it possible to generate annotation without requiring writing any code. While B2G can do an extensive array of functions, this analysis primarily utilizes the GO mapping and annotation functions.

##### a) Download BLAST2GO to personal computer and activate the Basic subscription plan.

The B2G application can be downloaded [here](https://www.blast2go.com/blast2go-pro/download-b2g). B2G is available for Mac, Windows, and Linux systems. 2GB of RAM is recommended. Additionally, Internet connection is required to use most application features.

Register for B2G Basic [here](https://www.blast2go.com/b2g-register-basic). B2G Basic is free and includes the necessary features for this analysis. Registering will generate an activation key, which be put into the B2G software. You must be a part of some research institution to obtain B2G Basic.

##### b) Load the XML files generated from DIAMOND BLAST.

While B2G has the ability to run BLAST, this analysis prefers to use results from DIAMOND because of its high performance capability and senstivity. The XML file generated from DIAMOND BLAST will be used here.

To load the file, go to File<Load<Load Blast results<Load Blast XML (Legacy)

Once the file is loaded, a table loads with info about those BLAST results (nr, Tags, SeqName, Description, Length, Hits, e-Value, and sim mean). All of the cells should be orange with Tags that say BLASTED. This indicates that these sequences have only been blasted, not mapped or annotated.

##### c) Map GO terms

Mapping is the process of retrieving GO terms associated with the Description obtained by the DIAMOND BLAST search. Several mapping {steps} occur:

- BLAST result accessions are used to retrieve gene names from NCBI and those gene names are used to search the GO database.
- BLAST result accessions are used to retrieve protein information (with GO terms already annotated) through UniProt, which uses the databases SD, UniProt, Swiss-Prot, TrEMBL, RefSeq, GenPept and PDB.
- BLAST result accessions are directly searched in GO database.

To map results, select the mapping icon (white circle with green circle inside) at the top of the screen. Then select Run Mapping. A box will open up; don't change anything, click run. Depending on the number of BLAST results, mapping could take hours to days. B2G will continue running if the computer is in sleep mode. Mapping status can be checked under the Progress tab in the lower left box. If mapping retrieved a GO term that may be related to a certain sequence, that sequence row will turn light green.

##### d) Annotate GO terms

Now that GO terms have been retrieved, the annotation process will select GO terms from the GO pool obtained through mapping and assign them to query sequences. Several annotation steps occur:

- For all found GO terms, an annotation rule (AR) is applied. The rule seeks to find the most specific annotations with a certain level of reliability and can be adjusted in the settings prior to running annotation.
- If a candidate GO term is found, an annotation score is calculated that weights the highest hit similarity of that candidate GO term by its evidence code (from mapping step). The candidate GO term will only be assigned to a sequence if it passes a certain threshold based on calculations above. For more info about the math behind annotation, go [here](http://docs.blast2go.com/user-manual/gene-ontology-annotation/).

To annotate, select the annot icon (white circle with blue circle inside) at the top of the screen. Then select Run Annotation. A box will open up; don't change anything unless thresholds need to be adjusted. If no changes are necessary, click next through the boxes until the final one and click run. Depending on the mapping results, annotating could take hours to days. B2G will continue running if the computer is in sleep mode. Annotation status can be checked under the Progress tab in the lower left box. If a GO term has been successfully assigned to a sequence, that sequence row will turn blue.

##### e) Export annotated sequences and info

To export the file with B2G annotated genes, go to File<Export<Export as Table. A box will open up,  and save as a text file. Re-save as csv as needed.

##### f) View output file

seqName | top_hit | length | evalue | simMean | GO.ID | GO_names
--- | --- | --- | --- | --- | --- | --- |
Pver_evm.model.Segkk0_pilon.11 | EDO40119.1 | 1779 | 1.9e-199 | 68.38 | P:GO:0007166; P:GO:0007275; C:GO:0016020; C:GO:0120025 | P:cell surface receptor signaling pathway; P:multicellular organism development; C:membrane; C:plasma membrane bounded cell projection

**Column names**

- seqName - query sequence
- top_hit - top match from DIAMOND BLAST
- length - alignment length
- evalue - number of expected hits of similar quality that could be found by chance (similar to a pvalue). The smaller the evalue, the better the match!
- simMean - mean similarity between query and top match
- GO.ID - gene ontology mapping results, GO terms and evidence codes
- GO_names - annotated GO terms


#### iii) Run Uniprot to obtain GO terms

[Uniprot](https://www.uniprot.org) (Universal Protein Resource) provides information for protein sequence and annotation data. It maintains several protein databases:

- UniProt Knowledgebase (UniProtKB) - central collection hub for functional protein info and annotation
- UniProt Reference Clusters (UniRef) - provides clustering from UniProtKB to ensure complete coverage of sequences while masking redundant sequences
- UniProt Archive (UniParc) - database of publicly available protein sequences

In this analysis, Uniprot uses BLAST input to search against its protein databases.

##### a) Make a list of identifiers found by DIAMOND BLAST.

Uniprot uses the .tab file generated from DIAMOND BLAST as input. The column 'sseqid' is the primary input to Uniprot, as it can use that BLAST input to find protein matches. Because UniProt is a website and lacks proper storage and RAM, it cannot handle an entire .tab DIAMOND BLAST file. To ensure that UniProt reads all inputs, subset files so that a file only has ~2000-3000 lines per file in Terminal.

```
# Check how many lines in the full file
wc -l Pver_annot.tab
    21606 Pver_annot.tab

#https://kb.iu.edu/d/afar (for instructions on split command)
#use split command to split into multiple files that have 2,000 designated lines each

split -l 2000 Pver_annot.tab tab

tabaa  tabad  tabag  tabaj
tabab  tabae  tabah  tabak
tabac  tabaf  tabai
```

##### b) Navigate to Uniprot [Retrieve/ID mapping page](https://www.uniprot.org/uploadlists/) and under Provide your identifiers, click 'upload your own file' and upload a subsetted tab file.

UniProt will be able to process the smaller/subsetted files. However, each subsetted .tab file will need to be put in as input one at a time.

##### c) Under 'Select options', choose the place where the identifiers were generated (From) and what database to compare to (To).

In this analysis, the identifiers were generated from the EMBL/GenBank/DDBJ CDS option (aka NCBI/BLAST) and they will be compared to (UniProtKB). Once that is finished, hit submit. It may take a few minutes). If the mapping fails and the page displays 'Service Unavailable', go back and try again.

##### d) Select appropiate columns to include in table.

If mapping was successful, the screen should have something like ```X out of Y EMBL/GenBank/DDBJ CDS identifiers were successfully mapped to X UniProtKB IDs in the table below```. There will be a table below with the results. Select the 'Columns' tab right above the table. This will open a window to pick more column options so more information can be included in the table. Many of the column options don't apply to this analysis.

Under 'Names & Taxonomy', Entry name, Gene names, Organism, and Protein Names should already be selected. Under 'Sequences', Length should already be selected. Under 'Miscellanous', the gold paper with a star and the blue paper symbols should already be selected. If they are not, select those columns. Under 'Gene Ontology', select Gene ontology (GO) and Gene ontology IDs. Under 'Genome Annotation', select KEGG. Once selections are completed, scroll back to the top of the page and hit Save in the upper-right hand corner.

##### e) Save files of interest

To save the main table, click Download in the top left of the table. Select Download all, Uncompressed, and Tab-seperated as Format. Click Go.

To save the unmapped identifiers, click on the 'Click here to download the 6281 unmapped identifiers'. To save UniParc results (if any are found), click on the UniParc link under the main header of ```X out of Y EMBL/GenBank/DDBJ CDS identifiers were successfully mapped to X UniProtKB IDs in the table below```. To finish downloading UniParc results, click Download in the top left of the table. Select Download all, Uncompressed, and Tab-seperated as Format. Click Go.


##### f) View output files

This is an example of Uniprot output .tab file:

my_list | top_hit | uniprotkb_entry | status | protein_names | gene_names | organism | length | go_ids | gene_ontology | kegg |
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
EDO36971.1 | A7SH22 | A7SH22_NEMVE | unreviewed | Predicted protein (Fragment) | v1g118173 | Nematostella vectensis (Starlet sea anemone) | 239 | GO:0004888; GO:0005230; GO:0005887; GO:0007165; GO:0007268; GO:0030594; GO:0034220; GO:0042391; GO:0043005; GO:0045202; GO:0050877 | integral component of plasma membrane [GO:0005887]; neuron projection [GO:0043005]; synapse [GO:0045202]; extracellular ligand-gated ion channel activity [GO:0005230]; neurotransmitter receptor activity [GO:0030594]; transmembrane signaling receptor activity [GO:0004888]; chemical synaptic transmission [GO:0007268]; ion transmembrane transport [GO:0034220]; nervous system process [GO:0050877]; regulation of membrane potential [GO:0042391]; signal transduction [GO:0007165]  | nve:5508437 |

**Column names**

- my_list - query id (generated from DIAMOND BLAST)
- top_hit - top match from Uniprot databases
- uniprotkb_entry - official Uniprot entry for top_hit
- status - entry status; indicates if entry has been manually annotated and reviewed (reviewed = SwissProt section of Uniprot, unreviewed = computer-annotated TrEMBL section)
- protein_names - list of all names for a particular protein (Predicted protein = unreviewed = omputer-annotated TrEMBL section)
- gene_names - name of the genes that code for particular protein sequence
- organism - name of organism that is source of protein information
- length - alignment length
- go_ids - gene ontology mapping results, GO terms and evidence codes
- gene_ontology - annotated GO terms
- ko - K number used to reconstruct biological pathways in KEGG mapper (outside scope of this analysis)
- kegg - gene identifier used to search biological pathways in KEGG mapper (outside scope of this analysis)

### Step 4: Merge all information for full annotation

Follow instructions/code [here](https://github.com/JillAshey/FunctionalAnnotation/blob/main/RAnalysis/acerv_annot_compile.R) to compile all annotation information in R.
