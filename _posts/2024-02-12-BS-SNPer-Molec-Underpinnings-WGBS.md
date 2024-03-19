---
layout: post
title: BS-SNPer on WGBS data
date: '2024-02-12'
categories: Analysis
tags: WGBS, BS-SNPer
---

Code is inspired by [this post](https://github.com/RobertsLab/project-gigas-oa-meth/blob/master/code/07-BS-SNPer.ipynb), [this post](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/_posts/2022-11-21-Testing-BS-SNPer-on-WGBS-data.md), and [this post](https://github.com/hputnam/Geoduck_Meth/blob/master/code/10-ct-snp.sh).

In this notebook post, I'll use [BS-Snper](https://github.com/hellbelly/BS-Snper) to call SNP variants from the *Pocillopora meandrina* and *grandis* WGBS data from the Molecular Underpinnings project. Adult colonies originated from a control and nutrient enriched reef site in Moorea, French Polynesia. DNA was extracted from whole tissue for WGBS. Reads were aligned to the *P. verrucosa* genome through the nf-ccore methylseq pipeline.

Other lab members, Kevin Wong and Emma Strand have tried to run BS-SNPer on their data but keep encountering this error: Too many characters in one row! Try to split the long row into several short rows (fewer than 1000000 characters per row).
Error! at /opt/software/BS-Snper/1.0-foss-2021b/bin/BS-Snper.pl line 110.

I am going to try an approach with my samples, with the methods used by Steven Roberts in this [post](https://github.com/hputnam/Geoduck_Meth/blob/master/code/10-ct-snp.sh).

# 1. Set working directory

`cd /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS`

`mkdir BS-SNPer`

# 2. Identify SNP variants

I will identify variants in individual files, as well as SNPs across all samples.

## 2a. Merge SNP variants

First I need to sort all the deduplicated bam files

`cd /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/bismark_deduplicated`

`nano /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/sort_bam.sh`

```bash
#!/bin/bash
#SBATCH --job-name="sort_bam"
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/bismark_deduplicated

# load modules needed

module load SAMtools/1.16.1-GCC-11.3.0

for f in *.deduplicated.bam
do
  STEM=$(basename "${f}" _R1_001_val_1_bismark_bt2_pe.deduplicated.bam)
  samtools sort "${f}" \
  -o "${STEM}".deduplicated_sorted.bam
done
```

`sbatch /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/sort_bam.sh`

`Submitted batch job 302471`

To identify SNPs across all samples, I need to merge my samples, then use that as the input file for BS-Snper.

`nano /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/merge_bam.sh`

```bash
#!/bin/bash
#SBATCH --job-name="merge_snp"
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/bismark_deduplicated

# load modules needed

module load SAMtools/1.16.1-GCC-11.3.0

# Merge Samples with SAMtools

samtools merge \
BS_SNPer_merged.bam \
*sorted.bam

```

`sbatch /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/merge_bam.sh`

`Submitted batch job 302931`


View output file header

## This feature is only in the older module of samtools

`module load SAMtools/1.9-foss-2018b `

`samtools view BS-SNPer_merged.bam | head `

```
GWNJ-1013:120:GW201202000:2:2318:21766:34992_1:N:0:GGCTTAAG+GGTCACGA	163	Pver_Sc0000000_size2095917	28	42	130M	=	240	342	ATATCTAACACATTAAAACTTCCTAAATCAAAAAAAAATTCCCTCTACATCAAACACAACAAACAACATTTTAATAAAAACACACATATATATTTCTTACCCCATTTTCAATAATAATAACAAATCAATA	FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFF:FFFFFFFFFFFFFFFFFFFFFF:FFF:FFFFFFFFFFFFF:F,FFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFF	NM:i:24	MD:Z:7G1G1G2G16G2G0G10G5G8G10G3G0G5G25G0G1G2G2G4G1T0G0G1G0	XM:Z:.......h.z.z..h................h..hh..........x.....x........x..........h...hh.....z.........................zx.h..h..h....h..hh.hXR:Z:GA	XG:Z:GA

GWNJ-1013:120:GW201202000:2:1535:13349:28714_1:N:0:TTACAGGA+GCTTGTCA	99	Pver_Sc0000000_size2095917	42	42	24M1D106M	=	378	465	GAAATTTTTTAAATTAAGAAGGAATTTTTTTGTATTAGATATAATAGATAATATTTTGATAGGAATATGTATATATATTTTTTGTTTTATTTTTGGTGATGATGATAAGTTGGTGTTGTGTTATTTTAGT	FFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFF:FFFF:FFFFFFFFFFFFFFFF:F:FFFFF,FFF::FF,FFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFF	NM:i:31	MD:Z:4C2C0C5C9^T1C0C0C1C2C2C3C1C2C3C2C13C1C1C10C2A0C0C0C0C5C11C15C1C2C2C0	XM:Z:....h..hh.....h..........hhh.x..h..x...h.h..x...h..h.............h.z.h..........h...hhhh.....z...........h...............h.h..x..h	XR:Z:CT	XG:Z:CT

GWNJ-1013:120:GW201202000:2:1117:2908:33927_1:N:0:TCTCTACT+GAACCGCG	163	Pver_Sc0000000_size2095917	68	6	130M	=	196	259	CCCTCTACATCAAACACAACAAACAACATTTTAATAAAAACACACATATATATTTCTTACCCCATTTTCAATAATAATAACAAATCAATATTATATCACTTCAACAAATACCTTTATCAACATAAAAATA	FFF,F,FFFFFFFF::FFF:FFFFFFFF:FFFFF,FFFFFFFFF:FFF,F,F,FF,,FFFFFFF,FFFFFFFFFFFF:FFF:FF,FFF,FFFFFFF:FF,F:FFF,FFFFFF,FFF,FFF,:,FFFFFFF	NM:i:25	MD:Z:6G5G8G10G3G0G5G25G0G1G2G2G4G1T0G0G1G2G1G8G5G9G4G2G1G0	XM:Z:......x.....x........x..........h...hh.....z.........................zx.h..h..h....h..hh.h..h.h........x.....h.........x....h..h.hXR:Z:GA	XG:Z:GA

GWNJ-1013:120:GW201202000:2:2622:23312:3787_1:N:0:TTACAGGA+GCTTGTCA	163	Pver_Sc0000000_size2095917	77	21	129M	=	346	399	TCAAACACAACAAACAACATTTTAATAAAAACACACATATATATTTCTTACCCCATTTTCAATAATAATAACAAATCAATATTATGTCACTTCAACAAATACCTTTATCAACATAAAAATAACAAAATT	FFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:	NM:i:24	MD:Z:3G8G10G3G0G5G25G0G1G2G2G4G1T0G0G1G2G10G5G9G4G2G1G2G5	XM:Z:...x........x..........h...hh.....z.........................zx.h..h..h....h..hh.h..h.H........x.....h.........x....h..h.h..z.....	XR:Z:GA	XG:Z:GA
```

## 2b. Identify SNPs

Options for the script are found [here](https://github.com/hellbelly/BS-Snper/blob/master/README.txt) and instructions are taken from BS-SNPer GitHub descriptions.

**Usage:**

You can run BS-SNPer in Linux or MAC OS, using the command:

`perl BS-Snper.pl <sorted_bam_file> --fa <reference_file> --output <snp_result_file> --methcg <meth_cg_result_file> --methchg <meth_chg_result_file> --methchh <meth_chh_result_file> --minhetfreq 0.1 --minhomfreq 0.85 --minquali 15 --mincover 10 --maxcover 1000 --minread2 2 --errorate 0.02 --mapvalue 20 >SNP.out 2>ERR.log`

**Attention:**

Both of the input and output file arguments should be passed to BS-SNPer in the form of absolute paths. Providing absolute paths ensures that BS-SNPer can accurately locate and process the specified input and output files, regardless of the current working directory or other factors.

**Options:**

	--fa: Reference genome file in fasta format
	--output: Temporary file storing SNP candidates
	--methcg: CpG methylation information
	--methchg: CHG methylation information
	--methchh: CHH methylation information
	--minhetfreq: Threshold of frequency for calling heterozygous SNP
	--minhomfreq: Threshold of frequency for calling homozygous SNP
	--minquali: Threshold of base quality
	--mincover: Threshold of minimum depth of covered reads
	--maxcover: Threshold of maximum depth of covered reads
	--minread2: Minimum mutation reads number
	--errorate: Minimum mutation rate
	--mapvalue: Minimum read mapping value
	SNP.out: Final SNP result file
	ERR.log: Log file

## Write script for BS-SNPer approach for our data

`nano /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/bs_snper_merged_steven.sh`

Following Stevens script:

```bash
#!/bin/bash
#SBATCH --job-name="BS_snper"
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output

module load BS-Snper/1.0-foss-2021b

perl /opt/software/BS-Snper/1.0-foss-2021b/bin/BS-Snper.pl \
/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/bismark_deduplicated/BS_SNPer_merged.bam \
--fa /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Pverr/Pver_genome_assembly_v1.0.fasta \
--output /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/SNP-candidates.out \
--methcg /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/CpG-meth-info.tab \
--methchg /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/CHG-meth-info.tab \
--methchh /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/CHH-meth-info.tab \
--minhetfreq 0.1 \
--minhomfreq 0.85 \
--minquali 15 \
--mincover 10 \
--maxcover 1000 \
--minread2 2 \
--errorate 0.02 \
--mapvalue 20 \
> /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/SNP-results.vcf 2>/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/merged.ERR.log

```

`sbatch /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/bs_snper_merged_steven.sh`

`Submitted batch job 302954`

Current issue:

```
Unknown option: input
FLAG: 1
refSeqFile = /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Pverr/Pver_genome_assembly_v1.0.fasta.
bamFileName = /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/bismark_deduplicated/BS_SNPer_merged.bam.
snpFileName = /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/SNP-candidates.txt.
methCgFileName = /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/CpG-meth-info.tab.
methChgFileName = /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/CHG-meth-info.tab.
methChhFileName = /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/CHH-meth-info.tab.
vQualMin = 15.
nLayerMax = 1000.
vSnpRate = 0.100000.
vSnpPerBase = 0.020000.
mapqThr = 20.
Too many characters in one row! Try to split the long row into several short rows (fewer than 1000000 characters per row).
Error! at /opt/software/BS-Snper/1.0-foss-2021b/bin/BS-Snper.pl line 110.
```

Possible reason for error: The limit for the number of characters per row is 1,00,000 and our input fasta file has a line that is  2,095,918. Research computing patched the module and BS-SNPer.pl script to accept 3,000,000. I reran the Steven approach script.

Python script used to determine number of characters:

`python /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/max_ln.py /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Pverr/Pver_genome_assembly_v1.0.fasta`

`sbatch /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/bs_snper_merged_steven.sh`

`Submitted batch job 303588`

Received another error:

```
Unknown option: input
FLAG: 1
refSeqFile = /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Pverr/Pver_genome_assembly_v1.0.fasta.
bamFileName = /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/bismark_deduplicated/BS_SNPer_merged.bam.
snpFileName = /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/SNP-candidates.txt.
methCgFileName = /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/CpG-meth-info.tab.
methChgFileName = /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/CHG-meth-info.tab.
methChhFileName = /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/CHH-meth-info.tab.
vQualMin = 15.
nLayerMax = 1000.
vSnpRate = 0.100000.
vSnpPerBase = 0.020000.
mapqThr = 20.
Too many chromosomes!
Error! at /opt/software/BS-Snper/1.0-foss-2021b/bin/BS-Snper.pl line 110.
```

Possible reason for this error: The .pl script had a limit of 10k chromosomes. Research computing patched it to accept 30k (same multiplier as the line length). I then reran the script.

`sbatch /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/bs_snper_merged_steven.sh`

`Submitted batch job 303600`

IT IS RUNNING YAY!

Another error after the script ran:

```
At the end of `merged.ERR.log` it says it failed to call a copy of samtools it expected to be local to the BS-Snper install.

```

Research computing modified the install recipe to copy it to where it's expecting to find it.


`sbatch /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/bs_snper_merged_steven.sh`

`Submitted batch job 303619`

The script ran without errors! The output .vcf files match the SNP output file standard VCF format.

The output files include an SNP output file and a methylation output file.
The SNP output file has a standard VCF format.
The methylation output file has a tab-separated format same as MethylExtract (https://sourceforge.net/projects/methylextract/files/MethylExtract_1.4/ManualMethylExtract.pdf/download):

1. CHROM: Chromosome.
2. POS: Sequence context most 5’ position on the Watson strand (1-based).
3. CONTEXT: Sequence contexts with the SNVs annotated using the IUPAC nucleotide ambiguity code (referred to the Watson strand).
4. Watson METH: The number of methyl-cytosines (referred to the Watson strand).
5. Watson COVERAGE: The number of reads covering the cytosine in this sequence context (referred to the Watson strand).
6. Watson QUAL: Average PHRED score for the reads covering the cytosine (referred to the Watson strand).
7. Crick METH: The number of methyl-cytosines (referred to the Watson strand).
8. Crick COVERAGE: The number of reads covering the guanine in this context (referred to the Watson strand).
9. Crick QUAL: Average PHRED score for the reads covering the guanine (referred to the Watson strand).`

View SNP-results.vcf file:

`head /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/SNP-results.vcf`

`##fileDate= 20240220
##bssnperVersion=1.1
##bssnperCommand=--fa /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Pverr/Pver_genome_assembly_v1.0.fasta   --input /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/bismark_deduplicated/BS_SNPer_merged.bam --output /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/SNP-candidates.out --methcg /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/CpG-meth-info.tab --methchg /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/CHG-meth-info.tab --methchh /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/CHH-meth-info.tab --minhetfreq  0.1 --minhomfreq  0.85 --minquali 15 --mincover 10 --maxcover 1000 --minread2 2 --errorate 0.02 --mapvalue 20
##reference=file:///data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Pverr/Pver_genome_assembly_v1.0.fasta
##Bisulfite=directional>
##contig=<ID=Pver_Sc0000000_size2095917,length=2095917>
##contig=<ID=Pver_Sc0000001_size2081954,length=2081954>
##contig=<ID=Pver_Sc0000002_size1617595,length=1617595>
##contig=<ID=Pver_Sc0000003_size1576134,length=1576134>`

Count number of lines in the SNP-results.vcf file

`wc -l /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/SNP-results.vcf`

`4,764,382`

## 2c. Intersect VCF with SNP locations and CG motif track

Use BEDtools intersectBED to determine which SNPs are present at CG sites:

You first need to make CG motif file from your origin genome to filter for just CG motifs using [EMBOSS fuzznuc](https://emboss.sourceforge.net/apps/cvs/emboss/apps/fuzznuc.html) function.

`nano /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Pverr/fuzznuc.sh`

```bash
#!/bin/bash
#SBATCH --job-name="fuzznuc"
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Pverr

module load EMBOSS/6.6.0-foss-2018b

fuzznuc \
-sequence /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Pverr/Pver_genome_assembly_v1.0.fasta \
-pattern CG \
-outfile 20240311_fuzznuc_pverr_CGmotif.gff \
-rformat gff

```

`sbatch /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Pverr/fuzznuc.sh`

`Submitted batch job 305435`

Then use BEDtools intersectBED to determine which SNPs are present at CG sites:

`nano /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/intersectbed.SNPs.sh`

```bash
#!/bin/bash
#SBATCH --job-name="BEDtools"
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=danielle_becker@uri.edu
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output

module load BEDTools/2.30.0-GCC-11.3.0

intersectBed \
-u \
-a /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/SNP-results.vcf \
-b /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Pverr/20240311_fuzznuc_pverr_CGmotif.gff \
| grep "C	T" \
> /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/CT-SNPs.tab

```

`sbatch /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/intersectbed.SNPs.sh`

`Submitted batch job 308874 `

Look at CT SNPs output file:

`head CT-SNPs.tab`

Components of output file:

1. Chromosome/Contig: The first column of the output file is likely to contain information about the chromosome or contig where the SNP is located. This information is often found in the VCF file (SNP-results.vcf).
2. Position: The second column represents the position of the SNP on the chromosome or contig.
3. ID/Name: The third column may contain an identifier or name for the SNP.
4. Reference Allele: The fourth column shows the reference allele at the specified position.
5. Alternate Allele: The fifth column shows the alternate allele at the specified position.
6. Quality Score (QUAL): There might be a column indicating the quality score of the SNP. This information is often found in the VCF file.
7. Filter Status (FILTER): This column could show whether the SNP passes certain quality filters.
8. INFO fields:
- DP: This field represents the total read depth at the position of the variant.
- ADF (Allelic Depth Forward) and ADR (Allelic Depth Reverse): These fields represent the number of reads supporting the variant allele in the forward and reverse directions, respectively.
- AD (Allelic Depth): The total number of reads supporting the variant allele (ADF + ADR).
9. FORMAT field (e.g., "GT:DP:ADF:ADR"):
- GT (Genotype): Represents the genotype of the variant. For example, "0/1" indicates a heterozygous variant.
- DP: Same as in the INFO field, representing the total read depth.
- ADF and ADR: Same as in the INFO field, representing allelic depths for forward and reverse reads, respectively.

``

Look at CT-SNPs specific tab file to see how many CT SNPs are in my dataset:

`Pver_Sc0000000_size2095917      1286    .       C       T       15      PASS    DP=144;ADF=0,0;ADR=122,22;AD=122,22;    GT:DP:ADF:ADR:AD:BSD:BSQ:ALFR   0/1:144:0,0:122,22:122,22:1,0,146,0,0,22,122,0:37,0,36,0,0,36,37,0:0.847,0.153
Pver_Sc0000000_size2095917      1610    .       C       T       62      PASS    DP=128;ADF=0,0;ADR=115,13;AD=115,13;    GT:DP:ADF:ADR:AD:BSD:BSQ:ALFR   0/1:128:0,0:115,13:115,13:0,0,21,0,0,13,115,0:0,0,37,0,0,37,37,0:0.898,0.102
Pver_Sc0000000_size2095917      5189    .       C       T       1000    PASS    DP=95;ADF=0,0;ADR=77,18;AD=77,18;       GT:DP:ADF:ADR:AD:BSD:BSQ:ALFR   0/1:95:0,0:77,18:77,18:0,0,73,0,0,18,77,0:0,0,36,0,0,37,36,0:0.811,0.189
Pver_Sc0000000_size2095917      6277    .       C       T       1000    PASS    DP=12;ADF=0,0;ADR=0,12;AD=0,12; GT:DP:ADF:ADR:AD:BSD:BSQ:ALFR   0/1:12:0,0:0,12:0,12:0,0,64,0,0,12,0,0:0,0,37,0,0,37,0,0:0.000,1.000`

`wc -l CT-SNPs.tab`

`151896 out of 4,764,382 = 3.18%`

# 3. Filter SNP variants from 10x.bed files

Download CT-SNPs.tab results file to Desktop from Andromeda

`scp -r danielle_becker@ssh3.hac.uri.edu:/data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/CT-SNPs.tab /Users/Danielle/Desktop/Putnam_Lab/Becker_E5/RAnalysis/Data/WGBS/BS-SNPer`

Filter SNPs from BCT-SNPs.tab output file from 10x .bed files created from this [workflow](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2022-12-10-P.verrucosa-WGBS-Workflow-Host.md) and remaining steps performed in this [Rscript](https://github.com/hputnam/Becker_E5/blob/master/RAnalysis/Scripts/WGBS/BS-SNPer.filter.Rmd)
