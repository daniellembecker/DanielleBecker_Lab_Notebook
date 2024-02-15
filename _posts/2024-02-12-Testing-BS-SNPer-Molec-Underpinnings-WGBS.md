---
layout: post
title: Testing BS-SNPer on WGBS data
date: '2024-02-12'
categories: Analysis
tags: WGBS, BS-SNPer
---

Code is inspired by [this post](https://github.com/RobertsLab/project-gigas-oa-meth/blob/master/code/07-BS-SNPer.ipynb), [this post](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/_posts/2022-11-21-Testing-BS-SNPer-on-WGBS-data.md), and [this post](https://github.com/hputnam/Geoduck_Meth/blob/master/code/10-ct-snp.sh).

In this notebook post, I'll use [BS-Snper](https://github.com/hellbelly/BS-Snper) to call SNP variants from the *Pocillopora meandrina* and *grandis* WGBS data from the Molecular Underpinnings project. Adult colonies originated from a control and nutrient enriched reef site in Moorea, French Polynesia. DNA was extracted from whole tissue for WGBS. Reads were aligned to the *P. verrucosa* genome through the nf-ccore methylseq pipeline.

Other lab members, Kevin Wong and Emma Strand have tried to run BS-SNPer on their data but keep encountering this error: Too many characters in one row! Try to split the long row into several short rows (fewer than 1000000 characters per row).
Error! at /opt/software/BS-Snper/1.0-foss-2021b/bin/BS-Snper.pl line 110.

Kevin also reached out to URI research computing and there is now an updated BS-SNPer module that may have solved this error that I will try.

As we have run the same pipeline and similar extractions from coral tissue samples, we think it may be an issue with the .bam output files from the nf-core methylseq pipeine. I am going to try an approach with my samples, first with the methods used by Steven Roberts in this [post](https://github.com/hputnam/Geoduck_Meth/blob/master/code/10-ct-snp.sh) and then an approach used by Kevin Wong [here](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/_posts/2022-11-21-Testing-BS-SNPer-on-WGBS-data.md).

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

`samtools view merged-sorted-deduplicated.bam | head `

```

```

## 2b. Identify SNPs

Options for the script are found [here](https://github.com/hellbelly/BS-Snper/blob/master/README.txt) and below.

```
--fa: Reference genome file in fasta format
--input: Input bam file (I'm using deduplicated sorted bams)
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
```

You can run BS-SNPer in Linux or MAC OS, using the command like:

```
perl BS-Snper.pl <sorted_bam_file> --fa <reference_file> --output <snp_result_file> --methcg <meth_cg_result_file> --methchg <meth_chg_result_file> --methchh <meth_chh_result_file> --minhetfreq 0.1 --minhomfreq 0.85 --minquali 15 --mincover 10 --maxcover 1000 --minread2 2 --errorate 0.02 --mapvalue 20 >SNP.out 2>ERR.log
```

**Attention**: Both of the input and output file arguments should be passed to BS-SNPer in the form of absolute paths.

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
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output

module load BS-Snper/1.0-foss-2021b

perl /opt/software/BS-Snper/1.0-foss-2021b/bin/BS-Snper.pl \
data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/bismark_deduplicated/BS_SNPer_merged.bam \
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




Following Kevins script:

`nano /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/bs_snper_merged_kevin.sh`

```bash
#!/bin/bash
#SBATCH --job-name="BS_snper"
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output

module load BS-Snper/1.0-foss-2021b

perl /opt/software/BS-Snper/1.0-foss-2021b/bin/BS-Snper.pl \
--fa /data/putnamlab/dbecks/Becker_E5/Becker_RNASeq/data/refs/Pverr/Pver_genome_assembly_v1.0.fasta \
--input /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/bismark_deduplicated/BS_SNPer_merged.bam \
--output /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/SNP-candidates.txt \
--methcg /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/CpG-meth-info.tab \
--methchg /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/CHG-meth-info.tab \
--methchh /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/CHH-meth-info.tab \
--mincover 5 \
> /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_outputSNP-results.vcf 2> /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/Becker_WGBS/BS-SNPer/merged_SNP_output/merged.ERR.log

```

`sbatch /data/putnamlab/dbecks/Becker_E5/WGBS_Becker_E5/scripts/bs_snper_merged_kevin.sh`
