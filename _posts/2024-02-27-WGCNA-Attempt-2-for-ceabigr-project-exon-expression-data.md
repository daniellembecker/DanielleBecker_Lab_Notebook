---
layout: post
title: WGCNA Attempt 2 for ceabigr project exon expression data
date: '2024-02-27'
categories: CEABIGR
tags: GeneExpression Oyster WGCNA R
---

Today I played around with another attempt at WGCNA analyses to look at exon expression data for the ceabigr project. 

This time, we are trying data formatted with `sample_foldchange` in rows and `gene` in columns.  

See my [previous post here](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/WGCNA-Attempt-1-for-ceabigr-project-exon-expression-data/) for more explanation.  

# Coding and results 

```
#library(kableExtra)
# library(DESeq2)
# library(pheatmap)
# library(RColorBrewer)
# library(data.table)
#library(DT)
# library(Biostrings)
#library(methylKit)
library(WGCNA)
library(data.table) # for data manipulation
library(tidyverse)

knitr::opts_chunk$set(
  echo = TRUE,         # Display code chunks
  eval = TRUE,         # Evaluate code chunks
  warning = TRUE,     # Hide warnings
  message = TRUE,     # Hide messages
  fig.width = 6,       # Set plot width in inches
  fig.height = 4,      # Set plot height in inches
  fig.align = "center" # Align plots to the center
)
```
 
### Read in data 

Read in data that now has sample_fold in rows and genes in columns with values representing fold change from the first exon relative to each subsequent exon.  

Running for just females for now.  

```
#setwd("Projects/ceabigr")

#datExpr <- fread("output/68-female-exon-fold/logfc.txt")
datExpr <- fread("output/72-exon-data-rfmt/female_exon_tf.csv")
```

Set row names and remove character column. 

```
rownames(datExpr) <- datExpr$SampleID_fold
sample_vector<-datExpr$SampleID_fold
datExpr <- datExpr[ , -1, with = FALSE]
```

Remove rows that are a sum of 0. These are the fold 1 rows.  

```
# Find rows with sum not equal to 0
nonZeroSumRows <- rowSums(datExpr, na.rm=TRUE) != 0

# Subset the matrix using logical indexing
datExpr <- datExpr[nonZeroSumRows, ]
```

Remove columns that have na's (removed genes not present in every sample).  

```
datExpr<-datExpr %>%
    select_if(~ !any(is.na(.)))
```
We had 13,280 genes before, and now have 11,270 genes. 

Choose a soft power. 

```
powers = c(1:35)
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], sft$fitIndices[,3], pch = 19, xlab="Soft Threshold (power)", ylab="scale free topology model fit", type="n")
text(sft$fitIndices[,1], sft$fitIndices[,3], labels=powers, cex=0.5)
abline(h = 0.9, col = "red")

```

It's having a problem estimating variance to generate the soft powers... It doesn't reach the 0.9 threshold. We will have to return to this... Unsure what it means.  

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/oysters/ceabigr/soft_power.png?raw=true)  

### Blockwise modules network construction 

Run blockwise modules with a signed network. 

```

picked_power =  10
temp_cor <- cor       
cor <- WGCNA::cor                                             # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(datExpr,                         # <= input here

                          # == Adjacency Function ==
                          power = picked_power,               # <= power here
                          networkType = "signed",

                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 1000,                  
                          maxBlockSize = 10000,

                          # == Module Adjustments ==
                          mergeCutHeight = 0.05,
                          reassignThreshold = 1e-6,
                          minCoreKME = 0.5,
                          minKMEtoStay = 0.3,

                          # == TOM == Archive the run results in TOM file (saves time) but it doesn't save a file
                          saveTOMs = F,
                          saveTOMFileBase = "ER",

                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

cor <- temp_cor     # Return cor function to original namespace

# Identify labels as numbers 
mergedColors = netwk$colors
# Plot the dendrogram and the module colors underneath

table(mergedColors)

membership<-as.data.frame(mergedColors)

membership$gene<-rownames(membership)

names(membership)<-c("module", "gene")
```

### Plot module eigengene level

Next plot eigengene levels 

```
#Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = mergedColors, softPower = 8)
MEs = MEList$eigengenes

ncol(MEs) #How many modules do we have now?

table(mergedColors)
table<-as.data.frame(table(mergedColors))
table
```

The table shows number of genes in each module.

  mergedColors Freq
1            0 5520
2            1 3511
3            2 2239

There are three modules.  

Plot module expression across exon location. 

```
head(MEs)
names(MEs)
Strader_MEs <- MEs
#Strader_MEs$exon <- c("2", "3", "4", "5", "6")
cleaned_vector <- sample_vector[!grepl("fold1", sample_vector)]
Strader_MEs$sample <- cleaned_vector
head(Strader_MEs)

Strader_MEs <- separate(Strader_MEs, col = sample, into = c("sample", "fold"), sep = "_", remove = FALSE)
head(Strader_MEs)
```

```
plot_MEs<-Strader_MEs%>%
  gather(., key="Module", value="Mean", ME0:ME2)
```

First, assign treatment by sample ID. 

```
meta<-read_csv("data/adult-meta.csv")

plot_MEs$treatment<-meta$Treatment[match(plot_MEs$sample, meta$OldSample.ID)]
```

Plot module expression across exon location. 

```
library(ggplot2)
library(tidyverse)

expression_plot<-plot_MEs%>%
  group_by(Module, fold) %>%
  
  ggplot(aes(x=fold, y=Mean, color=treatment)) +
  facet_wrap(~Module)+
  geom_hline(yintercept = 0, linetype="dashed", color = "grey")+
  geom_point() +
  scale_color_manual(values=c("gray", "darkred"))+
  #geom_line(group=1)+
  #ylim(-0.5,1) +
  ylab("Mean Module Eigenegene") +
  xlab("Exon fold change")+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA)); expression_plot

ggsave(plot=expression_plot, filename="output/69-wgcna/module-expression.png", width=10, height=6)
```

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/oysters/ceabigr/module-expression_attempt2.png?raw=true) 

It looks like Module 0 is genes that peak in the middle (exon 3-4); Module 1 are genes that have higher expression of earlier exons (exon 2-3) and Module 2 are genes that increase with higher expression of later exons (exons 5-6). 

# Next steps 

I will next need to figure out which genes change module assignment by treatment. Currently, I don't have a way to do this because each gene is assigned to one module. We do not have separate columns for the gene value from control and that same gene from treatment as we did last time.  



