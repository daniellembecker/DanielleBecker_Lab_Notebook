---
layout: post
title: Sanger Sequencing Training
date: '2023-06-02'
categories: Tutorial
tags: Bioinformatics Molecular Sanger
---

Introduction to Sanger sequencing and training for coral species identification. 

# Software 

We will be using [Geneious Prime software](https://www.geneious.com/), which I obtained through a University of Rhode Island site license. 

# What is Sanger sequencing? 

In order to conduct marker gene analyses, we are extracting DNA from a sample and we are interested if there are nucleotide base differences (single nucleotide polymorphism). We need to identify a region in the genome that we can amplify that is conserved enough to amplify with enough variation to tell the difference between species, genotypes, and populations. 

In Sanger sequencing, we are looking at one sequence at a time. First, in the PCR we will add a primer of the conserved region, DNA polymerase (Taq polymerase able to handle high temperature), the DNA sample, and fluorescent tagged nucleotides mixed with non-tagged nucleotides. We then QC our PCR product on a gel to confirm expected length. 

Next we move to Sanger sequencing, which uses fluorescent tagged nucleotides at the terminating end. The sequences of each length  runs through a capillary, smallest first. The laser beam and excitation then indicates which nucleotide is which and these are then lined up to give you a sequence! 

![](https://www.sigmaaldrich.com/deepweb/assets/sigmaaldrich/marketing/global/images/technical-documents/protocols/genomics/sequencing/sanger-sequencing_steps_process_diagram/sanger-sequencing_steps_process_diagram.png)

![](https://letstalkscience.ca/sites/default/files/styles/width_800px/public/2020-07/capillary_gel_electrophoresis.png?itok=vlzUQLl7)

We then get a chromatogram to tell us what the sequence is! 

# Looking at sequences 

We have two files for each sequence. We are working through an example sample for this tutorial. 

If the region of interest is longer than the sequencer can do in one run, we have to have both a forward and reverse for each to match them together to make a longer sequence. 

Here is what this looks like in Geneious.   

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/sanger/geneious_snapshot.png?raw=true)

The colors show up as peaks on the chromatogram with the intensity of the fluorescence of that color indicated as the height of the peak.  

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/sanger/geneious_twoseqs.png?raw=true)

We are using the "sequence view" and have "allow editing" enabled. 

The sequencer has interpretted the chromatograms to make the "right" call. But we need to QC these calls and make sure that the base calls are correct. 

We get a lot of noise at the beginning of the sequence in particular. There will be an N where a base call could not be made. 

For example, we deleted the first 25 bp where there were many N's. 

Sometimes, we can get "dye blobs" that result from the capillary of the sequencer. In that case, we may need to change the calls. 

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/sanger/dyeblob_before.png?raw=true)

We changed the N's here to T's! 

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/sanger/dyeblob_after.png?raw=true)

We can also check this information with the reverse sequence. If you don't feel you can change it, you leave it as an N. 

We also removed the noisy end of the sequence.
 
![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/sanger/deleted_end.png?raw=true)

Here is our edited sequence now. 

![](https://github.com/AHuffmyer/ASH_Putnam_Lab_Notebook/blob/master/images/NotebookImages/sanger/revised_sequence.png?raw=true)