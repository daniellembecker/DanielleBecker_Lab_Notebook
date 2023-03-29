---
layout: post
title: Mcap 2021 Stable Isotope Metabolomics Analysis Part 1
date: '2023-03-28'
categories: Mcapitata_LarvalTemp_2021
tags: Mcapitata Metabolomics Multivariate R StableIsotopes
---

This post details the initial stable isotope metabolomics analysis for the *Montipora capitata* 2021 larval temperature exposure. In this analysis I conduct initial visualizations and multivariate analysis for total ion counts for all metabolites.  

# Overview 

In this project, I exposed *M. capitata* larvae to high (+ ~2°C) and ambient temperatures in Hawai'i in 2021. At the end of this exposure, I incubated larvae at the treatment temperature with isotope labels and collected samples for stable isotope metabolomics. 

During incubations we labeled samples (n=1 per tank, n=6 per temperature) with carbon-13, carbon-12 (isotope control), and dark carbon-13 (photosynthesis control). 

For more details on this project, [see my notebook posts for the Mcap2021 project](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/categoryview/#mcapitata-larvaltemp-2021).   

The GitHub repo for this project [is here](https://github.com/AHuffmyer/SymbioticIntegration).  

The script for analysis today [can be found here on GitHub](https://github.com/AHuffmyer/SymbioticIntegration/blob/main/Mcap2021/Scripts/metabolomics_total_counts.Rmd). 

# Data processing 

There are two datasets for this analysis - one with metabolite ion counts (abundance) and one with proportion labeling fraction (proportion of carbon with C13 label). In this initial analysis, I am first looking at total metabolite ion counts to view treatment differences in metabolomic profiles. In the next step, I will analyze fraction labelling patterns.  

During processing, Eric C. (Rutgers) confirmed that labeling of C13 was not present in C12 samples and was greatly reduced in dark C13. This shows that the data are reliable and incubations and extractions were done correctly. 

First, metabolite ion counts were normalized to the median ion count in each sample. This accounts for any variation in the number of larvae that were in each sample. The code looks like this:  

```
counts <- counts %>%  mutate(across(X12C_Ambient.A1_I7:Dark.13C_High.Pool_I40, ~ .x / median(.x, na.rm = TRUE)))
```

I manipulated the dataframe to obtain tank, sample, and treatment information from the sample name. I then summed all ion counts detected across all isotopic labels. For all samples (C12 or C13 labeled), the analysis results in ion counts for each carbon (C12 + C13-carbon 1 + C13-carbon 2 and so on). In order to obtain total ion counts for each metabolite in each sample, I sumed all ion counts for all detected isotopes. Note that there will be no ion counts in the C13 detection for the C12 samples.  

```
total_counts<-counts%>%  group_by(compound, sample, isotope, treatment, tank, tube)%>%  summarise(total.ion.count=sum(ion.count.norm))

```

This resulted in a data frame with ion counts normalized to sample median for each compound for each sample.  

# PCA visualization 

I next generated PCA visualizations and conducted PERMANOVA analyses.  

### PCA: Isotope effects  

First, I looked at differences in metabolomic profiles between sample type (C12, C13, dark C13).  

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/Mcap2021/pca_isotope.jpg)

```
Permutation test for adonis under reduced modelTerms added sequentially (first to last)Permutation: freeNumber of permutations: 999adonis2(formula = vegan ~ isotope, data = counts_pca, method = "eu")         Df SumOfSqs      R2      F Pr(>F)    isotope   2  1537.28 0.61887 20.297  0.001 ***Residual 25   946.72 0.38113                  Total    27  2484.00 1.00000      
```

There are clearly differences (p=0.001) between the C13 samples and other samples. This could be an artifact of the way the data was collected or could be an effect of C13 exposure under light conditions. 

As expected, the dark C13 samples are different than the light C13 samples. 

### PCA: C12 samples  
 
Next, I looked at differences between high and ambient treatments within the C12 and C13 samples separately.  

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/Mcap2021/pca_treatment_12C.jpg)  

```
Permutation test for adonis under reduced modelTerms added sequentially (first to last)Permutation: freeNumber of permutations: 999adonis2(formula = vegan ~ treatment, data = counts_pca_12C, method = "eu")          Df SumOfSqs      R2      F Pr(>F)  treatment  1   204.18 0.20175 2.5275  0.016 *Residual  10   807.82 0.79825                Total     11  1012.00 1.00000         
```

There is a significant treatment difference in metabolomic profiles in the 12C samples.  

### PCA: C13 samples  

Finally, I looked at separation in high and ambient treatments within the C13 samples. These are the samples that will be most important to our results, since we will look at these results for fraction labeling.  

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/Mcap2021/pca_treatment_13C.jpg)  

```
Permutation test for adonis under reduced modelTerms added sequentially (first to last)Permutation: freeNumber of permutations: 999adonis2(formula = vegan ~ treatment, data = counts_pca_13C, method = "eu")          Df SumOfSqs     R2      F Pr(>F)   treatment  1   230.13 0.2274 2.9433  0.004 **Residual  10   781.87 0.7726                 Total     11  1012.00 1.0000 

```

There are clear differences by treatment in the C13 samples - this is interesting!  

# PLSDA and VIP analysis 

Similar to the above approach, I conducted supervised PLS-DA analyses of metabolite profiles between treatments in the C13 and C12 samples separately. This will allow us to compare metabolite profiles between the isotope treatments and between temperature treatments.  

In this analysis, I will generate PLS-DA visualizations and extract variables of importance of projection (VIP), which are the metabolites that drive the differences between treatments.  

## C13 samples  

I generated a PLS-DA by temperature treatment within the C13 samples.  

```
#assigning datasets X_C13 <- metabolites%>%  dplyr::select(compound, isotope, treatment, tank, log.ion.counts)%>%  filter(isotope=="13C")%>%  spread(compound, log.ion.counts) #generate spread data framelevels(as.factor(X_C13$treatment))Y_C13 <- as.factor(X_C13$treatment) #select treatment namesY_C13X_C13<-X_C13[4:95] #pull only data columns# run PLSDA MyResult.plsda_C13 <- plsda(X_C13,Y_C13) # 1 Run the methodplotIndiv(MyResult.plsda_C13)    # 2 Plot the samplesplotVar(MyResult.plsda_C13, cutoff = 0.7)    treatment_cols<-c("blue", "red")             pdf("Mcap2021/Figures/Metabolomics/PLSDA_C13.pdf", width=9, height=6)plotIndiv(MyResult.plsda_C13, col=treatment_cols, ind.names = FALSE, legend=TRUE, legend.title = "Treatment - C13 Samples", ellipse = FALSE, title="", style = "graphics", centroid=FALSE, point.lwd = 2, cex=2)dev.off() 
```

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/Mcap2021/plsda_c13.png)

As seen in the PCA, there is separation between high and ambient treatments.  

I then extracted the VIP's and filtered to VIP>1 to look at the metabolites driving this change.  

```
#extracttreatment_VIP <- PLSDA.VIP(MyResult.plsda_C13)treatment_VIP_df <- as.data.frame(treatment_VIP[["tab"]])treatment_VIP_df# Converting row names to columntreatment_VIP_table <- rownames_to_column(treatment_VIP_df, var = "Metabolite")#filter for VIP > 1treatment_VIP_1 <- treatment_VIP_table %>%   filter(VIP >= 1)write_csv(treatment_VIP_1, "Mcap2021/Output/Metabolomics/Treatment_C13_VIPs.csv")
```

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/Mcap2021/vip_c13.png)

I also viewed these VIP's by abundance and treatment. 

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/Mcap2021/VIP_Abundance_C13.jpg)

To see clearly which metabolites were higher in the High treatment, I calculated relative change between ambient and high treatments summarized at the treatment level. 

```
C13_vip_diff_plot<-C13_vip_metabolites%>%  arrange(compound)%>%  group_by(compound, treatment)%>%  summarise(mean=mean(log.ion.counts))%>%  group_by(compound)%>%  summarise(change=(mean-first(mean))/first(mean))%>% #calculate the difference of high from ambient  filter(change !=0) #remove zero values (comparing ambient to ambient)
```

Negative values indicate metabolite abundance high is lower than ambient, positive indicates high treatment is higher than ambient. 

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/Mcap2021/VIP_Difference_C13.jpg)

NADP+, NAD+, arginine-glutamine, phosphocholine, inositol, asparagine, and aspartate are higher at high temperature and the others are lower at high temperature.  

I submitted the VIP list to metaboanalyst to conduct KEGG/class enrichment. These were the significantly enriched pathways in this VIP list: KEGG Pathways: 
- Aminoacyl-tRNA biosynthesis (FDR p<0.001)- Arginine biosynthesis (FDR p<0.001)- Alanine, aspartate, and glutamate metabolism (FDR p<0.001)- Valine, leucine, and isoleucine biosynthesis (FDR p=0.004)- Nicotinate and nicotinamide metabolism (FDR p=0.026)Class Enrichment: 
- Amino acids and peptides (FDR p<0.001)- Monosaccharides (FDR p<0.001)- Cholines (FDR p=0.033)- Sulfonic acids (FDR p=0.033)- Short-chain acids and derivatives (FDR p=0.047)

## C12 samples  

I then generated a PLS-DA for the C12 samples.  

```
#assigning datasets X_C12 <- metabolites%>%  dplyr::select(compound, isotope, treatment, tank, log.ion.counts)%>%  filter(isotope=="12C")%>%  spread(compound, log.ion.counts) #generate spread data framelevels(as.factor(X_C12$treatment))Y_C12 <- as.factor(X_C12$treatment) #select treatment namesY_C12X_C12<-X_C12[4:95] #pull only data columns# run PLSDA MyResult.plsda_C12 <- plsda(X_C12,Y_C12) # 1 Run the methodplotIndiv(MyResult.plsda_C12)    # 2 Plot the samplesplotVar(MyResult.plsda_C12, cutoff = 0.7)                pdf("Mcap2021/Figures/Metabolomics/PLSDA_C12.pdf", width=9, height=6)plotIndiv(MyResult.plsda_C12, col=treatment_cols, ind.names = FALSE, legend=TRUE, legend.title = "Treatment - C12 Samples", ellipse = FALSE, title="", style = "graphics", centroid=FALSE, point.lwd = 2, cex=2)dev.off() 
```
![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/Mcap2021/plsda_c12.png)

There is also separation between metabolomic profiles of high and ambient larvae in the C12 labeled samples. It is confirmatory that the separation we see in C13 is also seen in the isotope control samples. 

I then extracted VIPs for the C12 samples. 

```
#extracttreatment_VIP_C12 <- PLSDA.VIP(MyResult.plsda_C12)treatment_VIP_C12_df <- as.data.frame(treatment_VIP_C12[["tab"]])treatment_VIP_C12_df# Converting row names to columntreatment_VIP_C12_table <- rownames_to_column(treatment_VIP_C12_df, var = "Metabolite")#filter for VIP > 1treatment_VIP_C12_1 <- treatment_VIP_C12_table %>%   filter(VIP >= 1)write_csv(treatment_VIP_C12_1, "Mcap2021/Output/Metabolomics/Treatment_C12_VIPs.csv")
```

These are the >1 VIPs for the C12 samples driving separation between ambient and high. 

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/Mcap2021/vip_c12.png)

I also viewed these VIP's by abundance and treatment. 

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/Mcap2021/VIP_Abundance_C12.jpg)

Finally, I calculated relative change in abundance using the approach described above. 

```
C12_vip_diff_plot<-C12_vip_metabolites%>%  arrange(compound)%>%  group_by(compound, treatment)%>%  summarise(mean=mean(log.ion.counts))%>%  group_by(compound)%>%  summarise(change=(mean-first(mean))/first(mean))%>% #calculate the difference of high from ambient  filter(change !=0)
```
  
![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/Mcap2021/VIP_Difference_C12.jpg)

UDP-D-Glucose, UDP-N-acetyl-glucosamine, NADP+, NAD+, Phenyllactic acid, glycine, CDP-choline, arginine-glutamine, AMP, inositol, adenosine, adenine, and a-ketoglutarate are increased at high temperature. There is overlap in this list from the C13 samples, but there are additional metabolites that show up in the C12 samples. 

I submitted the VIP list to metaboanalyst to conduct KEGG/class enrichment. These were the significantly enriched pathways in this VIP list: KEGG Pathways:  
- Aminoacyl-tRNA biosynthesis (FDR p<0.001)  - Arginine biosynthesis (FDR p<0.001)  - Alanine, aspartate, and glutamate metabolism (FDR p=0.003)  - Valine, leucine, and isoleucine biosynthesis (FDR p=0.006)  - Purine metabolism (FDR p=0.017)  
*In this set purine metabolism is significant whereas in C13 nicotinamide metabolism is significant. All others are the same.* Class Enrichment: 
- Amino acids and peptides (FDR p<0.001)- Pyrimidines (FDR p<0.001)- Purines (FDR p<0.001)- TCA Acids (FDR p<0.001)
*In C13 samples, Pyrimidines, purines, and TCA acids are not enriched, but cholines, sulfonic acids, monosaccharides, and short chain acids are*  

# Thoughts and relevant literature 

Overall, the key pathways separating high and ambient treatments deal with nitrogen and peptide/dipeptide metabolism. Here are some quick thoughts and references from the literature that may help generate some hypotheses for these results.  

*Radecker et al. 2020:* [paper link](https://www.pnas.org/doi/pdf/10.1073/pnas.2022653118)

- This study found lower expression of glutamate synthase (GS) and lower expression of glutamine oxoglutarate aminotransferase (GOGAT) under heat stress. Interestinly, they also found lower anabolic activity of glutamine dehydrogenase (GDH) and higher catabolic activity under stress. Glutamine synthase incorporates ammonium into glutamate to generate glutamine. The GOGAT pathway takes glutamine and alpha-ketoglutarate and generates glutamate. GDH converts glutamate to or from alpha-ketoglutarate by the addition or removal of ammonium. These three pathways are responsible for nitrogen cycling and availability in the host.  

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/Mcap2021/radecker_2020.png)  
Radecker et al. 2020

- There is evidence in our data that shifts in these pathways of nitrogen cycling between ambient and high. In the C12 dataset, alpha-ketoglutarate is higher slightly in high temperature larvae with a decrease in glutamine. This would suggest an increase in GOGAT activity at high temperature combined with increased rate of the TCA cycle. Glutamine was also reduced at high temperature in the C13 samples. We would not be able to pick up glutamate abundance because we are looking at only negative polarity metabolites. 
- The authors of this study suggest that this shift can result from nutrient deficit and starvation that leasds to increased ammonium availability in the host. 
- We also saw increased photosynthesis in high temperature larvae. It is possible this could relate to increase ammonium availability in the holobiont. 

*Chiles et al. 2022*: [paper link](https://www.frontiersin.org/articles/10.3389/fmars.2022.1035523/full)   

- The abundance of dipeptides increased with thermal stress - there was an accumulation of dipeptides (Williams et al. 2021). 
- In our data, we saw an increase in arginine-glutamine at high temperature, but a decrease in glutamine, suggesting incorporation of glutamine into dipeptides. 
- Glutamate is a precursor for lysine, valine, alanine, and arginine. Alanine, lysine, and valine are all higher in ambient conditions compared to high, suggesting lower availability of glutamate and/or lower production of amino acids from glutamate at high temperature. 
- Glutamate can also be used in the urea cycle to make citrulline, which is higher at ambient temperature, further suggesting lower availability of glutamate and/or urea cycle production. 
- In Chiles et al. 2022, arginine-glutamine had high nitrogen labeling and therefore high turnver of dipeptide pools. It is useful to compare labeling of the individual free amino acids to the dipeptide to determine how it was synthesized. 
- Turnover rates are determiend by production rate and pool size. 
- Dipeptides could also come from the symbiont and are translocated to the host. From Chiles et al. 2022: "It is possible that rapid synthesis of amino acids into dipeptides allows Symbiodiniaceae to diminish bioavailability of nitrogen to reduce the potential influence of host nitrogen mediated control." 

![](https://www.frontiersin.org/files/Articles/1035523/fmars-09-1035523-HTML/image_m/fmars-09-1035523-g006.jpg)
Chiles et al. 2022

*Williams et al. 2021:* [paper link](https://www.science.org/doi/10.1126/sciadv.abd4210)  

- There was differential accumulation of four dipeptides in symbiotic vs asymbiotic Aiptasia, which may suggest a role in symbiosis and heat stress 
- Corals divert a huge amount of energy to producing Montiporic acids (MAs). We saw high abundance of MA's in the ion counts in our data, but they were not contributing to differences in treatments, suggesting they are maintained at high levels. 
- "These compounds are disubstituted acetylenes with a carboxyl group linked to two alkyne (carbon-carbon triple bond) groups, followed by an unbranched alkane tail. The four known MAs (MA-A to MA-D; Fig. 1D) have antimicrobial activity, are cytotoxic against leukemia cells, and reduce the photosynthetic competency of coral symbionts." - Williams et al. 2021
- This study found decreases in glucose in heat stress, but we did not find that in our study. This could be due to the short term heat stress and the lower intensity of heat stress relative to previous work. This could also represent metabolic resilience to stress. 
- Arginine-glutamine (RQ) was acculumulated under heat stress in this study, as seen in ours. This could be from increased proteolysis and/or insufficient peptide clearance. The majority of metabolites accumulated in thermal stress are dipeptides, which also matches our data.  
- This study also found accumulation of lysine-glutamine, arginine-valine, and arginine-alanine. 
- Dipeptides can also be small molecule regulators, such as H+ buffers, antioxidants, and glucose regulation. Specifically, RQ can ameliorate impacts of oxygen imbalance and oxygen damage in lungs. 
- This study also found increases of methionine, methionine sulfoxide, CDP-choline. We also found increases in CDP-choline at high temperature (C12 samples). 
- CDP Choline is "a key metabolite involved in phospholipid metabolism, cell signaling, and glutamate transport and an intermediate in betaine synthesis". 
- Our work suggests RQ is host derived. In this study, the authors say, "We postulate that RQ may be a host response to oxygen stress resulting from redox imbalance in the coral or the algal symbionts caused by the high-temperature treatment."

*Hillyer et al. 2017a:* [paper link](https://link.springer.com/article/10.1007/s11306-017-1306-8)

- When we check for dark samples, we will look for something like the following statement, "No significant enrichment of either symbiont, or host pools was detected in the dark incubations, suggesting that label incorporation was via symbiont photosynthesis and translocation of organic products."
- Found significant enrichment in glucose and glycerol pools after 6 and 9 days of stress. We may not see this due to shorter term stress. "As per the findings of Hillyer et al. (2017), the degree of labelling of both symbiont and host glucose pools observed here suggests that despite photodamage, de novo biosynthesis and translocation of the compound from remaining in hospite symbionts was on-going and may actually increase per symbiont under thermal stress." No difference in glucose in our study may support this.
- Increased enrichment of galactose indicates rapid incorporation and turnover of the galactose pool and biosynthesis from a precursor such as glucose, translocation from symbiont, or decline in downstream metabolism. 
- Found increases in abundance of inositol after 6 and 9 days of stress. No label detected in inositol suggests this is a host-mediated response to stress. "The accumulation of compatible solutes observed in the intracellular pools of the heat stressed host, likely indicate exposure to osmotic stress and mechanisms to maintain internal osmolality. In addition to their functions as compatible solutes, the roles of inositol derivatives in the cnidarian-dinoflagellate symbiosis are currently receiving increased attention (Bertucci et al. 2013; Rosic et al. 2015), in part because of their highly conserved functions in cellular recognition, signalling and cell development (Boehning et al. 2005; Chakraborty et al. 2008)."
- [Hillyer et al. 2017b](https://link.springer.com/article/10.1007/s11306-017-1306-8) will also be a helpful paper for labeling data. 

*Pei et al. 2022:* [paper link](https://link.springer.com/article/10.1007/s00216-022-04294-y)  

- In thermal stressed P. decussata, there was a decrease in lipids and amino acids, indicating increased lipid hydrolysis and aminolysis contributing to increased gluconeogenesis to meet energy demand. 
- They found that trypophan, panthenol, pantothenate were upregulated in stress and others such as glutamine, glutamate, aspartate, inosine were downregulated. We saw down regulation of glutamine and aspartate in our samples in heat stress. 
- Similar to our results, arginine, animoacyl, alanine/aspartate, and other peptide metabolism was enriched in thermal stress. 
- Also saw accumulation of peptides under stress. 

*What about accumulation of NAD+ and NADP+?* 

- In high temperature there was a large accumulation of NAD+ and NADP+ and enrichment of nicotinate and nicotinamide metabolism in the C13 samples. 
- [Xie et al. 2020](https://www.nature.com/articles/s41392-020-00311-7) provides a review of NAD+ metabolism. Here is a figure from the paper:  

![](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41392-020-00311-7/MediaObjects/41392_2020_311_Fig1_HTML.png?as=webp)

- Accumulation of NAD+ and NADP+ in our samples could result from antioxidant defence, anabolic reations, and signaling that uses NADPH and generates NADP+. If the Pentose phosphate pathway activity is reduced, then the pool of NADP+ would increase, as we see in our samples. We should look for metabolite evidence of the pentose phosphate pathway to see if we can infer rates of this pathway. 
- Further, electron transport chain, lactic adic generation, and PUFA desaturation would use NADH and generate NAD+. If glycolysis, fatty acid oxidation, and the TCA cycle are slower or do not match the rate of accumulation, then the pool of NAD+ would increase. We can also look for metabolite information on rates of these pathways in our samples. 
- Further, NAD+ is degraded through DNA repair, ER stress, RNA processing, metabolism, signaling, etc. If the pool is increased, then either the production of NAD+ increases, the degradation pathways decrease, or both are happening at the same time. 

# Next steps 

I need to address a few questions with the Rutgers team: 

- Is summing ion counts across all isotope label detections the correct way to get to total metabolite abundance? If not I will need to revise these calculations. 
- How do we want to address differences in C12 and C13? 
- Are there treatment differences in C13 labeling in the dark C13 samples and how do we need to correct for this? 

As we test these hypotheses we can see if adding gene expression would be helpful (for example, to test GOGAT, GDH, GS enzyme expression). 
We could also add positive polarity analysis using the duplicate extract if required to test hypotheses based on these data. 

Next I will explore the label fraction datasets.  