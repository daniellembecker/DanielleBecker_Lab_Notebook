---
layout: post
title: F1F2 project data analysis Hawaii 2023
date: '2023-07-31'
categories: Mcap_F1F2_2023
tags: Mcapitata Photophysiology R Respirometry
---

This post details respirometry data analysis of the *Montipora capitata* F1F2 project at the Hawaii Institute of Marine Biology (HIMB) in June 2023. 

# Overview 

We collaborated with the [Coral Resilience Lab](https://www.coralresiliencelab.com/) at HIMB to measure metabolic rates in an ongoing project to examine the effects of selective breeding on coral thermal tolerance. In this project ("F1F2") we examined metabolic rates of coral larvae from F1 and F2 generations of wildtype corals (fertilized from gametes from the wild population in Kaneohe Bay) and selectively bred corals with non-bleaching history (fertilized from parents that did not bleach in previous bleaching events). This post details the respirometry data collection, analysis, and findings from these groups of corals. 

See the [Coral Resilience Lab website](https://www.coralresiliencelab.com/) for more information on coral selective breeding and resilience research. 

# Respirometry Protocol

We examined metabolic rates using respirometry protocols used in our examination of larvae with different symbiont species in June 2023. The full notebook post from our work in this project can [be found on my notebook website here](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/Hawaii-2023-Coral-Spawning-and-Field-Expedition-Daily-Entry/) and at the [GitHub repository for our project here](https://github.com/AHuffmyer/larval_symbiont_TPC). 

This protocol is adapted from the [Putnam Lab SDR equipment and respirometry protocol](https://github.com/Putnam-Lab/Lab_Management/blob/master/Lab_Resources/Equipment_Protocols/Respirometry_Protocol/SDR-Respirometry-Protocol.md). 

The Coral Resilience lab successfully collected larvae from parent colonies (F1 WT) of non-bleached corals and wildtype adults (F1 WT) as well as larvae spawned from 5-year-old WT (F2 WT) and NB offspring (F2 NB). 

Respirometry measurements were conducted with a dual Sensor Dish Reader (SDR, PreSens) in combination with a glass 24-well microplate (80 µL wells) in incubators with an AquaIllumination Prime 16HD light source. Oxygen concentration (µmol per litre) is read every 15 seconds in each well. Pools of 5 larvae of each type (F1/F2 WT/NB) are loaded into each well. Each plate contained 5 wells of each larval group along with 4 blanks. All runs were run on duplicate plates (n=10 per larval group per temperature). Blanks and larval wells were loaded with 1 µm filtered seawater. 

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/Hawaii2023/sdr_picurves.jpeg)

We ran photosynthesis (P) and respiration (R) measurements across a temperature x light profile. We first dark adapted larvae in incubators and then exposed them to 0 PAR (dark respiration), 50 PAR (low light photosynthesis), 100 PAR (medium light photosynthesis), and 500 PAR (high light photosynthesis), followed by another measurement at 0 PAR for light-enhanced respiration. These light values were chosen based on a PI curve we conducted in *M. capitata* larvae that can be found on [GitHub here](https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/figures/larval_pi_curves/nls_curves_all.png). Notice saturation of the curve at approximately 100-200 PAR. We therefore selected a light value that was lower, at, and above this saturating value. 

We did one set of measurements at ambient (27°C), high (33°C), and extreme temperature (36°C). All plates were calibrated with 100% air saturated FSW before each run at the respective temperature. This provided us with the ability to calculate dark respiration (Rd), light enhanced respiration (LER), and photosynthesis across light levels. From these measurements, we can calculate gross photosynthesis, net photosynthesis, Rd, LER, and P:R ratio at low, medium, and high light across temperatures. 

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/Hawaii2023/sdr2.jpeg) 

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/Hawaii2023/sdr.jpeg) 

Here is an example of the oxygen data from this temp x light protocol. You can clearly see the respiration and light phases in the order of dark - low light - medium light - high light - and dark. 

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/Hawaii2023/sdr3.jpeg) 

Data and metadata can be found [on GitHub here](https://github.com/AHuffmyer/larval_symbiont_TPC/tree/main/data/F1F2_sdr). 

Data analysis and figures can be found [on GitHub here](https://github.com/AHuffmyer/larval_symbiont_TPC/tree/main/figures/F1F2_sdr) and [scripts can be found here for rate extraction](https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/F1F2_sdr_extraction.Rmd) and [rate analysis here](https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/F1F2_sdr_analysis.Rmd).  

Steps in data analysis are described below. 

# Data Analysis 

Oxygen concentration during each light phase of each run was extracted using the [LoLinR package in R](https://github.com/colin-olito/LoLinR) using the percentile rank method and alpha set to 0.4. Metabolic rates were corrected by subtracting blank values averaged for each replicate plate and then normalized to well volume and the number of larvae to express rates as nmol O2 per larva per min. 

Here is an example of the code loop used to extract rates from each replicate in each run. This code reads in the oxygen data from the data file, splits the data into each light phase by a time stamp data frame, and then extracts a slope that is calculated as the metabolic rate. 

The script for slope extraction can be found on [GitHub here](https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/F1F2_sdr_extraction.Rmd).  

```
for(file in 1:length(file.names)) { # for every file in list start at the first and run this following function
```

This script produces a PDF file that shows the statistics and slope extraction. Here is an example of respiration measurements in the dark: 

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/F1F2/r_extract.png) 

And here is an example of measurements from a light phase (photosynthesis): 

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/F1F2/p_extract.png) 

You can see the negative slope in the respiration data indicating oxygen consumption and a positive slope in the photosynthesis data indicating oxygen production. 

After slope extraction, we then correct the values by substracting blank values and accounting for volume and number of larvae. 

Finally, we convert oxygen consuption values to the inverse such that data is expressed on a positive scale for respiration. In cases where respiration oxygen consumption values were <0 (indicating net oxygen production), value was set to 0 (indicating no respiration). In cases where photosynthesis oxygen production values were <0 (indicating net oxygen consumption), value was set to 0 (indicating no photosynthesis).

We are now ready for data analysis! 

# Findings 

The script for data analysis can be found on [GitHub here](https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/F1F2_sdr_analysis.Rmd). 

## Metabolic rates across light and temperature

I first plotted metabolic rates across the light curve for each temperature. The plot below shows oxygen values across light levels. Note that the values below 0 are oxygen consumption (respiration) and values above 0 are oxygen production (photosynthesis). 

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/F1F2/pr_dots.png) 

There are a few things to take away from this first look at the data: 

- The shape of the curve changes by temperature, with the curve becoming more shallow (less photosynthesis) at higher temperatures. 
- The shape of this curve is expected with higher photosynthesis with higher light values. 
- Respiration values are greater at the ending dark (0.2 PAR) than the beginning dark stage (0.1 PAR). *Note that 0.1 and 0.2 are truly 0 PAR, but I added the decimal to order the X axis in the correct order.* This is expected because light enhanced respiration is greater due to the stimulation of animal respiration from the production of nutritional products and oxygen from symbiont photosynthesis. This stimulation is not present in dark adapted larvae since photosynthesis was not active prior to respiration measurements. 
- There is not a large obvious difference in larval group in the shape of these curves, but we may detect differences by calculating metrics which we will do next. 

I tested for light x temperature x group interactions with an ANOVA. 

```
summary(aov(P.nmol.org.min~PAR*Group*Temperature, data=pr_data))

                       Df  Sum Sq Mean Sq F value   Pr(>F)    
```

Metabolic rates are affected by PAR, PAR x group interaction, and a PAR x temperature interaction. This indicates that larval groups have different P and/or R rates and that P and/or R rates are different by temperature. 

There is no interaction between group x temperature or group x PAR x temperature, indicating that differences in larval groups are driven by light, rather than strongly different by temperature. 

## Calculating metabolic metrics 

I then calculated metabolic rates to generate the following metrics: 

1. **Light enhanced respiration**: Respiration rate in the dark following exposure to light
2. **Dark respiration**: Respiration in the dark following dark adapt period 
3. **Net photosynthesis (low, medium, high light)**: Oxygen production in excess of host respiration
4. **Gross photosynthesis (low, medium, high light)**: Total oxygen production (P gross + R) 
5. **P:R ratio (low, medium, high light)**: Photosynthesis divided by light enhanced respiration. We use light enhanced respiration for this calculation because respiration will be light enhanced when photosynthesis is active. This is the most ecologically relevant metric to compare oxygen demand and production. A P:R of 1 indicates that P meets R demand. A P:R greater than 1 indicates oxygen production that is higher than demand. 

These metrics are calculated using the following script: 

```
pr_data_calc<-pr_data%>%
```

Note that for P gross, P net, and P:R ratio we will have values at 50 PAR (low light), 100 PAR (medium light), and 500 PAR (high light). 

### LER and RD metrics  

Here is a plot of LER and Rd. 

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/F1F2/ler_rd.png) 

I ran an ANOVA to examine the effect of temperature and group on LER and RD. 

Here are the results for LER: 

```
model<-pr_data_calc%>%

                   Df   Sum Sq   Mean Sq F value Pr(>F)   

$Group
```

Light enhanced respiration is different by larval group. From the plot we can see that F2 NB had the lowest LER. LER in F2 NB was lower than F1 NB and F2 WT. There is no effect of temperature on LER. This could indicate that F2 NB corals have lower energy demand to survive in the same conditions as other groups. We will look at other metrics to interpret this difference. 

Next I analyzed Rd (dark respiration): 

```
model<-pr_data_calc%>%

                   Df   Sum Sq   Mean Sq F value   Pr(>F)    

$Temperature

```
Dark respiration decreased with temperature where Rd was lower at 33 and 36 than at 27. This could indicate metabolic depression in larvae at elevated temperatures. This would make sense with the high temperatures we used for these measurements. 

It is interesting that LER did not have the same decrease across temperature. This is an indication that stimulation of respiration by photosynthesis counteracts metabolic depression that happens at high temperature measured by Rd. There was no difference in dark respiration by larval group. 

### P net and P gross 

I next looked at P net and P gross across light and temperature. 

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/F1F2/pr_metrics.png) 

I used an ANOVA to analyze the effect of light ("Metric"), temperature, and group for P gross. 
```
model<-pr_data_calc%>%

                          Df  Sum Sq  Mean Sq F value   Pr(>F)    
```

There was a significant effect of light, temperature, group, light x temperature, and group x temperature on P gross - lots of interesting things are happening! 

- P gross is higher as light increases. 
- P gross decreases as temperature increases, but the effect of temperature is more pronouced at high light. At high light, P gross is more sensitive to temperature than at low light. 
- Photosynthesis is lower in the F2 NB than the F1 NB and the F2 WT. This is interesting to note. 

The lower P gross in the F2 NB aligns with what we saw in LER. LER was also lower in this group and this is likely due to decreased photosynthetic output and therefore less stimulation of respiration as compared to other groups. This may come with interesting metabolic trade offs between photosynthetic production and thermal tolerance. 

I also analyzed P net. 

```
model<-pr_data_calc%>%

                          Df  Sum Sq   Mean Sq F value   Pr(>F)    

```

There were significant effects of light and temperature on P net, but no effect of group. This shows that differences in photosynthesis are driven by respiration of the host because differences are seen in P gross but not P net. 

### P:R ratio 

Finally, I looked at P:R ratios. Note that a P:R of 1 indicates that LER equals gross photosynthesis. A value of >1 indicates gross photosynthesis produces more oxygen than consumed through LER.  

![](https://raw.githubusercontent.com/AHuffmyer/ASH_Putnam_Lab_Notebook/master/images/NotebookImages/F1F2/pr_ratio.png)

I used an ANOVA to test for effects of light ("Metric"), temperature, and group on P:R ratios. 

```
model<-pr_data_calc%>%

                          Df Sum Sq Mean Sq F value   Pr(>F)    

```

There is an effect of light and temperature on P:R ratios, but no effect of group. 

- P:R ratios increase as light increases. This makes sense as we know photosynthesis increases with light. P:R ratios are highest at high light in ambient conditions.
- P:R ratios are more sensitive to temperature effects at medium to high light than low light. Especially at high light - there is greater P:R ratio at ambient with a steep decline at 33 and 36°C. This also tracks what we saw with gross photosynthesis declining at elevated temperatures especially under high light. 
- Because LER was unaffected by temperature, changes in P:R ratio are driven by changes in photosynthetic output. 

# Take Homes 

There are several interesting take homes to explore further with this data. 

1. Light enhanced resipration was lower in the F2 NB group and P gross was lower in the F2 NB group. This demonstrates that metabolic rates are strongly driven by lower symbiont productivity in the F2 NB group as compared to other groups. This suggests that there may be metabolic trade offs between photosynthetic activity and thermal tolerance in these larvae. 
2. F2 NB corals maintained the same P:R ratios as other groups, indicating that even though symbiont productivity was lower, they are maintaining the same energy balance as other groups. This would suggest that the F2 NB experience energetic savings under ambient and stress conditions. 
3. These metrics were lower in F2 NB corals than the same generation of WT (F2 WT) and the previous generation of NB (F1 NB). This suggests that there is a generational effect of bleaching history on metabolic rates and may have interesting implications for selective breeding. 

In order to fully understand the implications of these differences in metabolic rates, we need to consider phenotypic and performance information. Do the F2 NB larvae have higher or lower survivorship than other groups? Do they experience differences in growth or energy content? Do they have greater thermal tolerance? This information is critical to determine whether the differences in metabolic rates outlined here are due to energetic savings and advantageous metabolic strategies or are a signal of metabolic dysfunction.   

# Next Steps 

The next step in this project is to determine what data is of interest to include in the project manuscript and to interpret in the context of phenotypic or performance metrics. 
