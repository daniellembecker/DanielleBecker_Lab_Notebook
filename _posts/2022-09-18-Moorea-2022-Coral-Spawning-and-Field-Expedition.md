---
layout: post
title: Moorea 2022 Coral Spawning and Field Expedition
date: '2022-09-18'
categories: Acropora_Integration_Moorea_2022
tags: Acropora Husbandry Reproduction ExperimentalDesign
---
This post includes the full notebook for Sept-Oct 2022 *Acropora* coral spawning and field expedition in Moorea, French Polynesia. These notebook posts by A. Huffmyer are updated in order of most recent to oldest. 

# Project Overview  

The goal of this field expedition for A. Huffmyer is to conduct field work to invesigate metabolic integration of coral host and symbiont in early life history as part of funded National Science Foundation Ocean Sciences Postdoctoral Fellowship and the University of Washington eScience Data Science Postdoctoral Fellowship. More information on the abstract of this award can be [found here](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2205966&HistoricalAwards=false). A visual abstract of the project is below.  

![NSF Overview](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/NSF_OCE/NSF_overview.png)  

The goal for this expedition is to answer all or some of the following questions dependent on material, spawning success, and timing:   
1. How does thermal stress impact the timing and content of symbiotic nutritional exchange in early life history of *Acropora* sp.?   
2. How does parental exposure to elevated temperature impact the initiation of symbiotic nutritional exchange in offspring in *Acropora* sp.?   
3. How does exposure to elevated temperature impact metabolic exchange in *Acropora* sp. prior to and during reproduction?   

We will also continue in situ monitoring of *Acropora* and *Pocillopora* colonies on the reef in Moorea (Putnam Lab and E5 projects). 

*Field expedition members:*     
Ariana Huffmyer  
Pierrick Harnay  
Danielle Becker-Polinski  

# 18 September 2022  

### *Setting up temperature control*  

Ariana added heaters to the cooler systems in tanks A and C for high temperature control and connected to the Apex Neptune system. Temperatures programmed to be +4°C to see if the system can give us a high degree of elevated temperature without electrical limitations.  

![heaters](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/Moorea2022/heaters.jpeg)  

![heaters2](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/Moorea2022/heaters2.jpeg)  

![heaters2](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/Moorea2022/tanks2.jpeg)  

The script in the Apex program follows this temperature profile for each virtual stage.  

29C - 30.4C (with 1.5C daily fluctuation)   
Stage 1 = 29  
Stage 2 = 29.2  
Stage 3 = 29.4  
Stage 4 = 29.6  
Stage 5 = 29.8  
Stage 6 = 30.0  
Stage 7 = 30.2  
Stage 8 = 30.4    

This profile will be revised based on the diel fluctuation and the ambient temperature in Moorea. This profile is a test. The peak temperature occurs between 1500 and 1715 daily with this profile.  

The script in the Apex system is as follows:  

```
All master outlets 
Set OFF 
If Outlet High_Stage1 = ON Then ON 
If Outlet High_Stage2 = ON Then ON 
If Outlet High_Stage3 = ON Then ON 
If Outlet High_Stage4 = ON Then ON 
If Outlet High_Stage5 = ON Then ON 
If Outlet High_Stage6 = ON Then ON 
If Outlet High_Stage7 = ON Then ON 
If Outlet High_Stage8 = ON Then ON 

High_Stage1
Set OFF 
If High_1 < 29 Then ON 
If High_1 > 29 Then OFF 
If Time 01:00 to 04:30 Then OFF 
If Time 06:45 to 00:59 Then OFF

High_Stage2
Set OFF 
If High_1 < 29.2 Then ON 
If High_1 > 29.2 Then OFF 
If Time 01:00 to 03:00 Then OFF 
If Time 04:30 to 06:45 Then OFF
If Time 08:15 to 00:59 Then OFF 

High_Stage3
Set OFF 
If High_1 < 29.4 Then ON 
If High_1 > 29.4 Then OFF 
If Time 01:00 to 02:00 Then OFF 
If Time 03:00 to 08:15 Then OFF 
If Time 09:45 to 00:59 Then OFF 

High_Stage4
Set OFF 
If High_1 < 29.6 Then ON 
If High_1 > 29.6 Then OFF 
If Time 02:00 to 09:45 Then OFF 
If Time 11:15 to 00:30 Then OFF 

High_Stage5
Set OFF 
If High_1 < 29.8 Then ON 
If High_1 > 29.8 Then OFF 
If Time 01:00 to 11:15 Then OFF 
If Time 12:15 to 23:00 Then OFF 
If Time 00:30 to 00:59 Then OFF 

High_Stage6
Set OFF 
If High_1 < 30 Then ON 
If High_1 > 30 Then OFF 
If Time 01:00 to 12:15 Then OFF 
If Time 13:30 to 19:45 Then OFF 
If Time 23:00 to 00:59 Then OFF

High_Stage7
Set OFF 
If High_1 < 30.2 Then ON
If High_1 > 30.2 Then OFF 
If Time 01:00 to 13:30 Then OFF 
If Time 15:00 to 17:15 Then OFF 
If Time 19:45 to 00:59 Then OFF 

High_Stage8
Set OFF 
If High_1 < 30.4 Then ON 
If High_1 > 30.4 Then OFF 
If Time 01:00 to 15:00 Then OFF 
If Time 17:15 to 00:59 Then OFF

```

The temperature treatment was initiated at 1530 today. By 1800 the temperatures were at the desired level. We will let the profiles run for 24 h to view performance. Tomorrow we will look at the Apex data and download loggers.  

### *Setting up SDR and lights*  

Ariana set up the incubators, SDR, and lights that will be required for respirometry and stable isotope incubations.  

![incubators](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/Moorea2022/incubators.jpeg) 

Aqua Illumination Prime 16HD reef lights are used to provide light in the incubators for photosynthesis measurements and stable isotope incubations. 

Light settings were set to 40% intensity on all channels (UV, violet, royal, blue, green, deep red, cool white; no moonlight) positioned 25 cm above surface of the bottom of the incubator and 22 cm above the surface of the SDR plate. The goal is to obtain light at approx. 500 PAR at the level of the SDR plate.  

Light was measured at 6 positions at the level of the SDR plate. This measurement should be conducted at the start of every day of measurement.  

Incubators were labeled Incubator 1 and Incubator 2. Associated SDR 1 and 2 will be placed in the same number incubator. 

Light measurements were recorded in Google Drive datasheet. Average light intensity was 494 PAR in Incubator 1 and 487 PAR in Incubator 2.   

A HOBO Onset Pendant logger is placed in each incubator that can be launched at each respiration run to record temperature every 1 min. 

Logger 21002976 is in Incubator 1 and logger 2133598 is in Incubator 2.  

![resp1](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/Moorea2022/resp1.jpeg) 

![resp2](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/Moorea2022/resp2.jpeg) 

### *Calibrating SDR*  

We replaced oxygen sensor spots on the SDR microplates and need to calibrate since there are new spots from different batches than the old spots. This will also correct any drift of the spots.  

SDR systems were assembled and set up and connected to the PreSens SDR software. To calibrate we followed the following steps:  

1. Mixed solution of 100% oxygen seawater by bubbling filtered seawater and made 0% oxygen seawater by adding sodium sulfite.  
2. Set oxygen unit to % O2 and record every 3 minutes. Select batch PSt5-1624.  
3. Select log measurement to create data for calibration labeled with Date and "calibration". Saved on PC desktop in Moorea 2022 spawning folder. Will be uploaded to Google Drive and GitHub.  
4. Load half of wells with 0% solution and half with 100% oxygen solution in each plate. Record locations. Alternating wells with 0 and 100% water.  
5. Start measurements and let record for 2-3 h until stabilized recording every 3 min. 
6. Measure temperature with probe and loggers throughout.  
7. Export the data and calculate mean values. We will create a new calibration dataset from these measurements.  

INSERT OTHER STEPS HERE WHEN WE DO IT  


# 17 September 2022  

### *Histology time series sampling - backreef Acropora and Pocillopora*  

We conducted the September time point of histology time series sampling at Mahana (17°29'12.5"S 149°53'16.1"W) for *Acropora pulchra* and *Pocillopora* sp. Protocols and more information on this project can be found on [GitHub](https://github.com/daniellembecker/Gametogenesis/tree/main/protocols) and on the [E5 public notebook](https://urol-e5.github.io/).  

This time series tracks gametogenesis through histology associated with Danielle's PhD research. For this sampling we collected fragments of n=40 *A. pulchra* and n=10 *Pocillopora* tagged colonies at Mahana site. We collected fragments from 0900-1300. 

Fragments were then processed from 1400-1600 in the lab. After collecting the samples, Danielle, Ariana, and Pierrick followed the [Same Day Sample Processing Protocol](https://github.com/daniellembecker/Gametogenesis/blob/main/protocols/2021-12-26-Sample_Same_Day_Processing_Protocol.md) to prepare samples for histological analysis, molecular and physiological downstream processing. Clippings of each fragment were preserved for histology, molecular samples (flash frozen and in DNA RNA shield), and physiology.  

### *Sampling*  

We again checked for spawning by stopping the water and pumps in all tanks at 2000 to watch for setting in collected colonies. We did not see setting by 2300. There was no spawning. 

# 16 September 2022  

### *Colony collection*    

Parent colonies were collected from the Mahana site on Moorea's north shore (17°29'12.5"S 149°53'16.1"W).  

![field site](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/Moorea2022/mahana.png) 

We first collected n=24 *Acropora. pulchra* colonies for Danielle's project at 0800-0900 at the Mahana field site. We transported these colonies back to Gump station at 1000. Parental colonies had previously experienced a marine heatwave scenario or ambient diel conditions in March-April 2022. All information on the heatwave experiment and data can be found [here](https://github.com/daniellembecker/Gametogenesis/tree/main/heatwave). 

We then collected n=25 colonies (approx. 15-30cm in diameter) from large patches of *Acropora pulchra* from the Mahana field site between 1000-1200. We have a year-long project that has been tracking histology and physiology of colonies (Danielle's project) at this site and therefore we have temperature and light data from this site as environmental context. The colonies collected are considered "wildtype" and were not taken from tagged colonies. Colonies were collected from different patches to obtain a genetically diverse sample. 

Ariana and Danielle used a hammer and chisel to remove large fragments (15-30cm diameter) from these patches. Ariana Danielle and Pierrick transported them back to the boat where they were kept in coolers and buckets filled with seawater. Colonies were then transported back to Gump station at 1230 and placed in 2 flow through seawater tanks (n=12-13 colonies per tank). Tanks are equipped with 3 small recirculating pumps for water flow and had high water flow.   

![parents tank](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/Moorea2022/colonies_tank.jpeg) 

![parents tank 2](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/Moorea2022/parent_tank.jpeg) 

These colonies will be used to collect wildtype gametes to rear to larval stages for Ariana's project. Spawn will be collected from each of the 2 tanks with wildtype colonies.  

### *Launching loggers*  

Temperature loggers were launched to track temperature in the adult parent tanks (n=2 tanks for Danielle and n=2 tanks for Ariana; N=4 loggers) and larval conicals (n=4 conicals per blue tank; N=16 loggers).  

Black Onset Hobo loggers were launched in each parent tank at 1900. Pendant Onset Hobo loggers were launched in larval conicals at 2200. Loggers are set to collect data every 10 minutes. At the end of the experiment we will need to do a group calibration to calculate offsets for the loggers.  The loggers and location are as follows:  

| Logger SSN  | Tank |
| ----------- | ----------- |
| 20444040    | DBP Parent 2 |
| 20719656    | Wildtype Parent 2 |
| 20444033    | Wildtype Parent 1 |
| 20946644    | Wildtype Parent 1 |
| 21002982    | Conical A1 |
| 21002974    | Conical A2 |
| 21002979    | Conical A3 |
| 21002972    | Conical A4 |
| 21002981    | Conical B1 |
| 21335980    | Conical B2 |
| 21335987    | Conical B3 |
| 20937873    | Conical B4 |
| 21002977    | Conical C1 |
| 21335982    | Conical C2 |
| 20937872    | Conical C3 |
| 21002980    | Conical C4 |
| 21335986    | Conical D1 |
| 21335983    | Conical D2 |
| 21002975    | Conical D3 |
| 21335984    | Conical D4 |  

Loggers were placed at the bottom of tanks. Pendants were positioned upright to measure light as well as temperature.  

### *Apex system* 
 
Apex temperature probes were calibrated at 2230 using a digital thermometer (Traceable). Temperatures were between 26.1 and 26.4 for all tanks at the time of calibration.    

### *Tank set up*  

We have been having problems with the water flow not keeping up with demand in the spawning systems with low flow coming from the cooler float valve. Ariana fixed this by cleaning the float valves to remove shell fragments and cutting tubing to go straight from the outlet to the float valve with no slack. To keep the system clean, we will purge the line before starting water flow to the coolers and will regularly check the float valves for shells or debris.  

Tanks are now set up and running! We will add temperature control in the next few days.   

### *Spawning*  

Spawning equipment and workflows were set up. Danielle prepared supplies and equipment to measure eggs per bundle, take samples of sperm and eggs for molecular responses, measure fecundity on sample fragments, and fertilize by parent treatment history. See Danielle's [notebook](https://daniellembecker.github.io/DanielleBecker_Lab_Notebook/) and the [E5 public notebook](https://urol-e5.github.io/) for more information on her project. 

Ariana gathered supplies needed for wildtype spawning. The protocol for Ariana's spawning work will be as follows:  

See a detailed post on this [protocol here](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/Montipora-capitata-spawning-and-larval-rearing-protocol/). 

1. Use gravy separators and mesh bottom cups to collect floating bundles from each of the large parent tanks. Count the number of parents that released to contribute to the gamete pool. Gather gametes into a bucket/bowl. 
2. Once all bundles are released, pool all gametes and gently mix. 
3. Alloquot 5 mL of bundles into each of 16-20 50 mL falcon tubes. 
4. Gently add 40 mL of filtered (50 um) seawater to each falcon tube. 
5. Allow falcon tubes to rest gently sideways with occassional gentle rocking for approx. 1 hour or until bundles are broken apart and eggs and sperm are allowed to mix. 
6. Remove the sperm water from the falcon tube using a pipette and add 10 mL of new filtered seawater to rinse the eggs. Remove this water. 
7. Pour the eggs into the conicals. Use 1 falcon tube for each conical. Have low water flow in the conicals. 
8. Place extra material if available into the 5 gallon buckets with filtered seawater using the same proceedure. 
9. Check the concial cultures every hour. Remove dead material with a pipette, gently mix, and use kim wipes to remove lipid films.  
10. Gently increase the water flow the next morning and continue hourly checks and cleaning.  

We stopped the water and pumps in all tanks at 2000 to watch for setting in collected colonies. We did not see setting by 2300. There was no spawning. Red lights only were used to look at colonies.   
 
# 15 September 2022 

### *Set up*  

Today we continued setting up tanks and unpacking and organizing supplies. We continued assembling portable spawning systems in preparation for spawning. We made and purchased supplies for banjo filters for larval tanks and Gump staff assisted us with installing filtration systems for our water lines. 

This is what our spawninig sytems look like and a full protocol is [available here](https://docs.google.com/spreadsheets/d/1IPcQGoKwVJidsff5xrbiKTevRrtkKmqbUoALj0bd230/edit#gid=0). These systems (pictured below) contain a cooler that acts as a header tank and can be controlled with heaters in combination with a Neptune Apex system. We will set this up in the comiing days. The water in the cooler is sourced from a seawater output. The water in the cooler is then pumped out to a PVC manifold that then delivers water to each larval tank. The larval tanks (i.e., conicals or "squaricals") each have a banjo filter to allow for low thorugh seawater with large surface area to prevent embryos being damaged or lost on the drain. We set up four of these systems this week.   
 
![Tanks1](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/Moorea2022/tanks.jpeg) 

![Tanks2](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/Moorea2022/conicals.jpeg) 

![Filters](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/Moorea2022/filters.jpeg) 

We set up 4 portable spawning systems each with a cooler and 6 larval conicals in a blue tank. Two of these tanks are for Ariana (labelled A and B) and two are for Danielle (labelled C and D). This provides each of us with a high and ambient temperature system. Tanks A and C will be the high temperature tanks and tanks B and D will be ambient temperature.   

We also talked through plans for spawninig and discussed potential protocols for measuring reproductive metrics for Danielle. 

Finally, we prepared 6 extra buckets in a blue tank as extra space in case we have more gametes than can fit in the conicals. These buckets will contain seawater and sit in a water bath. They will not be flow through and will require water changes during embryo rearing.  

### *Apex system* 

Ariana set up two units of Apex systems in preparation for temperature control in larval conicals. There is one Apex unit for Ariana and one unit for Danielle. This will allow for temperature monitoring in n=3 in each tank (A, B, C, D). The temperature of the conicals will be controlled based on the temperature readings in one high temperature conical. Each Apex system consists of the following parts: 

1 Energy Bar AB8
1 Apex control unit 
5 PM1 modules 
6 Aquabus cables 
6 Apex temperature probes 

Apex units were connected to WiFi and created virtual outlets for each stage of temperature control. The Apex format looks like this:  

![Filters](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/images/NotebookImages/Moorea2022/apex_example.png) 

We have the Apex systems named "AH Moorea Sq" and "DBP Moorea Sq".  

# 14 September 2022  

### *Setup*  

Today we arrived in Moorea, French Polynesia. Goals for the day were to unpack supplies, assess space and materials and start experimental set up. 

We will be constructing portable spawning systems developed in the Putnam Lab by Ariana Huffmyer and Jill Ashey. A full protocol for these systems and parts list is [available here](https://docs.google.com/spreadsheets/d/1IPcQGoKwVJidsff5xrbiKTevRrtkKmqbUoALj0bd230/edit#gid=0).  

To prepare for spawning we cleaned and prepared tanks and started seawater flow. We will be using 4 large tanks to hold adult/parent colonies and 4 smaller tanks to hold the spawning systems. Ariana and Pierrick cleaned and filled tanks and constructed drain PVC systems for each tank. Danielle prepared the bins and tanks in water tables that she is using for her project. We unpacked an dorganized the lab equipment.  

