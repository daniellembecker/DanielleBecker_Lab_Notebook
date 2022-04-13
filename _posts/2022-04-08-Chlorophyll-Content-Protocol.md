---
layout: post
title: Chlorophyll Content Protocol
date: '2022-04-08'
categories: Protocol
tags: Protocol Chlorophyll Physiology
---
This protocol details chlorophyll measurement in coral larval tissue samples. This protocol is modified from the [E5 chlorophyll protocol](https://github.com/urol-e5/protocols/blob/master/2020-01-01-Chlorophyll-Protocol.md).  

**Note that this protocol needs to be optimized prior to use for ASH larval samples**  

# Quantifying Chl-a and Chl-c2 Concentration in Symbiodiniaceae from coral tissues

Original: 20200101  
Last Revised: 20220413 ASH  

Contents  
- [**Materials**](#Materials)   
- [**Protocol**](#Protocol)  
- [**References**](#References)  

1. <a name="Materials"></a> **Materials**
    - 	100% acetone
    - 	flammable safe fridge 4°
    - 	**quartz** 96 well plate
    -   Microcentrifuge (Mo'orea Gump Molecular Lab)
    - 	1ml pipette and tips
    -   1.5 ml microfuge tubes
    -   [Synergy HTX Multi-Mode Microplate Reader](https://www.biotek.com/products/detection-multi-mode-microplate-readers/synergy-htx-multi-mode-reader/) (PPP lab plate reader)  
    -   [Gen5 Software](https://www.biotek.com/products/software-robotics-software/gen5-microplate-reader-and-imager-software/) (PPP Lab computer)  

**2. <a name="Protocol"></a> Protocol**  
  
**Sample Preparation**  
1. Thaw homogenate aliquot. These samples are previously separated using the [fraction separation protocol described here](https://ahuffmyer.github.io/ASH_Putnam_Lab_Notebook/Montipora-capitata-2021-larval-temperature-experiment-sample-processing/).   
2. Centrifuge the aliquot of larval homogenate at 13,000 rpm for 3 minutes to separate the host and Symbiodiniaceae cells.  
3. Remove the supernatant and add to a new tube labeled "Holo-Supernatant".  
4. Add 1 mL of 100% acetone to the pellet in the 1.5 mL microcentrifuge tube and vortex the tubes for 15 sec.  
5. Place the tubes in a fridge in the dark at 4°C for 24 hours.  
6. Vortex the tubes for 15 sec.  
7. Spin the tubes down at 13,000 rpm for 3 minutes in the microcentrifuge to pellet any debris.  
8. Pipette 200µl of sample to duplicate wells of 96-well quartz plate.   
9. Pipette 200µl of acetone blank to duplicate wells.  
10. Cover the plate with silicone pad every 5th sample or so to reduce evaporate as samples are added.  
11. Remove silicon pad.   
12. Follow steps below to measure the extract Absorbance on the [Synergy HTX Multi-Mode Microplate Reader](https://www.biotek.com/products/detection-multi-mode-microplate-readers/synergy-htx-multi-mode-reader/) at 630, 663, and 750 nm in a 96-well quartz plate.
13. Standardize for path length in 200µl of sample in 96-well quartz plate.

**Measure the Absorbance**  

Follow the [Synergy HTX Operating Manual](https://github.com/urol-e5/protocols/blob/master/synergy_htx_manual.pdf) and [Gen5 Software Manual](https://github.com/urol-e5/protocols/blob/master/Gen5_software_manual.pdf) to install the software on your host computer and general operating instructions.

1. Open the Gen5 software on your computer.

![software](https://raw.githubusercontent.com/urol-e5/protocols/master/images/GEN5_chloro.jpeg)

2. Whenever you start Gen5, the Task Manager opens. Load the presaved chlorophyll protocol on the Gen5 software on the HP Putnam Lab computer.  
3. To create a new protocol, follow the 'Getting Started' section in the [Gen5 Software Manual](https://github.com/urol-e5/protocols/blob/master/Gen5_software_manual.pdf).
4. The plate loader should automatically open. Load your plate following the distinctions in the plate loader and you will be all set.
5. To measure Absorbance, select the 'Chlorophyll Protocol' option.

![task.manager](https://raw.githubusercontent.com/urol-e5/protocols/master/images/GEN5_software.jpeg)

6. Press the green start protocol button in the GEN5 software and allow the protocol to run through the 630, 663, and 750 nm.
7. Export the raw GEN5 software data file and a CSV to the Desktop in a designated folder. Also, immediately upload to a google drive folder to have for future use. 
8. Standardize for path length in 200µl of sample in 96-well quartz plate.


**Calculating Chlorophyll Concentration**  

Chlorophyll a and c2 concentrations are calculated from the equations in [Jeffrey and Humphrey 1975](https://reader.elsevier.com/reader/sd/pii/S0015379617307783?token=0937035D38C07F29ADF00F1F2A21F20F221219B1CC11A444A4F84D16B98EC3A6AD941D191BA2135A68C98BA62A0B69FE) after substracting A750nm from all measurements.  

Need to correct for differences in path length of the volume in the 96 well plate compared to the 1cm path length of a cuvette.
[Warren 2007](https://www.tandfonline.com/doi/full/10.1080/01904160802135092?casa_token=RqeUl1Ccg7AAAAAA%3A6SyNAs848qrRk1-Tf1g088xWD10z1Xngb8cmcgRvC3jYSYPugr2cL8QG9wFvrFj7xZF-pqqUozonRg)

4. <a name="References"></a> **References**

    1.  [Jeffrey and Humphrey 1975](https://reader.elsevier.com/reader/sd/pii/S0015379617307783?token=0937035D38C07F29ADF00F1F2A21F20F221219B1CC11A444A4F84D16B98EC3A6AD941D191BA2135A68C98BA62A0B69FE)
    2. [Warren 2007](https://www.tandfonline.com/doi/full/10.1080/01904160802135092?casa_token=RqeUl1Ccg7AAAAAA%3A6SyNAs848qrRk1-Tf1g088xWD10z1Xngb8cmcgRvC3jYSYPugr2cL8QG9wFvrFj7xZF-pqqUozonRg)
    3. [Synergy HTX Operating Manual](https://github.com/urol-e5/protocols/blob/master/synergy_htx_manual.pdf)
    4. [Gen5 Software Manual](https://github.com/urol-e5/protocols/blob/master/Gen5_software_manual.pdf)
