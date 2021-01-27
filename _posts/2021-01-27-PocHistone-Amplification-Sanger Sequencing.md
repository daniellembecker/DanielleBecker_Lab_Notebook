---
layout: post
title: PocHistone Amplification and Sanger Sequencing Prep Protocol
category: [ PCR , Protocol ]
tag: [ Pocillopora, PocHistone ]
---
# Protocol for PocHistone Amplification and Sanger Sequencing for Determining *Pocillopora* Species

#### Goal:
Amplify the histone 3 region region  of _Pocillopora_ coral DNA using PCR (polymerase chain reaction). Then use amplification product for Sanger Sequencing and determination of _Pocillopora_ coral species between _P. eydouxi_ and _P. meandrina_. Primer sequences and concept developed in [Johnston et. al 2018](https://peerj.com/articles/4355/).

Information on Sanger sequencing at the URI GSC is [here](https://web.uri.edu/gsc/sanger_sequencing/). Sequencing takes place after 10am on Tuesdays and Thursdays.

### Process

**Followed the [PocHistone Protocol](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/_posts/2021-01-25-pochistone-protocol.md) exactly. See that for in depth details for this protocol.**

#### DNA Dilutions 2020-08-27

- Arranged samples in a 96 well plate randomly with 3 spots for control reactions
- Allocated DNA for 10ng/ul in 10ul of nuclease free water, so 100ng of DNA. If samples had concentrations below 10ng/ul to start with, just 10ul of those samples was added to each well

![](https://raw.githubusercontent.com/meschedl/MESPutnam_Open_Lab_Notebook/master/images/dilution-plate.png)

#### Histone Amplifications 2021-01-25

- 32 samples plus three negative controls is 35 reactions, use an additional 3 for error
- Made a master mix for 38 samples and 3 reactions each:
  - 1900ul Phusion master mix
  - 49.4ul FatP6.1 10uM Primer
  - 49.4ul RORF 10uM Primer
  - 1672ul nuclease free water
- Added 97ul of master mix into 35 wells in a new plate
- Used a multichannel to add 3ul of DNA from the dilution plate in to the same orientation wells in the plate with the master mix
- Covered plate and vortexed and spun down
- Separated plate out 2 times into 3 separated reaction mixes each with 33ul
- Covered plates, spun down, and placed in three thermocyclers FatP6.1 RORF program

#### 1X Bead Cleanup and Quantification 2021-01-25

- Combined triplicate reactions back together
- Added 1X (100ul) beads to each well
- Followed bead cleanup protocol
- Resuspended and eluted DNA in 50ul ultra-pure water and removed into a new plate (same orientation)
- dsDNA broad range Qubit assay for 35 samples (n# 40)

**Sample.ID**|**Qubit reading 1 (ng/ul)**|**Qubit reading 2 (ng/ul)**|**average DNA (ng/ul)**
:-----:|:-----:|:-----:|:-----:
C22|38.2|38.4|38.3
E10|49|49|49
E2|49|49.4|49.2
C19|44.6|44.6|44.6
C25|49.6|49.6|49.6
C29|47.8|47.8|47.8
C21|51.8|50.8|51.3
E4|48.6|49.2|48.9
E12|52.8|53.8|53.3
E1|48|47.8|47.9
C27|52.6|53.8|53.2
E6|49.2|49.2|49.2
E8|44.8|45.4|45.1
C20|54|54|54
C28|52.6|52.2|52.4
E9|54.4|54.4|54.4
E14|52.2|51.8|52
E15|49.6|48.4|49
C17|52.6|53|52.8
E3|47.6|47.6|47.6
C30|55.2|55.2|55.2
C24|45.8|45.6|45.7
E16|48|47.6|47.8
C31|38.4|38|38.2
E7|44.8|45.4|45.1
C23|43.6|43.6|43.6
E11|45.8|45.8|45.8
E13|40|40|40
C18|35|35.2|35.1
E5|43.6|43.8|43.7
C32|33|33.4|33.2
C26|32.8|32.8|32.8


#### 1% Gel To Confirm Bands 2021-01-26

- Made a 1% gel in the medium gel box
  ![](https://raw.githubusercontent.com/daniellembecker/DanielleBecker_Lab_Notebook/master/images/20210126_PocHistone_Gel.png)
- Bands are strong but some smearing is noticeable which could be an issue. It’s not contamination because our controls are blank. So it’s either an imaging artifact or something that amplified non-specific in the sample. We will be sequencing four samples from a sample size of each treatment and one with smearing and one without to see if there is more of a problem that needs troubleshooting.

#### Dilution and Sequencing Prep 2021-01-27

- Diluted each DNA amplification by 1:10 (2ul of DNA and 18ul of ultra-pure water)
- Created plate with 16.7ng of DNA for each sample and ultra-pure water up to 10ul
 ![](https://raw.githubusercontent.com/daniellembecker/DanielleBecker_Lab_Notebook/master/images/20210126_16.7ng_H20_10ul.png)
- Added 2ul of 3.2uM RORF primer to each well
- Covered, vortexed, spun down plate
- Brought upstairs for sequencing the next day

**Spreadsheet for Sequencing**

| Sample IDa | Well  (GSC use only) | Template Typeb | A. Template Size (bases) | B. Template Stock Conc. (ng/µl) | C. PCR template: ng needed = ((A ÷ 100) x 1.25)  x 2 | D. PCR template: Volume = (C ÷ B) µl | F. Volume PCR-H20 needed (10 minus D  or E) µl | G. Volume primer needed 1  µl per reaction |
|------------|----------------------|----------------|--------------------------|---------------------------------|------------------------------------------------------|--------------------------------------|------------------------------------------------|--------------------------------------------|
| HPD33      |                      | PCR            | 699                      | 3.83                            | 17.475                                               | 4.36                                 | 5.64                                           | 2                                          |
| HPD34      |                      | PCR            | 699                      | 4.9                             | 17.475                                               | 3.41                                 | 6.59                                           | 2                                          |
| HPD35      |                      | PCR            | 699                      | 4.92                            | 17.475                                               | 3.39                                 | 6.61                                           | 2                                          |
| HPD36      |                      | PCR            | 699                      | 4.46                            | 17.475                                               | 3.74                                 | 6.26                                           | 2                                          |
| HPD37      |                      | PCR            | 699                      | 4.96                            | 17.475                                               | 3.37                                 | 6.63                                           | 2                                          |
| HPD38      |                      | PCR            | 699                      | 4.78                            | 17.475                                               | 3.49                                 | 6.51                                           | 2                                          |
| HPD39      |                      | PCR            | 699                      | 5.13                            | 17.475                                               | 3.26                                 | 6.74                                           | 2                                          |
| HPD40      |                      | PCR            | 699                      | 4.89                            | 17.475                                               | 3.42                                 | 6.58                                           | 2                                          |
| HPD41      |                      | PCR            | 699                      | 5.33                            | 17.475                                               | 3.13                                 | 6.87                                           | 2                                          |
| HPD42      |                      | PCR            | 699                      | 4.79                            | 17.475                                               | 3.49                                 | 6.51                                           | 2                                          |
| HPD43      |                      | PCR            | 699                      | 5.32                            | 17.475                                               | 3.14                                 | 6.86                                           | 2                                          |
| HPD44      |                      | PCR            | 699                      | 4.92                            | 17.475                                               | 3.39                                 | 6.61                                           | 2                                          |
| HPD45      |                      | PCR            | 699                      | 4.51                            | 17.475                                               | 3.70                                 | 6.30                                           | 2                                          |
| HPD46      |                      | PCR            | 699                      | 5.4                             | 17.475                                               | 3.09                                 | 6.91                                           | 2                                          |
| HPD47      |                      | PCR            | 699                      | 5.24                            | 17.475                                               | 3.19                                 | 6.81                                           | 2                                          |
| HPD48      |                      | PCR            | 699                      | 5.44                            | 17.475                                               | 3.07                                 | 6.93                                           | 2                                          |
| HPD49      |                      | PCR            | 699                      | 5.2                             | 17.475                                               | 3.21                                 | 6.79                                           | 2                                          |
| HPD50      |                      | PCR            | 699                      | 4.9                             | 17.475                                               | 3.41                                 | 6.59                                           | 2                                          |
| HPD51      |                      | PCR            | 699                      | 5.28                            | 17.475                                               | 3.16                                 | 6.84                                           | 2                                          |
| HPD52      |                      | PCR            | 699                      | 4.76                            | 17.475                                               | 3.51                                 | 6.49                                           | 2                                          |
| HPD53      |                      | PCR            | 699                      | 5.52                            | 17.475                                               | 3.03                                 | 6.97                                           | 2                                          |
| HPD54      |                      | PCR            | 699                      | 4.57                            | 17.475                                               | 3.65                                 | 6.35                                           | 2                                          |
| HPD55      |                      | PCR            | 699                      | 4.78                            | 17.475                                               | 3.49                                 | 6.51                                           | 2                                          |
| HPD56      |                      | PCR            | 699                      | 3.82                            | 17.475                                               | 4.37                                 | 5.63                                           | 2                                          |
| HPD57      |                      | PCR            | 699                      | 4.51                            | 17.475                                               | 3.70                                 | 6.30                                           | 2                                          |
| HPD58      |                      | PCR            | 699                      | 4.36                            | 17.475                                               | 3.83                                 | 6.17                                           | 2                                          |
| HPD59      |                      | PCR            | 699                      | 4.58                            | 17.475                                               | 3.65                                 | 6.35                                           | 2                                          |
| HPD60      |                      | PCR            | 699                      | 4                               | 17.475                                               | 4.18                                 | 5.83                                           | 2                                          |
| HPD61      |                      | PCR            | 699                      | 3.51                            | 17.475                                               | 4.76                                 | 5.24                                           | 2                                          |
| HPD62      |                      | PCR            | 699                      | 4.37                            | 17.475                                               | 3.82                                 | 6.18                                           | 2                                          |
| HPD63      |                      | PCR            | 699                      | 3.32                            | 17.475                                               | 5.03                                 | 4.97                                           | 2                                          |
| HPD64      |                      | PCR            | 699                      | 3.28                            | 17.475                                               | 5.09                                 | 4.91                                           | 2                                          |
