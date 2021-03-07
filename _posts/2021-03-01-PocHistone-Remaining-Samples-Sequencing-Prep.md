---
layout: post
title: PocHistone Sanger Sequencing Prep Protocol
category: [ PCR , Protocol ]
tag: [ Pocillopora, PocHistone ]
---
## Protocol for PocHistone Sanger Sequencing Prep for the Remaining 28 Samples

#### Goal:
Amplify the histone 3 region of _Pocillopora_ coral DNA using PCR (polymerase chain reaction). Then use amplification product for Sanger Sequencing and determination of _Pocillopora_ coral species between _P. eydouxi_ and _P. meandrina_ on the remaining 28 samples out of 32 samples. Primer sequences and concept developed in [Johnston et. al 2018](https://peerj.com/articles/4355/).

Information on Sanger sequencing at the URI GSC is [here](https://web.uri.edu/gsc/sanger_sequencing/). Sequencing takes place after 10am on Tuesdays and Thursdays.

### Process

**Followed the [PocHistone Protocol](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/_posts/2021-01-25-pochistone-protocol.md) exactly. See that for in depth details for this protocol.**

#### DNA Dilutions 2020-08-27

DNA dilutions for samples completed by Maggie on 2020-08-27. See the [notebook post](https://github.com/meschedl/MESPutnam_Open_Lab_Notebook/blob/master/_posts/2020-08-27-mtORF-protocol.md) for in depth details for this protocol.

#### Histone Amplifications 2021-01-25

- 32 samples plus three negative controls is 35 reactions, use an additional 3 for error
- Made a master mix for 38 samples and 3 reactions each:
  - 1900ul Phusion master mix
  - 49.4ul working stock PocHistoneF 10uM Primer
  - 49.4ul working stock PocHistoneR 10uM Primer
  - 1672ul nuclease free water
- Added 97ul of master mix into 35 wells in a new plate
- Used a multichannel to add 3ul of DNA from the dilution plate in to the same orientation wells in the plate with the master mix
- Covered plate and vortexed and spun down
- Separated plate out 2 times into 3 separated reaction mixes each with 33ul
- Covered plates, spun down, and placed in three thermocyclers PocHistone program

#### 1X Bead Cleanup and Quantification 2021-01-25

- Combined triplicate reactions back together
- Added 1X (100ul) beads to each well
- Followed bead cleanup protocol
- Resuspended and eluted DNA in 50ul ultra-pure water and removed into a new plate (same orientation)

#### Broad Range Qubit 2020-01-25

- dsDNA broad range Qubit assay for 35 samples (n# 40)
- Followed [qubit protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/Qubit-Protocol/)
- Included 32 samples, 3 controls, and 2 standards

**Sample.ID**|**Qubit reading 1 (ng/ul)**|**Qubit reading 2 (ng/ul)**|**average DNA (ng/ul)**
:-----:|:-----:|:-----:|:-----:
std. 1|176.98|NA|NA
std. 1|20368.25|NA|NA
E10|49|49|49
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

#### Dilution (2021-01-26) and Sequencing Prep for Remaining 28 Samples (2021-03-01)

###### All samples were diluted at the same time on 2021-01-26 
- Prepped the remaining 28 samples for sequencing 
- For sequencing I needed 57ul (28 * 2 +1 for error) of 3.2uM FORF primer
- The mtORF product should be ~669 bp, so 669 / 100 * 1.25 * 2 = 16.7ng of each sample
- 16.7ng is still less than 1ul for every sample, so a 1:10 dilution was done for each sample on 2021-01-26
- Added 18ul of ultra pure water to the 96-well plate for 32 samples, and added 2ul of each cleaned DNA to their corresponding well
- Pipetted 16.7ng of each sample from diluted stock to the 96-well plate for 32 samples
- Increased volume in tubes to 10ul with ultra pure water

| sample | GSC sample code | dilution.conc | ng.needed | volume diluted PCR product | volume ultra pure H2O |
|--------|-----------------|---------------|-----------|----------------------------|-----------------------|
| C22    | HPD37           | 3.83          | 16.7      | 4.36                       | 5.64                  |
| E10    | HPD38           | 4.9           | 16.7      | 3.41                       | 6.59                  |
| E2     | HPD39           | 4.92          | 16.7      | 3.39                       | 6.61                  |
| C19    | HPD40           | 4.46          | 16.7      | 3.74                       | 6.26                  |
| C25    | HPD41           | 4.96          | 16.7      | 3.37                       | 6.63                  |
| C29    | HPD42           | 4.78          | 16.7      | 3.49                       | 6.51                  |
| C21    | HPD43           | 5.13          | 16.7      | 3.26                       | 6.74                  |
| E4     | HPD44           | 4.89          | 16.7      | 3.42                       | 6.58                  |
| E1     | HPD45           | 4.79          | 16.7      | 3.49                       | 6.51                  |
| C27    | HPD46           | 5.32          | 16.7      | 3.14                       | 6.86                  |
| E6     | HPD47           | 4.92          | 16.7      | 3.39                       | 6.61                  |
| C20    | HPD48           | 5.4           | 16.7      | 3.09                       | 6.91                  |
| C28    | HPD49           | 5.24          | 16.7      | 3.19                       | 6.81                  |
| E9     | HPD50           | 5.44          | 16.7      | 3.07                       | 6.93                  |
| E14    | HPD51           | 5.2           | 16.7      | 3.21                       | 6.79                  |
| E15    | HPD52           | 4.9           | 16.7      | 3.41                       | 6.59                  |
| E3     | HPD53           | 4.76          | 16.7      | 3.51                       | 6.49                  |
| C24    | HPD54           | 4.57          | 16.7      | 3.65                       | 6.35                  |
| E16    | HPD55           | 4.78          | 16.7      | 3.49                       | 6.51                  |
| C31    | HPD56           | 3.82          | 16.7      | 4.37                       | 5.63                  |
| E7     | HPD57           | 4.51          | 16.7      | 3.7                        | 6.3                   |
| C23    | HPD58           | 4.36          | 16.7      | 3.83                       | 6.17                  |
| E11    | HPD59           | 4.58          | 16.7      | 3.65                       | 6.35                  |
| E13    | HPD60           | 4             | 16.7      | 4.18                       | 5.83                  |
| C18    | HPD61           | 3.51          | 16.7      | 4.76                       | 5.24                  |
| E5     | HPD62           | 4.37          | 16.7      | 3.82                       | 6.18                  |
| C32    | HPD63           | 3.32          | 16.7      | 5.03                       | 4.97                  |
| C26    | HPD64           | 3.28          | 16.7      | 5.09                       | 4.91                  |

###### Sequencing prep on 2021-03-01
- Added 16.7 ng of DNA following the volume of the diluted PCR product for each sample above to each corresponding well in a 96-well plate
- Added the volume of ultra pure H2O corresponding to each sample above into the 96-well plate wells

![](https://raw.githubusercontent.com/daniellembecker/DanielleBecker_Lab_Notebook/master/images/28_samples_sequence_layout.png)


- Made 57ul of diluted 10uM stock FORF to 3.2uM
- Added 38.76ul of ultra pure water to 1 PCR tube
- Added 18.24ul of 10uM FORF primer to tube
- Vortexed and centrifuged to mix
- Added 2ul of 3.2uM primer to each of the 28 wells of the 96-well plate

#### Sequencing Submission 2021-03-01
- Submitted remaining 28 samples for processing in the GSC lab at 3:00 pm

**Spreadsheet for Sequencing**

| Sample IDa | Well           | Template | A. Template Size (bases) | B. Template Stock Conc. (ng/µl) | C. PCR template: ng needed = ((A ÷ 100) x 1.25) x 2 | D. PCR template: Volume = (C ÷ B) µl | F. Volume PCR-H20 needed (10 minus D or E) µl | G. Volume primer needed 1 µl per reaction |
|------------|----------------|----------|--------------------------|---------------------------------|-----------------------------------------------------|--------------------------------------|-----------------------------------------------|-------------------------------------------|
|            | (GSC use only) | Typeb    |                          |                                 |                                                     |                                      |                                               |                                           |
| HPD37      |                | PCR      | 669                      | 3.83                            | 16.725                                              | 4.36                                 | 5.64                                          | 2                                         |
| HPD38      |                | PCR      | 669                      | 4.9                             | 16.725                                              | 3.41                                 | 6.59                                          | 2                                         |
| HPD39      |                | PCR      | 669                      | 4.92                            | 16.725                                              | 3.39                                 | 6.61                                          | 2                                         |
| HPD40      |                | PCR      | 669                      | 4.46                            | 16.725                                              | 3.74                                 | 6.26                                          | 2                                         |
| HPD41      |                | PCR      | 669                      | 4.96                            | 16.725                                              | 3.37                                 | 6.63                                          | 2                                         |
| HPD42      |                | PCR      | 669                      | 4.78                            | 16.725                                              | 3.49                                 | 6.51                                          | 2                                         |
| HPD43      |                | PCR      | 669                      | 5.13                            | 16.725                                              | 3.26                                 | 6.74                                          | 2                                         |
| HPD44      |                | PCR      | 669                      | 4.89                            | 16.725                                              | 3.42                                 | 6.58                                          | 2                                         |
| HPD45      |                | PCR      | 669                      | 4.79                            | 16.725                                              | 3.49                                 | 6.51                                          | 2                                         |
| HPD46      |                | PCR      | 669                      | 5.32                            | 16.725                                              | 3.14                                 | 6.86                                          | 2                                         |
| HPD47      |                | PCR      | 669                      | 4.92                            | 16.725                                              | 3.39                                 | 6.61                                          | 2                                         |
| HPD48      |                | PCR      | 669                      | 5.4                             | 16.725                                              | 3.09                                 | 6.91                                          | 2                                         |
| HPD49      |                | PCR      | 669                      | 5.24                            | 16.725                                              | 3.19                                 | 6.81                                          | 2                                         |
| HPD50      |                | PCR      | 669                      | 5.44                            | 16.725                                              | 3.07                                 | 6.93                                          | 2                                         |
| HPD51      |                | PCR      | 669                      | 5.2                             | 16.725                                              | 3.21                                 | 6.79                                          | 2                                         |
| HPD52      |                | PCR      | 669                      | 4.9                             | 16.725                                              | 3.41                                 | 6.59                                          | 2                                         |
| HPD53      |                | PCR      | 669                      | 4.76                            | 16.725                                              | 3.51                                 | 6.49                                          | 2                                         |
| HPD54      |                | PCR      | 669                      | 4.57                            | 16.725                                              | 3.65                                 | 6.35                                          | 2                                         |
| HPD55      |                | PCR      | 669                      | 4.78                            | 16.725                                              | 3.49                                 | 6.51                                          | 2                                         |
| HPD56      |                | PCR      | 669                      | 3.82                            | 16.725                                              | 4.37                                 | 5.63                                          | 2                                         |
| HPD57      |                | PCR      | 669                      | 4.51                            | 16.725                                              | 3.7                                  | 6.3                                           | 2                                         |
| HPD58      |                | PCR      | 669                      | 4.36                            | 16.725                                              | 3.83                                 | 6.17                                          | 2                                         |
| HPD59      |                | PCR      | 669                      | 4.58                            | 16.725                                              | 3.65                                 | 6.35                                          | 2                                         |
| HPD60      |                | PCR      | 669                      | 4                               | 16.725                                              | 4.18                                 | 5.83                                          | 2                                         |
| HPD61      |                | PCR      | 669                      | 3.51                            | 16.725                                              | 4.76                                 | 5.24                                          | 2                                         |
| HPD62      |                | PCR      | 669                      | 4.37                            | 16.725                                              | 3.82                                 | 6.18                                          | 2                                         |
| HPD63      |                | PCR      | 669                      | 3.32                            | 16.725                                              | 5.03                                 | 4.97                                          | 2                                         |
| HPD64      |                | PCR      | 669                      | 3.28                            | 16.725                                              | 5.09                                 | 4.91                                          | 2                                         |

#### Sequences Recieved 2021-03-03

- Many of the sequences showed peaks and overlap in their chromatograms, but some were good. During troubleshooting, I organized the sequences that did not exhibit overlap and peaks in their chromatograms. 

Example of a successful and unsuccessful chromatogram for two samples in the same sequence run:

![](https://raw.githubusercontent.com/daniellembecker/DanielleBecker_Lab_Notebook/master/images/20210303_chromatogram_successful.png)


![](https://raw.githubusercontent.com/daniellembecker/DanielleBecker_Lab_Notebook/master/images/20210303_chromatogram_unsuccessful.png)

- I was able to align and delinate to species for samples: E10, C25, E6, C24, E7, and E13. 
- Upon discussion, Maggie and Hollie suggested running a gel for longer could show if our gel was masking bacterial, symbiont, or coral gene duplication DNA that the primer worked on as well as the targeted region.
    - Disscussion points: since the gel showed very bright bands it could mean that it was probably a little too much DNA in each well. Which may have been masking something that is a little larger in size. If we had a lot of DNA, it could be too much to all be moving at the same pace, she suggested that could be why we got some fuzzy/glowy bands. Adding less DNA to the gel (like using the 1:10 dilution) and running it for longer, maybe 1.5 hours, will make less bright bands and spread them out more if there are two. 
- We then reached out to Dr. Erika Johnston who developed and worked heavily with these primers, she suggested that maybe the forward primer could be picking up non-specific DNA that could have been amplified in the smear, possibly from bacteria, etc. She suggested using the reverse primer instead beacuse it is much more specific than the forward. 
- I will be prepping the remaining 22 samples for sequencing using the reverse primer to see if this was the main issue from this round of sequencing.

Example chromatograms
