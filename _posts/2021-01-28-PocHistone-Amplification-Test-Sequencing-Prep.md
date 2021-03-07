---
layout: post
title: PocHistone Amplification and Test Sanger Sequencing Prep Protocol
category: [ PCR , Protocol ]
tag: [ Pocillopora, PocHistone ]
---
## Protocol for PocHistone Amplification and Test Sanger Sequencing Prep

#### Goal:
Amplify the histone 3 region of _Pocillopora_ coral DNA using PCR (polymerase chain reaction). Then use amplification product for Sanger Sequencing and determination of _Pocillopora_ coral species between _P. eydouxi_ and _P. meandrina_ on 4 test samples out of 32 samples. Primer sequences and concept developed in [Johnston et. al 2018](https://peerj.com/articles/4355/).

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

#### Dilution and Sequencing Prep 2021-01-27

- Used samples E12 (small band underneath); E8 (no smearing); C17 (small band above); C30 (smearing)
- For sequencing I need 9ul (4 * 2 +1 for error) of 3.2uM RORF primer
- The mtORF product should be ~669 bp, so 669 / 100 * 1.25 * 2 = 16.7ng of each sample
- 16.7ng is still less than 1ul for every sample, so I did a 1:10 dilution for each one
- Added 18ul of ultra pure water to 4 new strip tubes, and added 2ul of each cleaned DNA to their corresponding well
- Pipetted 16.7ng of each sample from diluted stock to new strip tubes
- Increased volume in tubes to 10ul with ultra pure water

|Sample|GSC sample code|Dilution conc.|ng needed|Volume diluted PCR product|volume ultra pure H2O|
|---|---|---|---|---|---|
|E12|HPD33|5.33|16.7|3.13|6.87|
|E8|HPD34|4.51|16.7|3.70|6.30|
|C17|HPD35|5.28|16.7|3.16|6.84|
|C30|HPD36|5.52|16.7|3.03|6.97|

- Made 9ul of diluted 10uM stock RORF to 3.2uM
- Added 6.12ul of ultra pure water to 1 PCR tube
- Added 2.88ul of 10uM RORF primer to tube
- Vortex and spin down to mix
- Added 2ul of 3.2uM primer to each of the 8 strip tubes diluted for sequencing
- Stored prepared samples in -20 freezer

#### Sequencing Submission 2021-01-28
- Submitted four test samples for sequencing at 9:40am

**Spreadsheet for Sequencing**

| Sample IDa | Well  (GSC use only) | Template Typeb | A. Template Size (bases) | B. Template Stock Conc. (ng/µl) | C. PCR template: ng needed = ((A ÷ 100) x 1.25)  x 2 | D. PCR template: Volume = (C ÷ B) µl | F. Volume PCR-H20 needed (10 minus D  or E) µl | G. Volume primer needed 1  µl per reaction |
|------------|----------------------|----------------|--------------------------|---------------------------------|------------------------------------------------------|--------------------------------------|------------------------------------------------|--------------------------------------------|
| HPD33      |                      | PCR            | 669                      | 5.33                            | 16.725                                               | 3.13                                 | 6.87                                           | 2                                          |
| HPD34      |                      | PCR            | 669                      | 4.51                            | 16.725                                               | 3.70                                 | 6.30                                           | 2                                          |
| HPD35      |                      | PCR            | 669                      | 5.28                            | 16.725                                               | 3.16                                 | 6.84                                           | 2                                          |
| HPD36      |                      | PCR            | 669                      | 5.52                            | 16.725                                               | 3.03                                 | 6.97                                           | 2                                          |

#### Sequences Recieved 2021-01-29

- All sequences looked good, was able to align and delinate to species using Geneious Prime.
