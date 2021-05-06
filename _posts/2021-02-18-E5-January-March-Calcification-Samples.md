---
layout: post
title: E5 Timepoint 1 and 2 Instantaneous Calcification Sample Processing (Samples Collected in 2020)
category: [ Processing, Protocol ]
tag: [ Physiology, Calcification Samples, Mercuric Chloride Fixed Samples ]
---
## Processing the E5 January and March 2020 Instantaneous Calcification Samples

### **Goal**
Process calcification samples collected for the [E5 uROL Epigenetics Project](https://urol-e5.github.io) from January and March 2020 (timepoint 1 and 2). The samples were collected in January (n = X) and March (n = X) 2020 and fixed with 75uL of saturated mercuric chloride for transportation back to URI on November, 22nd 2020.

### **Process**

**Samples were collected, experimental run, and transported back to URI following the [E5 Instantaneous Calcification Protocol](https://github.com/urol-e5/protocols/blob/master/2020-01-01-Instantaneous-Calcification-Protocol.md) in Mo'orea, French Polynesia.**

Upon the return to URI, we wanted to test the total alkalinity signal we would receive from our initial seawater samples before the experimental procedure and after to calculate delta total alkalinity. To do so, titrations will need to be completed for all samples to calculate total alkalinity following the alkalinity anomaly technique ([Chisholm and Gattuso, 1991](https://aslopubs.onlinelibrary.wiley.com/doi/pdf/10.4319/lo.1991.36.6.1232)).

All samples will be processed following the established [Putnam Lab Titrator Protocol](https://github.com/Putnam-Lab/Lab_Management/blob/master/Lab_Resourses/Equipment_Protocols/Titrator_Protocols/Titrator_Protocol.md).

All data processed and collected on the Putnam Lab Titrator will be available on the [Putnam Lab Titrator GitHub Repository](https://github.com/Putnam-Lab/Titrator).


#### **Data and notes for testing pH probe and buffers on 20210221**

**Notes**
- In the initial titration set-up I noticed that the pH calibration was a bit off due to the zero point and slope values we received [here](https://github.com/Putnam-Lab/Titrator/blob/main/Data/pHCalibration.csv). I ran a CRM and found that is was within the < 1% accuracy that we would expect but at a value of -0.82% was a little high especially with the floating pH calibration data. I began to troubleshoot the problem and thought that maybe the pH buffers could be off from prior knowledge. A new round of buffers was ordered to test the pH probe. pH calibration and CRM accuracy information can be found in the [Putnam Lab Titrator Protocol](https://github.com/Putnam-Lab/Lab_Management/blob/master/Lab_Resourses/Equipment_Protocols/Titrator_Protocols/Titrator_Protocol.md).

**Data**

| SampleID | TA               | Mass   | Salinity |
|----------|------------------|--------|----------|
| JUNK1    | 1791.44566335087 | 59.81  | 35       |
| CRM 1    | 2215.55476501538 | 59.561 | 33.417   |



#### **Data and notes for testing variability in CRMs and pH probe 20210228**

**Notes**
- Continued the troubleshooting process and still saw that the pH calibration values were off and the probe was floating in its mV values when measured manually for 40 seconds in each of the 4, 7, and 10 pH buffers. The values were small but still seemed inconsistent. To test how much variability we were seeing in CRMs due to the pH probe, I ran 8 CRMs. There was a lot of variability between the CRMs so a new pH probe was ordered and received on 20210304.

pH buffer 4 manually measured values (mV)

![](https://raw.githubusercontent.com/daniellembecker/DanielleBecker_Lab_Notebook/master/images/pH.probe.4.jpg)

pH buffer 7 manually measured values (mV)

![](https://raw.githubusercontent.com/daniellembecker/DanielleBecker_Lab_Notebook/master/images/pH.probe.7.jpg)

pH buffer 10 manually measured values (mV)

![](https://raw.githubusercontent.com/daniellembecker/DanielleBecker_Lab_Notebook/master/images/pH.probe.10.jpg)

**Data**

CRM Accuracy Data

| Date     | CRM value   | Batch value | % off        | Batch # | Notes               |
|----------|-------------|-------------|--------------|---------|---------------------|
| 20210218 | 2215.547681 | 2234.07     | -0.829084093 | 141     | CRM opened 20210216 |
| 20210221 | 2215.554    | 2234.07     | -0.828801246 | 141     | CRM opened 20210221 |
| 20210228 | 2211.139264 | 2234.07     | -1.026410826 | 141     | CRM opened 20210216 |
| 20210228 | 2213.599263 | 2234.07     | -0.916297939 | 141     | CRM opened 20210216 |
| 20210228 | 2218.811032 | 2234.07     | -0.683012077 | 141     | CRM opened 20210228 |
| 20210228 | 2208.830425 | 2234.07     | -1.129757569 | 141     | CRM opened 20210228 |
| 20210228 | 2203.497378 | 2234.07     | -1.368472003 | 141     | CRM opened 20210228 |
| 20210228 | 2213.354976 | 2234.07     | -0.927232556 | 141     | CRM opened 20210228 |
| 20210228 | 2230.94523  | 2234.07     | -0.13986896  | 141     | CRM opened 20210228 |
| 20210228 | 2222.545765 | 2234.07     | -0.515840388 | 141     | CRM opened 20210228 |

Total Alkalinity Output

| SampleID | TA          | Mass   | Salinity |
|----------|-------------|--------|----------|
| JUNK1    | 932.4357056 | 59.701 | 35       |
| CRM1     | 2211.139264 | 59.861 | 33.417   |
| CRM2     | 2213.599263 | 60.436 | 33.417   |
| CRM3     | 2218.811032 | 60.178 | 33.417   |
| CRM4     | 2208.830425 | 60.307 | 33.417   |
| CRM5     | 2203.497378 | 60.315 | 33.417   |
| CRM6     | 2213.354976 | 60.201 | 33.417   |
| CRM7     | 2230.94523  | 60.413 | 33.417   |
| CRM8     | 2222.545765 | 60.105 | 33.417   |


#### **New pH probe received on 20210303**

#### **New pH probe assembled, installed and tested on 20210307**

#### **Data and notes for testing variability in CRMs and new pH probe 20210314**

**Notes**
- Continued the troubleshooting process and had great pH calibration values with the new probe and buffers before any runs today. Purged all bubbles and checked every box before running one junk and 8 CRM samples. After the run with 8 CRMS, I noticed that the variability from the CRM value was still >1% in a few of the samples. I attributed this to the older CRM bottle that was opened on 20210228. I opened a brand new CRM bottle and ran four samples but there was still variability to the original CRM value but less than in the previous run. I called Mettler Toledo support and they said this could be an offset with the titrant (acid) and they have had this happen. We could measure the titrant against tris buffer to calculate the offset and apply that to our equations. Noticed that the titrant was replaced in May 2019, which seemed very old and think evaporation could have effected its composition. Next round I am going to change out the titrant with a new bottle.

**Data**

CRM Accuracy Data

| Date     | CRM value   | Batch value | % off        | Batch # | Notes               |
|----------|-------------|-------------|--------------|---------|---------------------|
| 20210314 | 2182.69197  | 2234.07     | -2.29975022  | 141     | CRM opened 20210228 |
| 20210314 | 2183.070705 | 2234.07     | -2.282797543 | 141     | CRM opened 20210228 |
| 20210314 | 2179.989317 | 2234.07     | -2.420724647 | 141     | CRM opened 20210228 |
| 20210314 | 2180.118916 | 2234.07     | -2.414923637 | 141     | CRM opened 20210228 |
| 20210314 | 2194.262889 | 2234.07     | -1.781820203 | 141     | CRM opened 20210228 |
| 20210314 | 2194.808799 | 2234.07     | -1.757384566 | 141     | CRM opened 20210228 |
| 20210314 | 2192.134088 | 2234.07     | -1.877108222 | 141     | CRM opened 20210228 |
| 20210314 | 2192.253611 | 2234.07     | -1.871758246 | 141     | CRM opened 20210228 |
| 20210314 | 2205.731211 | 2224.47     | -0.842393439 | 180     | CRM opened 20210314 |
| 20210314 | 2194.005521 | 2224.47     | -1.369516288 | 180     | CRM opened 20210314 |
| 20210314 | 2192.044753 | 2224.47     | -1.457661693 | 180     | CRM opened 20210314 |
| 20210314 | 2196.378938 | 2224.47     | -1.262820467 | 180     | CRM opened 20210314 |

Total Alkalinity Output

| SampleID | TA          | Mass   | Salinity |
|----------|-------------|--------|----------|
| JUNK1    | 1981.003448 | 60.113 | 35       |
| CRM1     | 2182.69197  | 60.225 | 33.417   |
| CRM2     | 2183.070705 | 59.774 | 33.417   |
| CRM3     | 2179.989317 | 60.107 | 33.417   |
| CRM4     | 2180.118916 | 60.161 | 33.417   |
| CRM5     | 2194.262889 | 59.659 | 33.417   |
| CRM6     | 2194.808799 | 59.719 | 33.417   |
| CRM7     | 2192.134088 | 59.651 | 33.417   |
| CRM8     | 2192.253611 | 59.892 | 33.417   |
| CRM1     | 2205.731211 | 60.297 | 33.623   |
| CRM2     | 2194.005521 | 60.319 | 33.623   |
| CRM3     | 2192.044753 | 59.881 | 33.623   |
| CRM4     | 2196.378938 | 59.971 | 33.623   |


#### **Data and notes for testing variability in CRMs with new titrant 20210328**

**Notes**
- Continued the troubleshooting process and changed out the older titrant with a new Batch #A16 titrant. Purged the lines four times to make sure to run through any of the older titrant from the lines. Ran a pH calibration that had sound numbers and then I ran a junk sample before five CRMs. IT WORKED! The CRM values were very accurate and everything looked good from the run. Can continue with samples moving forward.

**Data**

CRM Accuracy Data

| Date     | CRM value   | Batch value | % off        | Batch # | Notes               |
|----------|-------------|-------------|--------------|---------|---------------------|
| 20210328 | 2766.991908 | 2225.47     | 24.33292331  | 180     | CRM opened 20210314 |
| 20210328 | 2226.220045 | 2226.47     | -0.011226509 | 180     | CRM opened 20210314 |
| 20210328 | 2224.425627 | 2227.47     | -0.136674041 | 180     | CRM opened 20210314 |
| 20210328 | 2224.201971 | 2228.47     | -0.191522852 | 180     | CRM opened 20210314 |
| 20210328 | 2225.549637 | 2229.47     | -0.175842846 | 180     | CRM opened 20210314 |
| 20210328 | 2228.292488 | 2230.47     | -0.097625722 | 180     | CRM opened 20210314 |


Total Alkalinity Output

| SampleID | TA               | Mass   | Salinity |
|----------|------------------|--------|----------|
| JUNK1    | 2766.9919084924  | 60.137 | 35       |
| CRM1     | 2226.22004515256 | 60.18  | 33.623   |
| CRM2     | 2224.4256267384  | 60.274 | 33.623   |
| CRM3     | 2224.20197069824 | 60.29  | 33.623   |
| CRM4     | 2225.54963650217 | 60.077 | 33.623   |
| CRM5     | 2228.2924875591  | 60.135 | 33.623   |

#### **Data and notes for testing E5 samples TA signal 20210504**

**Notes**
- The CRM values were very accurate and everything looked good from the previous run, testing a batch of samples from March (TP2) run 1 on 20200304 to see if there is a processable deltaTA signal from the initial sample to the blank and other coral sample chambers. All went well and we were able to see > 30 for deltaTA in almost all samples which is a good value moving forward. Going to run more samples from another run to make sure this occurs in other replicates.

**Data**

CRM Accuracy Data

| Date     | CRM value | Batch value | % off    | Batch # | Notes               |
|----------|-----------|-------------|----------|---------|---------------------|
| 20210504 | 2224.823  | 2224.47     | 0.015869 | 180     | CRM opened 20210418 |


Total Alkalinity Output

| SampleID  | TA               | Mass   | Salinity |
|-----------|------------------|--------|----------|
| JUNK1     | 2030.80445223374 | 60.035 | 35       |
| INITIAL-2 | 2325.79665116253 | 59.525 | 36.03    |
| ACR-225   | 2229.92073279084 | 60.488 | 36.03    |
| POC-215   | 2271.34957005683 | 59.695 | 36.02    |
| BLANK-1   | 2298.63819200678 | 60.297 | 36.03    |
| POR-224   | 2206.22121318086 | 60.266 | 36       |

#### **Data and notes for testing E5 samples TA signal 20210506**

**Notes**
- The CRM values were very accurate and everything looked good from the previous run, testing a batch of samples from March (TP2) run 2 on 20200304 to see if the deltaTA values are also representative of what we saw in yesterdays runs. Most of the samples again had deltaTA's > 30 and all looked good. 

**Data**

CRM Accuracy Data

| Date     | CRM value | Batch value | % off    | Batch # | Notes               |
|----------|-----------|-------------|----------|---------|---------------------|
| 20210506 | 2225.094  | 2224.47     | 0.028052 | 180     | CRM opened 20210418 |


Total Alkalinity Output

| SampleID  | TA               | Mass   | Salinity |
|-----------|------------------|--------|----------|
| JUNK1     | 2018.39880477983 | 60.226 | 35       |
| INITIAL-1 | 2309.66768581353 | 60.083 | 37.38    |
| INITIAL-2 | 2310.6341993138  | 59.885 | 37.38    |
| BLANK-2   | 2312.66284694853 | 59.917 | 37.4     |
| POC-201   | 2302.40007114124 | 59.624 | 37.4     |
| ACR-244   | 2277.31062313817 | 59.786 | 37.39    |
| POR-235   | 2288.18858541504 | 59.525 | 37.41    |
| ACR-229   | 2272.13489556652 | 59.84  | 37.38    |
| POR-245   | 2268.37232119263 | 60.317 | 37.43    |
