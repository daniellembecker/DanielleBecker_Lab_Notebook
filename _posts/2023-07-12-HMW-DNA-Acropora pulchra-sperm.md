---
layout: post
title: High Molecular Weight DNA for Acropora pulchra sperm
category: [ HMW , DNA]
tag: [ Acropora pulchra, HMW ]
---
## Testing how to extract High Molecular Weight DNA for Acropora pulchra sperm

#### Goal:
Use *Acropora pulchra* sperm samples to extract high molecular weight DNA for E5 project *A. pulchra* genome. Testing methods between the Genomic Tip High Molecular Weight DNA Extraction Protocol and NEB High Molecular Weight DNA Extraction.

#### Example Protocols

[Maggie Schedl Qiagen Genomic Tip HMW DNA Extraction Protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/HMW-Tip-Protocol/)

[Maggie Schedl Qiagen Genomic Tip HMW DNA Extraction of Fire Coral](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/firecoral-HMW/)

[Maggie Schedl HMW DNA Extraction of Pocillopora acuta with the Qiagen Genomic Tip](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/pacuta-HMW/)

[New England Biolabs High Molecular Weight DNA Extraction Guidelines](https://www.neb.com/monarch/high-molecular-weight-dna-extraction?gclid=Cj0KCQjwnrmlBhDHARIsADJ5b_kHBvc274KwNljTVwN6r2KJbjBQjnWFzkyJLI6GpmEZw7UNU-HPjmEaAt_tEALw_wcB)

[New England Biolabs High Molecular Weight DNA Cells and Blood Kit](https://www.neb.com/-/media/nebus/files/manuals/manualt3050.pdf?rev=41e2c417f75c4889b2d0caeaa6746419&hash=1A3E03456DCB79B426622B85DFB85111)


**Important Notes about HMW DNA** from [Maggie](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/HMW-Tip-Protocol/
)

- HMW DNA is very fragile, pipetting or vortexing can sheer it easily
- Do not vortex after extraction, to mix the samples flick gently then spin down
- Whenever you are pipetting the actual sample (ie, not buffers) you should use wide bore pipette tips, these are next to the kit on the Putnam shelf. They have a larger opening to the tip and are gentler
- **Do not freeze the DNA**, freeze thaws can break up the long fragments, keep HMW DNA in the 4 degree fridge for storage
- HMW DNA can really clump, especially in this extraction method that pellets it in a precipitation. For quantifying with the Qubit, take a top and bottom of the tube 1ul: do 2 Qubits for each sample and average them together for the whole tube concentration

#### Materials

- Monarch® HMW DNA Extraction Kit for Cells & Blood
NEB #T3050S/L
- Thermal mixer containing a 2 ml tube block (if not available, use a 1.5 ml block).
- Isopropanol, 275 μl per sample (Low Input: 100 μl/sample).
- Ethanol (≥ 95%)
- 1.5 ml DNase-free, low DNA binding microfuge tubes (e.g., Eppendorf® DNA LoBind®, #0030108051) are recommended for elution and storage (1 per sample)
- Recommended: vertical rotating mixer (e.g., Thermo Scientific® HulaMixer® Sample Mixer).
- Wide-bore pipette tips.

#### Protocol, Using the [**New England Biolabs/Monarch High Molecular Weight DNA Cells and Blood Kit**](https://www.neb.com/-/media/nebus/files/manuals/manualt3050.pdf?rev=41e2c417f75c4889b2d0caeaa6746419&hash=1A3E03456DCB79B426622B85DFB85111)

MATERIALS REQUIRED BUT NOT SUPPLIED

- Thermal mixer containing a 2 ml tube block (if not available, use a 1.5 ml block).
- Isopropanol, 275 μl per sample (Low Input: 100 μl/sample).
- Ethanol (≥ 95%)
- 1.5 ml DNase-free, low DNA binding microfuge tubes (e.g., Eppendorf® DNA LoBind®, #0030108051) are recommended for elution and storage (1 per sample)
- Recommended: vertical rotating mixer (e.g., Thermo Scientific® HulaMixer® Sample Mixer).
- Wide-bore pipette tips.

IMPORTANT NOTES BEFORE YOU BEGIN
- Review the complete protocol before beginning.
- Add ethanol (≥ 95%) to the gDNA Wash Buffer as indicated on the bottle label.
- Cool the Nuclei Prep Buffer to 4°C.
- Preheat thermal mixer with 2 ml block to 56°C.

Protocol completed on Wednesday, July 19th 2023 by DMBP and ZD

1. We used *A. pulchra* sperm sample (tube number: 207, ambient treatment) preserved in Zymo DNA/RNA shield from October 2022 spawning in Mo'ore'a, French Polynesia. These samples were preserved by mixing 1 mL of sperm with 1 mL of DNA/RNA shield. The concentration of sperm used was 10^6 sperm/mL, so this 2mL volume of sperm + DNA/RNA shield should contain 10^6 sperm, or 500,000 cells/mL. The recommended input for this kit is 10^6 cells, so one full tube of sperm + DNA/RNA shield should be roughly the desired input. We need to figure out how to reduce the volume of this for the input.

2. According to the Monarch protocol, the cells (sperm) stored in DNA/RNA shield are already sufficiently lysed:

    > "For samples stored in 2X Monarch DNA/RNA Protection Reagent, cells are already lysed. Dilute the reagent to 1X with nuclease-free water and use this in place of the Nuclei Prep and Lysis Buffer combination for sample lysis. Following the 10-minute incubation with Proteinase K, add RNase A and incubate the sample at 56°C for 5 minutes. Subsequently, follow the standard protocol."

    Therefore, we will thaw the sperm in the DNA/RNA shield, and vortex to homogenize/lyse the sperm. Then, since the DNA/RNA shield is already at 1X concentration, we will use this mixture as input for the following steps.

3. Transfer **300 uL of the sperm** in DNA/RNA shield volume to the  Monarch 2 ml Tube. Add 5uL of RNase A and 10 uL of Proteinase K to this. Invert 10 times to mix. Avoid introducing air bubbles. Do not vortex or pipette.

    > This will be approximately 15,000 or 1.5 x 10^4 cells, which is below the recommended input of cells. We may need to increase the input to 600 uL and increase the RNAse A and Proteinase K accordingly (as well as the Precipitation Enhancer and Isopropanol below).

4. Incubate at **56°C for 10 minutes** in a thermal mixer with agitation at the desired speed to control the shearing and tune the size of gDNA. The speed of the thermal mixer influences fragment length; higher speeds reduce overall size. For the standard ligation based Oxford Nanopore Technologies (ONT) sequencing protocol, agitation at **2,000 rpm is recommended**. If desired, samples can be stored at 4°C overnight after the incubation.

Note: at this step we set the agitation speed at 1,500 rcf because 2,000 rpm seemed too fast.

5. Prepare and label the plastics for Part 2: HMW gDNA Binding and Elution.
Each sample will require:
    - 1 Monarch Collection Tube II (no need to label).
    - 1 Monarch Bead Retainer inserted into the collection tube; this will be used to remove the wash buffer from the gDNA bound to the beads.
    - 1 Monarch 2 ml Tube; this will be used for eluting the gDNA from the beads.
    - 1 1.5 ml microfuge tube (DNA low bind is recommended, not provided); this will be used to collect the eluate.

6. Add **75 µl of Precipitation Enhancer** after the 10-minute incubation and mix by inverting 8–10 times.

Note: at this step, you should see precipitate in the tube

![image1](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/images/HMW_precip.jpeg)

##### Protocol Part 2: HMW gDNA BINDING AND ELUTION

1. Using clean forceps, **add 2 DNA Capture Beads** to each sample, which should be contained in a Monarch 2 ml Tube.

2. Add **275 μl isopropanol**, close the cap, and mix on a vertical rotating mixer at 10 rpm for 4 minutes to attach DNA to the beads. When working with ≥ 5 x 10^6 cells, double the inversion time to 8 minutes; this is especially important if low agitation speeds were used during lysis.

    If a vertical rotating mixer is not available, invert slowly and gently by hand 30 times. A manual inversion is complete when the tube returns to the upright position. Slow inversion is critical for the DNA to bind to the beads; each full inversion should take ~5–6 seconds. If necessary, flick the tube to release any beads that stick to the bottom of the tube.

    After a few inversions, the solution becomes more viscous and the DNA will wrap loosely around the beads. During the following inversions, precipitation of gDNA may be visible, especially with sample inputs ≥ 5 x 10^6 cells. The DNA complex will often contain small air bubbles. With more inversions, the DNA will completely wrap around the beads, often causing the beads to stick together.

Note: at this step should see DNA precipitate on the DNA beads

![image2](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/images/HMW_dna_ball.jpeg)


3. **Remove and discard liquid by pipetting**. Avoid removing any of the gDNA wrapped around the glass beads. For optimal DNA solubility, avoid letting the bound DNA dry out on the beads during this and the following steps; **add the next buffer quickly**. There are two suggested options for carrying out this step:
    - Keeping tube upright, insert pipette tip and gently push beads aside to remove liquid.
    - Angle tube so that beads remain at the bottom, and liquid reaches toward tube opening. Pipette from the liquid surface and continue to angle as liquid is removed (tube will be almost horizontal at the end).

<img src="https://github.com/zdellaert/ZD_Putnam_Lab_Notebook/blob/master/images/protocols/Monarch_HMW_1.png?raw=true"  width="400" />

4. Add **500 μl gDNA Wash Buffer**, close the cap, and mix by inverting the tube 2–3 times. Remove the gDNA Wash Buffer as described in Step 3. The loose gDNA complex will condense around the beads more tightly.

5. Repeat the wash in Steps 3–4. Remove the gDNA Wash Buffer by pipetting. It is not necessary to remove all the gDNA Wash Buffer at this point.

6. Place a labeled bead retainer into a Monarch Collection Tube II. Pour the beads into the bead retainer and close the cap. Discard the used Monarch 2 ml Tube. When working with multiple samples, be sure to close the cap of the bead retainer after each transfer of beads.

<img src="https://github.com/zdellaert/ZD_Putnam_Lab_Notebook/blob/master/images/protocols/Monarch_HMW_2.png?raw=true"  width="200" />

7. Pulse spin (1 second or less) the sample in a benchtop minicentrifuge to remove any residual wash buffer from the beads.

8. Separate the bead retainer from the collection tube, pour the beads into a new, labeled Monarch 2 ml Tube, and **insert the used bead retainer into the labeled 1.5 ml microfuge tube** (_DNA low bind recommended_, not provided) for later use during elution. Discard the used collection tube.

<img src="https://github.com/zdellaert/ZD_Putnam_Lab_Notebook/blob/master/images/protocols/Monarch_HMW_3.png?raw=true"  width="200" />

9. Immediately add **100 μl** (use 200 μl if working with > 2 x 10^6 cells) **Elution Buffer II** onto the glass beads and **incubate for 5 minutes at 56°C in a thermal mixer with agitation at the lowest speed (300 rpm)**. Halfway through the incubation, ensure the beads are not stuck to the bottom of the tube by *tilting the tube almost horizontally and gently shaking*.

    This ensures that the beads can move freely, allowing for optimal release of the DNA from the beads. It also ensures that the lower bead does not stick to the bottom of the tube during the following transfer step.

<img src="https://github.com/zdellaert/ZD_Putnam_Lab_Notebook/blob/master/images/protocols/Monarch_HMW_4.png?raw=true"  width="200" />

10. Ensure the bead retainer is inserted into the 1.5 ml microfuge tube. **Pour the eluate and the glass beads into the bead retainer and close the cap**. When working with more than 1 sample, it is important to close the cap after each transfer of beads. Typically, all the eluate flows into the bead retainer upon pouring. If any volume remains in the 2 ml tube, spin briefly and transfer.

<img src="https://github.com/zdellaert/ZD_Putnam_Lab_Notebook/blob/master/images/protocols/Monarch_HMW_5.png?raw=true"  width="200" />

11. **Centrifuge for 30 seconds** (1 minute if working with ≥ 5 x 10^6 cells) at 12,000 x g to separate the eluate from the glass beads. *Discard the beads and retainer.*


12. **Pipette eluate up and down 5–10 times with a wide bore pipette tip and ensure any visible DNA aggregates are dispersed.** Before analysis or downstream use, HMW DNA must be homogeneously dissolved. **After pipetting, incubate at 37°C for 30-60 minutes, overnight at room temperature, or for > 24 hours at 4°C**.

    *Pipette up and down 5-10 times again before analyzing or using the HMW DNA*. Samples processed using low agitation speeds during lysis will require additional time to fully dissolve. See additional guidance in “Homogenization of HMW DNA Samples”, page 24. Samples can be stored at 4°C for short term use (weeks), or at -20°C for long term storage. The elution buffer (10 mM Tris, pH 9.0, 0.5 mM EDTA) is formulated for long term storage of gDNA.

Note: at this step should see DNA precipitate at bottom of tube

![image3](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/images/HMW_end_dna.jpeg)

## QC on samples 20230720 DMBP and ZD

#### Qubit

- Use the broad range DNA assay and [protocol](https://github.com/meschedl/PPP-Lab-Resources/tree/master/Protocols_and_Lab_Resources/DNA_Quality_Control)
- Quantify both the top and bottom of each tube and then use the averages of each tube for the final concentration
- You should be getting very high DNA quants for this, above 100ng/ul is good

Qubit results:

| Sample.ID | Qubit Reading (ng/ul) |
|-----------|---------------|
| S1        | 214.67        |
| S2        | 24040.78      |
| Top1      | 13.9          |
| Top2      | 15.3          |
| Bottom1   | 45.6          |
| Bottom2   | 64.8          |

Notes: These data show that there is a clear separation in the sample tube and is not fully homogenized, we need to better disperse the clumps of DNA without too much disruption.

#### TapeStation

- For accurate sizing, we use the TapeStation instead of a gel. However it still isn't the absolute best, the best thing for HMW DNA is a pulse-phase gel that's run super slowly and alternates the position of the current to separate out the large fragments. You can look at [pictures](https://international.neb.com/monarch/high-molecular-weight-dna-extraction), pretty cool. However, we don't have the type of gel box for that!
- Use the **Genomic DNA** screentapes and reagent/ladder. This reads from 48,000bp to 100bp
- Follow the [protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/DNA-Tapestation/)
- Your goal is to have a huge peak at the 60,000bp size or around there, and very minimal smearing elsewhere

Overall tapestation results:

![image4](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/images/HMW_tape.jpeg)

Bottom sample:

DIN = 8.2 Concentration = 13.3

![image5](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/images/HMW_bottom.jpeg)

Top sample:

DIN = 8.1 Concentration = 14.5

![image6](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/images/HMW_tape_top.jpeg)

Notes: Traces of DNA are good for the tapestation with good gDNA band @ >60,000bp.

Next steps: Since we wanted to homogenize the DNA more if possible, Zoe placed the sample on a stir plate for an hour and then back into the fridge overnight at 4°C.

Qubit results on 20230721 by ZD:

| Sample.ID | Qubit Reading (ng/ul) |
|-----------|---------------|
| S1        | 221.02        |
| S2        | 25790.01      |
| Top1      | 18.2          |
| Top2      | 17.7          |
| Bottom1   | 19.6         |
| Bottom2   | 19.3         |

Notes: These data show that the samples is still not homogenized, samples placed back in 4°C fridge.

Next steps: Danielle will QC the sample again on Wednesday the 26th to see if anymore homogenizing has happened, if not, she is going to try another sample without the change to the thermo mixer speed. Could also increase the input amount of sample from 300 uL to 600 uL.

Qubit results on sample completed on 20230721, read on 20230727:

| Sample.ID | Qubit Reading (ng/ul) |
|-----------|---------------|
| S1        | 193.03        |
| S2        | 21519.79      |
| Top1      | 19.5         |
| Top2      | 17.6         |
| Bottom1   | 262        |
| Bottom2   | 312        |

Next steps: On 20230810 Danielle increased the input amount of sample from 300 uL to 600 uL on 20230810 for HMW. Adjusted ProteinaseK and RNAse volume additions by x2 to match new volume input of sample.

Qubit results on 600 uL sample:

| Sample.ID | Qubit Reading (ng/ul) |
|-----------|---------------|
| S1        | 189.38       |
| S2        | 19206.16      |
| Top1      | 19.6         |
| Top2      | 19.7         |
| Bottom1   | 20.8       |
| Bottom2   | 21.4        |

Next steps: On 20231012 Danielle and Hollie looked into troubleshooting steps. It was suggested to reduce the material (used 150 uL of sperm sample) from the original input volume of 300 uL. Upon using half of the input volume, step 3 of the protocol was adjusted to use half the ProteinaseK and RNase as well. Also, Danielle used the vertical rotating mixer now set up by the molecular Putnam bench so that this step will be more precise and hopefully help to dissolve/re-suspend the HMW DNA. The last adjustment was to reduce the elution buffer by 1/2 as well for step 9 of the protocol.

Qubit results on 150 uL sample:

| Sample.ID | Qubit Reading (ng/ul) |
|-----------|---------------|
| S1        | 208.21       |
| S2        | 21991.61     |
| Top1      | 18.38        |
| Top2      | 18.28         |
| Bottom1   | 18.70       |
| Bottom2   | 18.60       |

Next steps: Was waiting to confirm which approach and what quality and quantity was required for HI-C approach or Illumina. On 20231027 I used the original input volume of 300 uL. I also used the vertical rotating mixer now set up by the molecular Putnam bench so that this step will be more precise and hopefully help to dissolve/re-suspend the HMW DNA. The last adjustment was to reduce the elution buffer by 1/2 as well for step 9 of the protocol.

Qubit results on 300 uL sample 20231027:

| Sample.ID | Qubit Reading (ng/ul) |
|-----------|---------------|
| S1        | 209.32      |
| S2        | 21501.81     |
| Top1      | 30.08       |
| Top2      | 29.98         |
| Bottom1   | 30.12       |
| Bottom2   | 30.05       |

Qubit results on 300 uL sample 20231102:

| Sample.ID | Qubit Reading (ng/ul) |
|-----------|---------------|
| S1        | 196.68      |
| S2        | 21113.09     |
| Top1      | 39.6      |
| Top2      | 39.4         |
| Bottom1   | 39.2       |
| Bottom2   | 39.6       |

**NanoDrop**

To check the purity of the DNA sample, I used a NanoDrop on 20231103:

- I followed the [NanoDrop Protocol](https://github.com/Putnam-Lab/Lab_Management/blob/master/Lab_Resources/DNA_RNA-protocols/NanoDrop-2000c-Protocol.md) exactly using the NEB Elution Buffer II as my blank and the NanoDrop2000c in Dr. Jenkins lab CBLS 260.


Values:

| Sample ID    | Qubit concentration | NanoDrop Concentration | 260/280 | 260/230 |
|--------------|---------------------|------------------------|---------|---------|
| HMW_20231027 | 39.5 ng/uL          | 57.1 ng/uL             | 1.74    | 1.24    |

Composite Figure:

![firstnanodrop](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/images/NanoDrop/firstnanodrop.jpg)

We realized that the 260/280 and 260/230 ratios were a bit low for purity for the recommendations for submitting a sample for library prep and Illumina sequencing at Quick Biology Inc.
- Recommendations:
    - Illumina Need = DNA amount: ≥ 200 ng for normal input (PCR-free); minimum 50 ng (PCR option); DNA concentration: ≥20 ng/μL, Optimal concentration: 20 ng/μl to 500 ng/μl in Tris-HCL (pH8.0); Purity: A260/280:1.8-2.0, A260/230: ≥ 1.8, Intact DNA, RNA-free, and non-degraded.

Zoe suggested using an Ethanol Precipitation Protocol to clean up the sample and improve the purity. So, we followed the [Ethanol Precipitation Protocol](https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/PCR-Protocol/) section to clen up the eluted DNA.

On 20231103 we cleaned the sample.

Protocol steps:
1. Add 9.6 uL of 3M sodium acetate to the eluted DNA suspended in ~96 uL of NEB elution buffer II. (1/10th the total volume).
2. Add 316.8 uL of **ice cold** 100% ethanol to each tube (3x the total volume after adding sodium acetate). Equation = 96 uL eluted sample + 9.6 uL of 3M sodium acetate = 105.6 uL x 3 = 316.8 uL of ethanol.
3. Pipette up and down to mix and transfer entire volume to a new labeled 1.5ml tube.
4. Place tubes in -80 ºC freezer for 1 hour or in the -20 ºC freezer overnight.
5. Spin for 30 minutes at 15,000rcf (room temp). Make sure all the tubes are facing the same way (hinges out). The DNA will pellet on the back wall.
6. Pipette off the liquid without disturbing the pellet. It will most likely be invisible.
7. Add 200µL of **ice cold** 70% Ethanol.
8. Spin for 15 minutes at 15,000rcf. This is to wash the pellet. Again make sure each tube is facing the same direction.
9. Pipette off the ethanol without disturbing the pellet.
10. Let dry at room temperature with the caps open for about 15 minutes or until all ethanol has evaporated (put in fume hood).
11. Add 30µL of Tris to the pellet and pipette up and down to resuspend, go straight to nanodrop/qubit quantification or store cleaned PCR product in the -20 ºC freezer.

After cleaning the sample, we checked for quantity, quality, and purity of the sample and ran another Quibit, a NanoDrop reading, and a tapestation.

**Qubit**

We used the broad range DNA assay and [protocol](https://github.com/meschedl/PPP-Lab-Resources/tree/master/Protocols_and_Lab_Resources/DNA_Quality_Control) like above.

Qubit results after cleaning sample on 20231103:

| Sample.ID | Qubit Reading (ng/ul) |
|-----------|---------------|
| S1        | 206.31      |
| S2        | 21891.51     |
| HMW1     | 69.5     |
| HMW2     | 70.1         |

**NanoDrop**

Followed the same [NanoDrop Protocol](https://github.com/Putnam-Lab/Lab_Management/blob/master/Lab_Resources/DNA_RNA-protocols/NanoDrop-2000c-Protocol.md) but with Tris as my buffer.

Values:

| Sample ID    | Qubit concentration | NanoDrop Concentration | 260/280 | 260/230 |
|--------------|---------------------|------------------------|---------|---------|
| HMW_20231027 | 70.0 ng/uL          | 101.9 ng/uL             | 1.80    | 1.55    |

Composite Figure:
![aftercleaningnanodrop](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/images/NanoDrop/aftercleaningnanodrop.png)

The cleaning of the eluted DNA improved by 0.06 for the 260/280 ratio and by 0.31 for the 260/230 ratio! Yay!

I checked with Quick Biology and they said these results were acceptable to move forward with library prep.

**Tapestation**

To be sure that we still had HMW DNA in our sample after these extended steps, we followed the [TapeStation Protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/DNA-Tapestation/) exactly.

![tapestation](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/images/NanoDrop/tapestation_image_1105.PNG)

The quality looks very good. And we see peaks are in the >60,000bp which is what we want for HMW DNA :) Yay! Sample is ready to go and be shipped to Quick Biology Inc.
