---
layout: post
title:  Histological Image Analysis
category: [histology, gametogenesis, timeseries]
tag: [ histological, imaging, QuPath software]
---

# Histological image analysis for *Acropora pulchra* decalcified tissues

#### Goal:

Use QuPath software and ImageJ to observe and analyze all histological images for *A. pulchra* gametogenesis 12-month time series data collected in Moorea, French Polynesia monthly from December 2021 to November 2022. Calculating abundance and size of gonads (spermaries and oocytes) throughout spermatogenesis and oogenesis. 

#### Previous steps and metadata for samples:

Details on project overview, original experimental design, and detailed collection protocols can be found [here](https://github.com/urol-e5/urol-e5.github.io/blob/master/_posts/2022-03-02-March-April-Fieldwork-Overview.md). After collection each month, all samples were processed following the [same day sample processing protocol](https://github.com/daniellembecker/Gametogenesis/blob/main/protocols/2021-12-26-Sample_Same_Day_Processing_Protocol.md) and then decalcified following this [protocol](https://github.com/daniellembecker/Gametogenesis/blob/main/protocols/2022-04-16-Histological-Processing.md).

#### Resources:

[QuPath Software](https://qupath.github.io)
[ImageJ Software](https://imagej.net/ij/download.html)

----------------
### QuPath Image Protocol Steps and Link to ImageJ

1. Download QuPath and ImageJ software for your desktop from [QuPath Link](https://qupath.github.io) and [ImageJ Link](https://imagej.net/ij/download.html).

2. Open QuPath on your desktop and navigate to the project tab in the QuPath viewer.

3. Make a folder on your desktop where your QuPath project will be saved.

4. For large .vsi images, you will need to have your hard drive with your images plugged into your computer whenever you are accessing the project and images moving forward. QuPath needs to know the path to the original files whenever you are using the software.

5. Click on create project in the toolbar under the project tab and navigate to your saved folder on your desktop, click it, and select the open button to save. 

6. QuPath will automatically put folders labeled: classifiers, data, project.qpproj, and projectqpproj.backup into this folder for future use. Anytime you close out of QuPath, you can use the project.qpproj file to look at your images again, as long as your hard drive is plugged in. 

7. Click the add images button also within the project tab and a box will open where you can drag and drop any selected .vsi files from your hardrive into the project. Select the import button and you will then see your images show up in the left side of the program with labels underneath the tabbed section of your project.

8. Double-click any image to open it in the viewer section of QuPath and select what type of image this is, for our images they are Brightfield H&E due to the imaging process at Brown Bioimaging Center and the type of staining.

9. To zoom in and out on the image, use your shift and + key to zoom in and your - key to zoom out.

10. If you want to send a region of the image to ImageJ for measuring, you can use the ImageJ tool that is located if you click the >> button and press the send region to ImageJ option, set the resolution to 0.5 and the resolution unit to um. This will send your image to ImageJ for further processing.

11. To save a .vsi file as a .jpg with high resolution for uploading to an online database like OSF or for sending a copy, press the File selection in the QuPath program, not from the toolbar at the top of your screen. Select Export images > Original pixels which will bring up a Export image region box. Select JPEG from the Export format section, change the Downsample factor to the highest resolution you can. As you change this number, an error will show: Output image size XX by XX pixels (too big!) in red. As you keep changing it to a lower resolution (go up by 0.1 decimal point at a time), it will eventually show: Output image size XX by XX pixels with no warning after and it will be grey. Then click OK and save to an output folder.

12. Once you are done, you can save any changes to your images in the project. If you want to open your project again, plug in your hard drive with the raw .vsi files, open QuPath, navigate to the folder on your desktop where you have stored your project.qpproj file, select it and drag it into the project section of QuPath.






