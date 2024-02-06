---
layout: post
title: Workflow for using Unity remote desktop through Research Computing at URI to analyze histology images
date: 2024-02-06
category: [ Protocol ]
tag: [ remote desktop, Unity, Histology analysis]
projects: Gametogenesis, Acropora pulchra
---

### Objective:

Develop a way to use URIs Research Computing [Unity Remote Desktop](https://ood.unity.rc.umass.edu/pun/sys/dashboard/batch_connect/sys/bc_desktop/session_contexts/new) for analyzing histological images (.VSI format) from the [Olympus VS200 Slide Scanner](https://www.olympus-lifescience.com/en/solutions-based-systems/vs200/?creative=651058250123&keyword=vs200&matchtype=p&network=g&device=c&campaignid=19808814716&adgroupid=153723253824&gad_source=1&gclid=Cj0KCQiAzoeuBhDqARIsAMdH14G2xPd8sI9oVtAmzihQT06OLMU9ZTR4dJPopi_AG-fuMgdWLPsQDoQaAlKkEALw_wcB) in the imaging applications [FIJI](https://fiji.sc) and [ImageJ](https://imagej.net/ij/). FIJI is an additional plug-in to ImageJ to process .VSI images and larger memory files. We currently are having issues on normal computer operating systems because the .VSI files use ~300MB just to open in FIJI, which leads to very slow operating times and errors with RAM usage for the high resolution image options we need. By using the remote desktop, we can expedite the speed of opening and processing the images without having to rely on personal computer storage.

### Resources to look over before use

-  [Unity OpenOnDemand Desktop](https://ood.unity.rc.umass.edu/pun/sys/dashboard/batch_connect/sys/bc_desktop/session_contexts/new)
- [Storing image files on Unity](https://docs.unity.uri.edu/documentation/managing-files/)

### Workflow for opening and processing histology image files through the Unity remote desktop:

##### Step 1: Make sure you are connected to the URI VPN or on URIs campus to access the remote desktop

- [URI IT services virtual private network set-up instructions](https://its.uri.edu/services/94530c3f3f00b35e6d294546faa3667f75f63fff22/)

##### Step 2: Log into the [Unity OpenOnDemand Desktop](https://ood.unity.rc.umass.edu/pun/sys/dashboard/batch_connect/sys/bc_desktop/session_contexts/new) with your URI credentials on your computer desktop

- When you see the Unity Single Sign-On page, click on the University of Rhode Island organization icon

single sign on image

- Log in with your URI email credentials
- If you are using the URI VPN, you will be directly to the [Cisco DUO mobile verification page](https://duo.com/resources/ebooks/the-multi-factor-authentication-evaluation-guide?utm_source=google&utm_medium=paid_search&utm_campaign=DUO_AMER_NA_GS_Branded_General_T1&utm_content=General&gad_source=1&gclid=Cj0KCQiAzoeuBhDqARIsAMdH14HikwqO9eUWvR37ewGtuHWGonw3uZOpZG2uyI_-7Bl8Nip2bleyN2UaAnszEALw_wcB)
- Accept the Cisco DUO mobile notification on your phone application and you will be logged into the [Unity OnDemand Desktop](https://ood.unity.rc.umass.edu/pun/sys/dashboard/batch_connect/sys/bc_desktop/session_contexts/new) on your computer desktop

Unity desktop image

- hf

##### Step 3: Upload your images or files to a Unity directory that you can set up with URI research computing

- On the top left of the Unity OnDemand page, select the Files drop down and select the directory where you would like your files stored
- I navigated to the Putnam work directory called: / work / pi_hputnam_uri_edu /

putnam stoarge image

- I then selected the 'New Directory' option and made a directory to store my *Acropora pulchra* histological images: 'Apul_histo_images'
- Select the new directory you have made and select the 'Upload' option
- From your hard drive or computer, drag in the images or files you would like stored there
- A panel will pop up that says how mmany images it will be uploaded, confirm the images and then it will begin to load them into your new directory.

- directory image





##### Step 3: Launch an interactive desktop environment with the specifications you need for your project

- On the Unity desktop page, you will be able to choose specific parameters for your session. You will be able to choose partitions, maximum job duration, CPU core count, memory (in GB), GPU count, modules, and any extra arguments you would like.
- Descriptions for each of these options are on the Unity desktop homepage and you can change the settings for each new session
- For my image analysis, I used the specific settings I will outline below since I have images that are ~400MB each to open and I was advised from research computing to request 120 GB of memory or more to ensure there was enough RAM to work with when requesting the desktop.
- I was also advised to request 8-16 cores or more for CPU core count to help it run faster since the defaults would be too low

Settings:

unity desktop settings image


- Make sure to select the I would like to receive an email option and press 'Launch'
- You will then be redirected to a page that will show you that your Unity desktop session is queued

queued imagej

- It will take ~1 minute for this to change to running and then you can select 'Launch Unity Desktop'

launch desktop imagej


##### Step 4: Use FIJI/ImageJ on your Unity remote Desktop

- When you begin your session, a new tab will open with the homepage for your remote desktop that will look like a normal computer screen

desktop unity imagej

- Select the black and white terminal icon at the bottom left to open a linux/command line terminal

temrinal operation

- Depending on the modules you need for your project, contact research computing at itrcs-group@uri.edu
- For our FIJI/ImageJ analysis, we needed the updated Fiji/2.14.0-Java-1.8, which research computing updated on February 5th
- To open FIJI on the remote desktop, you need to load the module in terminal
- Type into the terminal: 'module load uri/main Fiji/2.14.0-Java-1.8'
- After this loads, then type in the 'ImageJ-linux64' command to open the applications

ImageJ image

- After this is opened, you are now using the FIJI/ImageJ application as you would with any other computer systems
- To open your stored files from Unity, select the 'File' option in FIJI, then 'Open'. This will show you your personal desktop operating folders.
- Select the hard drive icon at the top left to the right of the pencil icon and you will see a list of directories
- Our putnam lab directory is in the work directory, so I double clicked into the 'work' directory and then double clicked into the 'pi_hputnam_uri_edu' directory to find my uploaded folders


photo directory image 1
photo directory image 2
