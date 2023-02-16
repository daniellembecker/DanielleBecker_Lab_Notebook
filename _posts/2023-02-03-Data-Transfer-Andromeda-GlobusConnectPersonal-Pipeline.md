---
layout: post
title: Workflow for Transferring Genome Feature Files from Andromeda to GlobusConnectPersonal endpoint
date: 2023-02-03
category: [ Protocol ]
tag: [ Data Transfer, RNASeq, WGBS ]
projects: Molecular Underpinnings, Chronic Low-Level Nutrient Enrichment
---

### Objective:

#### Develop a streamlined way to transfer large genome feature files (e.g., BAMs, bisulfite-converted genomes, etc) from [Andromeda](https://its.uri.edu/research-computing/using-andromeda/) (URI HPC server) to [Globus Connect Personal](https://app.globus.org/file-manager) desktop endpoint and [Globus Connect Server](https://docs.globus.org/globus-connect-server/v5.4/) online databases.

##### To develop this protocol, I followed the general methods and protocol for this approach by the [Roberts Lab](https://robertslab.github.io/resources/code_Snippets/#transfer-files-tofrom-mox-using-globus-connect-personal). Their HPC server is called Mox and they have used this to transfer large files in the past.

#### Important Note:
##### To be able to transfer data files from the Andromeda HPC server for free, you need to have [Globus Connect Server](https://docs.globus.org/globus-connect-server/v5.4/) on Andromeda and [Globus Connect Personal](https://docs.globus.org/how-to/globus-connect-personal-windows/) on your designated desktop. If you use [Globus Connect Personal](https://docs.globus.org/how-to/globus-connect-personal-windows/) on the Andromeda server you will not be able to transfer for free. Unfortunately, Andromeda does not have a public IP from which to run the Globus Server. Looking into options for a subscription to Globus Connect Personal.


### Workflow for connecting GlobusConnectPersonal on the Andromeda server:

#### Step 1: Log into the URI Andromeda server
- ssh Andromeda login here
- instead of the usual login '**3' you need to use '**4' so that the transfers do not affect anything interactively on **3

#### Step 2: Setup Globus collection
- Download Globus Connect Personal onto the putnamlab server folder:
    - /data/putnamlab/
    - $ wget https://downloads.globus.org/globus-connect-personal/linux/stable/globusconnectpersonal-latest.tgz
- Extract the files from the downloaded tarball.
    - $ tar xzf globusconnectpersonal-latest.tgz
    - this will produce a versioned globusconnectpersonal directory
    - replace `x.y.z` in the line below with the version number you see
    - $ cd globusconnectpersonal-x.y.z
- Start Globus Connect Personal. Since this is the first time you are running it, you must complete setup before you can run the full application.
    - $ ./globusconnectpersonal

#### Step 3: Enter Globus Connect Personal login and authentication code information
  - After entering '$ ./globusconnectpersonal', your terminal will say:
      - Will now attempt to run globusconnectpersonal -setup
      - Globus Connect Personal needs you to log in to continue the setup process.
      - We will display a login URL. Copy it into any browser and log in to get a single-use code. Return to this command with the code to continue setup.
  - A login URL will appear starting with: https://auth.globus.org
  - Copy this URL in any web browser.
  - Login to your Globus Connect Personal using your ORC-ID, SSH for URI, or your email.
  - Click allow from the options seen below.

  ![globus](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/images/globusconnect.png)

  - You will then be directed to your Native App Authorization Code which is ~30 characters long.
  - Copy and paste this into your terminal window where the prompt says enter the auth code.

#### Step 4: Add desired Andromeda directory to config file and set permissions, example for putnamlab:

  - When the prompt in terminal asks: Input a value for the Endpoint Name:
  - Enter: ~/data/putnamlab/
  - You will then need to give the folder permissions.
    - Go into the global connect folder on Andromeda:
      - cd globusconnectpersonal-3.2.o
    - Make your home directory readable/writeable by Globus by using the command below.
      - echo  /data/putnamlab/,0,1 >> ~/.globusonline/lta/config-paths
      - ~/,0,1: Makes your home directory readable/writeable by Globus.
      - Makes the putnam lab directory on /data/putnamlab/ readable/writeable by Globus.
  - Start Globus Connect Personal: ./globusconnectpersonal -start. Nothing will happen after you hit enter. The cursor will simply flash - this is good.
  - Login to your Globus Connect Personal account via a web browser.
  - Click on Collections and you should now see your collection (name provided in endpoint Step 5), and it should have a green stack of papers next to it; the green indicates that the connection is activate.
  - You must be running " ./globusconnectpersonal -start" for the paper stack to be green. If it's red or yellow, you will need to re-run it. I recommend only running while you have an active transfer.
  - Click on the collection name.
  - Click on "Open in File Manager" (on the right side of the screen).
  - Navigate to the directory you setup in Step 5. NOTE: You'll have to navigate up a directory out of your home directory in order to get to the /putnamlab partition.
  - Transfer data from other Globus Endpoint to Andromeda and from Andromeda to the Globus endpoint!



#### Step 5: If you need to ever delete the globusconnectpersonal environment to start over with new auth code:

  - Go to the directory above your globusconnectpersonal-3.2.0 folder and globusconnectpersonal-latest.tgz file.
  - Use command rm globusconnectpersonal-latest.tgz and rm -r globusconnectpersonal-3.2.0 to delete directory and file.
  - Restart from Step 3.

### References:

[Roberts Lab Protocol](https://robertslab.github.io/resources/code_Snippets/#transfer-files-tofrom-mox-using-globus-connect-personal)

[Globus Connect Personal Protocol](https://docs.globus.org/how-to/globus-connect-personal-linux/)

[Globus Connect Server Youtube](https://www.youtube.com/watch?v=8ILtsSRiML8&feature=youtu.be)

[Globus Connect Server Protocol](https://docs.globus.org/globus-connect-server/v5.4/)
