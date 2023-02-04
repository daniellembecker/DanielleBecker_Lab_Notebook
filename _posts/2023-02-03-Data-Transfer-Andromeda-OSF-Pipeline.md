---
layout: post
title: Workflow for Transferring Genome Feature Files from Andromeda to OSF.io
date: 2023-02-03
category: [ Protocol ]
tag: [ Data Transfer, RNASeq, WGBS ]
projects: Molecular Underpinnings, Chronic Low-Level Nutrient Enrichment
---

### Objective:
#### Developing a streamlined way to transfer large genome feature files (e.g., BAMs, bisulfite-converted genomes, etc) from [Andromeda](https://its.uri.edu/research-computing/using-andromeda/) (URI HPC server) to [OSF.io](https://osf.io/dashboard) through the [Globus Connect Personal](https://app.globus.org/file-manager) online database.

#### Following the general methods and protocol for this approach by the [Roberts Lab](https://robertslab.github.io/resources/code_Snippets/#transfer-files-tofrom-mox-using-globus-connect-personal). Their HPC server is called Mox and they have used this to transfer large files in the past.

### Workflow:

#### Step 1: Log into the URI Andromeda server
- ssh username@ssh4.hac.uri.edu
- instead of the usual login 'ssh3.hac.uri.edu ' you need to use 'ssh4.hac.uri.edu' so that the transfers do not affect anything interactively on ssh3

#### Step 2: Load the package for the Globus Connect Personal database
- module load GlobusConnectPersonal/3.2.0

#### Step 3: Setup Globus collection
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
