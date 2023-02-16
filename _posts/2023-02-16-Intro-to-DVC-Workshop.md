---
layout: post
title: Intro to DVC Workshop
date: '2023-02-16'
categories: Analysis
tags: DVC
---

Today, I attended a workshop through the UW eScience Data Science postdoc group for an introduction to DVC, a version-control platform. 

**Long story short, it seems like a really useful tool to version control large files (images, sequences, bioinformatics files) integrated with GitHub!**   

# What is DVC? 

The [DVC website](https://dvc.org/doc) has a lot of information about this platform.  

DVC = Data Version Control  

 - Command line tool 
 - Open source
 - Integrates with Git
 - Language platform storage agnositc (not dependent on language, platform, or storage - R, Python, SSH, Drive, all of it is supported!) 
 - Managing experiments
 - Machine learning pipelines 

# How does versioning work in DVC? 

In Git, you are pushing and pulling from remote Git server (GitHub, GitLab) to local machine.  

What happens when you have to track really big files or directories?  

DVC can help with this! It creates a reference file from the original file with a .dvc extension that is a metadata file unique to that version of the dataset. These files are very small. You can then version the DVC metadata file in Git. The original file then goes into a remote storage (SSH, Drive, cloud, etc). 

You then use `dvc push` and `dvc pull`. This gets you the exact versions of data using the metadata files as a key. 

*How does this work with access to proprietary storage (e.g., Google Drive permissions)?*  

- Storage specifics are found on the website 
- You will need to give permissions in folders in Drive that you are sharing the data with, using tokens. 
- The same type of approach will work with SSH. 


# Working with DVC 

Check out the [get started page on DVC](https://dvc.org/doc/start) and [commands help page](https://dvc.org/doc/command-reference).  

1. In a repository, you start with `dvc init`, similar to `git init`. This creates a `.dvc` directory. 
2. Then we need to let DVC know where the storage is using `dvc remote add`. 
3. DVC will update the .dvc configuration file that will have information on the storage location. 
4. You can then add data with `dvc add <directory>`. This will create a .dvc file for this directory that tracks size, number of files, etc. 
5. Automatically adds these .dvc directories to .gitignore. This prevents large data that is tracked by DVC by being tracked by Git, since files may be too large to track. 
6. You can then add changes using `dvc push`. 
7. Check status with `dvc status` to see changes that have been made. Then proceed with `dvc add <directory>` followed by `dvc dvc push`. You can commit these changes in Git to track in GitHub! 
8. If you want a previous version, you can go to commit history in GitHub. Copy the commit link and then do `git checkout <commit link>`. This will change the DVC metadata file to previous state. 
9. Then you would do `dvc pull` to pull that previous version of files. Finally, commit your changes on Git. 

There are options for optimizing for large dataset [that are detailed here](https://dvc.org/doc/user-guide/data-management/large-dataset-optimization#large-dataset-optimization).   

Currently, DVC tracks individual files. But if you have a .zip file and change a file within that, DVC still treats the .zip as one file. So track at the individual level. 

Dealing with multiple datasets and locations managed by    using a Data Registry that uses different types of data storage. 

# Pipelines  

Check out the YouTube talk, "I don't like notebooks" about Jupyter notebooks. 

With Jupyter and R projects, you sometimes will load things or run things in different orders - and the order of this isn't tracked. Model training is particularly difficult in this regard. How do you keep track of what you have run and tried before and reproduce the results? 

How do you define dependencies between codes and model steps? Which stages does a change affect?  

DVC pipelines help with this and allow you to experiment quickly. It can answer questions like:  

- What was used to produce a model? 
- Can you compare model versions? 
- Will you be able to reproduce them later? 

The goals we want to achieve: 

1. Achieve best performance 
2. Ensure reproducibility 
3. Minimal set up and dependency - we want open source! 

VSCode is an easy user interface to run DVC to accomplish these goals.  

All free and open source! 

A DVC pipeline is defined as `dvc.yaml` file. Within this file, you can specify different pipeline stages. This will list the command you are running with dependencies (script, data), parameters, and output. Also works for running Jupyter notebooks (using `papermill` tool).  

VSCode then has an extension for DVC that shows these pipelines and experiments to show you the results that come from experiments. It tracks all the code changed for each iteration. 

`dvc dag` shows you a visualization of the pipeline, the order of stages, and the dependencies and sequence of tasks.  

The VSCode extension will show plots and results of the experiment pipelines. You can essentially change one parameter at a time easily with a line of code and it will re run the entire pipeline at the relevant steps. Outputs are then versioned in DVC and can be pushed to Git. Pretty cool! 

More information on [pipelines is on their website](https://dvc.org/doc/user-guide/pipelines).  

# Backing up data 

*Allows for additional back up of data!* 

DVC can allow you to have multiple remote storage locations - Storage 1, for example, is the default and then you can have as many additional remotes as you want. For example, you could have SSH and Drive both assigned. `dvc pull` and `dvc push` goes to default (storage 1), but you could specify to push and pull to a different storage location as a method for duplication and backup. This would allow changes to be pushed and pulled from multiple locations as backup. 

*What if you are using an HPC to do the computing for you?*  

If you are using cloud/HPC/slurm computing, you can track file changes that occur after the computing is done on slurm/HPC separately.  

# Take Homes 

We can definitely use DVC to track bioinformatic files, sequences, and other large files (images, etc.) that you can't add to GitHub! This will help us version control our files. We can also use multiple remote storage options to back up data. This would help us version control our shell scripts and bioinformatic files along the way (e.g., BAM files, trimmed reads). It will also track changes to images or other large files. 

# Resources  

The best resource will be the [website](https://dvc.org/doc) or checking out the [GitHub page](https://github.com/iterative/dvc) and submitting an issue. 



