---
layout: post
title: 2023 UW Software Carpentry Workshop in Terminal Git and R
date: '2023-02-27'
categories: Analysis Tutorial
tags: R GitHub
---

# Workshop: Software Carpentry in R - Terminal, GitHub, and R 

Workshop lessons for Terminal, GitHub, and R from UW Data Science Software Carpentry workshop.  

## Day 1: Terminal  

### Links for today and set up

[Live terminal review](https://staff.washington.edu/naomila/live/terminal.html.html)   
*Download, unzip, and put on desktop*   

[Course Website](https://uwescience.github.io/2023-02-27-uw-online/)  

[Today's lesson plan](https://swcarpentry.github.io/shell-novice/)   

Today we are learning the basics terminal and Git and will get into R specifics on Wednesday and Thursday. I will include the major points in this notebook post. Most of this is review information. 

*Note: In Mac, you may be using the `zsh` shell rather than `bash`. To start bash you can enter `bash` in Terminal.* 

### Using Terminal 

We can use Terminal to do work in our computer than can be automated - including organizing or renaming files, running scripts, searching/finding things, etc. 

**Navigating your computer:** 

- Commands are case sensitive and you can't use your mouse to move around the text 
- You can call back previous commands with the up arrow and move through commands with up and down arrows 
- Enter the command with `enter`
- `ls` to list commands, `ls Desktop` to list desktop contents without changing the directory 
- What if you have spaces in your directory? You need to surround add different slashes: `ls Desktop/Folder\ With\ Spaces\"
- List all files in subfolders directory: `ls -R Desktop`
- List help for command with `man command`and `q` to exit this window. Use these help windows to help with commands. 
- Use `ls -a` to list all/hidden files 
- Use `ls -l -h` to print in long format in human readable format for memory (e.g., print in MB or GB) 
- Use tab to auto complete file/folder names. Tip: press tab twice to have it show you all the possible options for autocomplete! 
- `pwd` to print present working directory. To move directories, use `cd Folder`.
- Use `clear` to clear terminal 
- A ~ is shorthand for your home directory e.g., `cd ~` moves you to your home directory. On Mac, you can also just do `cd` to get back to home. 
- `cd ..` moves back one directory 
- `cd -` moves you to the directory you were in just before the one you are in now - can help you bounce back and forth between two directories. 
- Absolute path = starting from the root directory of your computer (e.g., /Users/ashuffmyer/Desktop/shell-lesson-data). Relative path = starting from current working directory (e.g., Desktop/shell-lesson-data). Absolute paths always start with a slash. 
- A dot refers to the current directory. `cd .` doesn't move you to anywhere new. For example, `cd north-pacific-gyre` and `cd ./north-pacific-gyre` are the same. 

**Working with files**  

- `cat file.txt` will display contents of this file with only as much as will fit the screen
- `mkdir` to make a new folder. You can create multiple directories at the same time, `mkdir project data thesis` creates three folders (project, data, and thesis). 
- To create a new file and edit text, use `nano file`
- In nano, you can move throughout the file like a text editor, but limited. 
- To exit, use control+x and then save the file. 
- You can also use `less` to scroll through file in terminal with `q` to exit. 
- To move and rename a file, you can use `mv file.txt otherfolder/newname.txt`
- To copy a file, use `cp file.txt otherfolder/newname.txt`. This will keep the original file in the original location. 
- You can do this with multiple files all in the same command separated by spaces with the last thing specifying the destination.
- To copy a folder, you have to use a few other options, `cp -R folder newlocation'. 
- Use `rm file.txt` to delete a file and `rm -R folder` to remove a folder. 

### Batching multiple actions together: pipes and loops   

- We can use piping and filtering to take output of one command right into another. 
- To get length/number of lines of file, use `wc file`. Use `man wc` for more information. You can list multiple files by spaces in this command. 
- To select all files without listing them manually, use a "wildcard". E.g., `wc -l *.pdb`. Selects all files matching a pattern. E.g., `wc -l *th*.pdb` select any file with th in the name and .pdb file extension. 
- You can use loops to repeat functions across many files. [More info in the lesson plans](https://swcarpentry.github.io/shell-novice/05-loop/index.html).   
- If you want to output information into a file, you use "output redirection" to take the output of a command and save into a text file. E.g., `wc -l *.pdb > lengths.txt` puts file lengths for .pdb files into a new .txt file. 
- Use `sort -n file.txt` to display sort file by numeric value. There are other flags you can use for dates, reverse order, etc. To save, add > new.txt to save sorted file to new file. 
- To print, use `head -n 1 file.txt` to print 1 line of file. Use `tail` to view the end of a file in the same way.  
- In more complex cases, you are creating intermediate files that need to go into a final output. This is called "piping". This doesn't save intermediate files. A pipe is "|". For example, `wc -l *.txt | sort -n | head -n 1` counts number of lines, sorts them numerically, and displays the first line. 
- Obtain unique lines with `uniq -c file.txt`. Before this you need to `sort` because it will identify unique as something different than adjacent lines. 
- `echo` prints specified output to terminal. You can `>` after echo to put the output into a file. To add something new to a file without overwriting it, you can do `echo new text >> file.txt`. Otherwise, whatever you echo into a file will overwrite the contents. 
- The `cut` command can separate parts of a file. For example, `cut -d , -f 2 text.csv` will separate components by a comma delimiter and select only the second column. This can be piped to a new file if you want. `cut -d , -f 2 animals.csv | sort | uniq` will cut by commas, keep the second column, sort by order, and keep only unique lines. The `cut -d, -f 2 text.csv | sort | uniq -c | wc -l` command will cut by commas, keep the second column of the text.csv file, sort and find unique entries, then list the number of unique entries in terminal.  
 
*Example from workshop data (download link above):*   

```
cd Desktop/shell-lesson-data/exercise-data/proteins

cat cubane.pdb 
#shows output of file 

wc cubane.pdb
#shows number of lines and words in file 

wc -l cubane.pdb
#shows only number of line

wc -l *.pdb
#shows line numbers for all files with .pdb

wc -l *.pdb > lengths.txt
#put file lengths in a new file

sort -n lengths.txt
#display sort by file length (numeric)

sort -n lengths.txt > sorted.txt
#save sorted list 

head -n 3 sorted.txt
#view top 3 lines 

tail -n 3 sorted.txt
#view bottom 3 lines

sort -n lengths.txt | head -n 1
#output sort contents into head command, doesn't save an intermediate sorted file! 

#put it all togther with pipes! 
wc -l *.pdb | sort -n | head -n 1

echo Hello > greeting.txt 
#new file that says "Hello"

echo Goodbye > greeting.txt 
#file now only contains "Goodbye"

echo Hello Again >> greeting.txt
#file now has "Hello" and "Goodbye"

```

Some example exercises are [found here using other datasets and pipes](https://swcarpentry.github.io/shell-novice/04-pipefilter/index.html#tools-designed-to-work-together). 

Loops perform actions over multiple files. An example structure of a loop looks like this: 

```
for thing in list_of_things
do
    operation_using $thing    # Indentation within the loop is not required, but aids legibility
done
```

As an example, 

```
$ for filename in basilisk.dat minotaur.dat unicorn.dat
> do
>     echo $filename
>     head -n 2 $filename | tail -n 1
> done

This returns the following output: 

basilisk.dat
CLASSIFICATION: basiliscus vulgaris
minotaur.dat
CLASSIFICATION: bos hominus
unicorn.dat
CLASSIFICATION: equus monoceros
```

In a loop, `$` designates the variables value and for the shell to treat the variable as a name and substitute a value in its place, rather than a text or command. For example, adding a `$filename` in the loop above substitutes in the file names from the file name list we created.  

Also, `$filename` is equivalent to `${filename}`.  

You can call your variable anything! Here is another example calling the variable "x" rather than "filename". 

```
$ for x in basilisk.dat minotaur.dat unicorn.dat
> do
>     head -n 2 $x | tail -n 1
> done
```

This loop below would list all files ending with .pdb and add them into a new file called all.pdb: 

```
for datafile in *.pdb
do
    cat $datafile >> all.pdb
done
```

If you want to create a back up copy for lots of files, for example, you can loop the `cp` command: 

```
$ for filename in *.dat
> do
>     cp $filename original-$filename
> done
```

This creates new files that have "original" in the name. 

As a tip, use `echo` to make sure its doing what you think it is!

You can separate commands by semicolons in a loop as well. 

If you want to add a shell script as the command for a loop, it would look like this: 

```
for datafile in NENE*A.txt NENE*B.txt; do echo $datafile;
bash goostats.sh $datafile stats-$datafile; done
```

In this above example, it is taking datafiles that start with NENE and end with either A or B, echoing the name of that file, and running a bash script. This script takes two inputs: the file input (`$datafile`) and an output file (`stats-$datafile`). 

To look at the history of commands and tasks, you can run: `history | tail -n 5`. THis hows the recent history. You can then re-run tasks by specifying `!` followed by the task number from the history list. 

Nested loops are also possible and would look something like this:  

```
for species in cubane ethane methane
do
for temperature in 25 30 37 40
do
mkdir $species-$temperature
done
done
``` 
To enter these in shell, use the above notation and copy paste. 

### Shell Scripts 

A shell script can be used to store commands in a file and then run in an easy command. 

Make a shell script by running `nano script.sh`, paste/write in your script (best practice is to prepare script in text editor outside terminal then paste in), and then run. 

For example: 

```
nano script.sh

head -n 15 "$1" | tail -n 5

bash script.sh file.txt
```

This will take the first (`$1`) input (file.txt) and return the first 15 lines and the last 5 lines of the file. 

You can specify as many arguments as you want. For example, here is a shell script that takes three inputs: 

```
# Select lines from the middle of a file.
# Usage: bash script.sh filename end_line num_lines

nano script.sh

head -n "$2" "$1" | tail -n "$3"

bash script.sh file.txt 15 5
```

This will return the first 15 lines of file.txt and the last 5 lines of this file. 

If you have lots of files to run a script over, you can use the `$@` to indicate a list of files. For example, 

```
nano do-stats.sh

#put this in the nano editor, note the script here refers to some other script to calculate stats 

for datafile in "$@"
do
    echo $datafile
    bash script.sh $datafile stats-$datafile
done

#then run this script
bash do-stats.sh *.txt

```

To streamline this a bit, you can also do: 

```
nano do-stats.sh

#put this in the nano editor, note the script here refers to some other script to calculate stats 

for datafile in *.txt
do
    echo $datafile
    bash script.sh $datafile stats-$datafile
done

#then run this script
bash do-stats.sh 
```

The trade off is that you have a simpler command to run this script, but to change the files you want to run this over, you have to edit the .sh file. If you use the `$@` operator instead, you can just specify the files when you run the script in the command line. 

### Finding things in Terminal 

The `grep` command helps you find things. It is a contraction of global/regular expression/print terms. It will find and print lines in files that match a pattern. 

Use the command like this:  

```
grep not haiku.txt
```
This will return all lines with "not" in the document haiku.txt. 

To limit the searches to word boundaries (for example, searching for "the" could return lines with "the" and "thesis") use `grep -w pattern file.txt`.  

You can also search for a phrase with `grep -w "is not" haiku.txt`.  

Adding `-n` will return line numbers with the matching phrases, which can be helpful to subset data. 

Adding `-i` will make the search case sensitive.  

Adding `-v` will reverse the search! So it will look for lines that do NOT contain the phrase of interest.  

As done with other commands, adding -R or -r will search through a directory of files.  

**Using regular expressions:**  

Regular expressions are complex and powerful! There is a full lesson on this [available here](https://v4.software-carpentry.org/regexp/index.html). 

Here is an example:   

```
grep -E "^.o" haiku.txt
```

In this example, we use `-E` to tell the computer we are searching for a pattern in quotes (regular expression). The `^` anchors the start of the line. The `.` matches a single character. Alpha numeric match to themselves.   

So here, we are looking to return any lines with an "o" in the second position in a word.   

**Using find function:**  

- `find` will search for files and directories 
- To show files in current directory, use `find .`, which will return all things in directory. 
- To find directories, we can use `find . -type d`, which finds all directories in the current location. If we use `-type f`, it will return all files. 
- You can also search for files that match by name. For example, `find . -name "*.txt"` will return all files in current location that end with `.txt`. You want to put the search pattern in quotes because it will otherwise try to expand the first file name it comes across. 
- `ls` and `find` are different because ls lists everything it can, while find searches with certain properties. 

We can combine this with other commands. For example:   

```
wc -l $(find . -name "*.txt")
```

This list the line numbers for all files that have .txt as an extension. Cool!   

It’s very common to use find and grep together. The first finds files that match a pattern; the second looks for lines inside those files that match another pattern. Here, for example, we can find txt files that contain the word “searching” by looking for the string ‘searching’ in all the .txt files in the current directory:  

```
grep "searching" $(find . -name "*.txt")
```

This will return lines containing the word "searching" in each file.  

## Day 2: Git  

### Links for today and set up  

[Lesson plan for today](https://swcarpentry.github.io/git-novice/)  

We will also use GitHub today, which I already have installed and set up on my computer and online account.   

### What is GitHub and Git? 

- Git is a program for automated version control and GitHub is a platform/server that hosts Git to version control repositories. 
- Provides permanent backup and storage location for files that are openly accessible
- GitHub is the most commonly used Git server  

### Configuring and Setting Up Git  

*Note: I'm not doing these steps because its already set up on my computer*  

Set up Git configurations on your computer: 

- First, configure Git with `git config --global user.name "My Name"`, followed by `git config --global user.email "email@email.com"`. 
- The double dashes `--global` indicate to Git that you are making these settings for every time you use Git.  
- Then, use `git config --global core.autocrlf input`. This setting has to do with line returns, which is different between mac/windows.  
- Use `git config --global core.editor "nano -w"` to set nano as the text editor. 
- Finally, `git config --global init.defaultBranch main` to set branches default name as "main". 
- You only need to run these commands once on your computer. 
- Use `git config --help` for more information. 

We next need to tell Git how to communicate with GitHub. Passwords aren't very secure, so we can set up SSH authentication. Set this up with the following commands: 

- Create a public and private key: 
- `ssh-keygen -t ed25519 -C "email@email.com"`
- `-t` specifies the type of key to create, in this case we are using ed25519. 
- Save the file with `enter`. 
- Leave the passphrase empty. As long as you are the only one accessing the computer, this file will only be readable by you anyways. You can use a passphrase if there is someone else that uses the computer. 
- Basically, messages are encrypted with the public key, and it can only be decrypted with your private key. 
- Then navigate to `cd ~/.ssh` and use `ls` to view files. `cat id_ed25519.pub` and copy the entire line that comes up. 
- In your GitHub account, go to settings and SSH/GPG keys. Create a new key with whatever title you want. Then paste in the line that you copied. Now the key is stored in GitHub. 
- Now you can log into GitHub when you need to push and pull. 

### Set up a Repository  

- Create a repository in GitHub.com. 
- You can now copy the SSH repository link. 
- In Terminal, navigate to the directory you want the repository to go into. 
- Type `git clone <insert link>` to clone the respository to the current directory on your computer. 
- Use `ls -a` to view all files that start with a . (hidden files). 
- You can now add files into your repository. 

### Using Git to track and commit files 

- `git status` shows the status of your repository and changes made that have not been committeed 
- `git add <file>` adds a specific file or directory, `git add -A` adds all changes. This file is now added to the queue to be committed. 
- Commits are taking the work you have done and permanently move it to the Git repository. It will track these changes. 
- `git commit -m <message or description>` will commit the added files/changes 
- `git log` will show you recent activity with author, date, and message of commits. 
- `git diff <file>` will show the changes in one file comparing what is currently in GitHub vs the local copy. `git diff` without a file will show all changes. 
- `git restore <file>` will pull the original file from the repository (we will talk about this more later) 
- You can add as many changes as you want into a single commit 
- `git log --oneline` shows you the history with each commit on a single line! Super useful. 
- Git will not track a new folder/directory until there is something in it. It tracks files, not directories. 
- Use `git diff <time> <file>` to compare two states. For example, you can enter `git diff HEAD file.txt` will compare the file to the current state of the respository. You can also use `git diff commit# file.txt` to compare this file to a previous commit state. 
- `git diff HEAD~1 file.txt` compares the file to two commits ago. 
- `HEAD` just means the current repository state. 
- `git log` gives a longer form list of previous commits 
- Commit ID's are long, you will often see shorter alias commit ID's and you can use these for comparing commits 

### What if we make a mistake? 

- If we save a file with a mistake, we can use `git checkout` to replace a local file with the most recent one from the repository. This is useful to pull whatever is current in the repository and replace the local version. 
- E.g., `git checkout file.txt`
- To go back to the state of a previous commit (2 ago in this case), you would use, `git checkout HEAD~2 file.txt`. This will replace the file with the version from 2 commits ago. 
- You will then need to commit this change, since it is different from the current state. 

### How do we ignore files? 

- Its ok to have files that Git isn't tracking if not needed 
- You can tell Git that its ok to ignore files 
- You create a new file (in nano or other editor) called `.gitignore` 
- Each line of this file is a file you want it to ignore. Enter files or patterns (e.g., *.csv) to have it ignore files. You MUST do this before you commit and track a file. If you want to ignore files that have already been tracked or added, you have to delete those files then add them to ignore before adding them back in (e.g., DS_Store objects are notorious for this). 
- Commit and push this file. 

### Sending files to/from GitHub

- Once you have committed changes, use `git push` to send things up to GitHub. 
- This will now be reflected on your GitHub webpage. 
- Use `git pull` to pull any changes that have been sent to GitHub either from webpage or from a collaborator. 
- Commit and push at the end of major sessions or you accomplish a task/step that you want recorded. 

### Collaborating on GitHub - Forking and Pull Requests

- You can add collaborators to allow them to push and pull to your repository. Then will accept an invitation and clone the respositorty to their computer. 
- Another way is for the other user to go to your GitHub repository on GitHub.com. They can see files since its public, but they wont have permission to edit and commit. 
- In order to contribute, the user can "fork" the repository. 
- This will copy the respository to your own account. 
- Then use `git clone <linked to your forked repo>` to clone this repository fork. 
- The user can then edit contents of the respository in their forked repo. 
- You can see the number of forks on the original repository. 
- When changes are pushed, they are pushed to your forked repository. 
- You can then have your collaborator with the original repository pull changes from your fork. 
- You can do this with "Pull Requests" and "new pull request". 
- At the top there is a bar with the arrow showing the direction of the changes to be made. For example, you would want the change from the forked repo to be pulled into the original repo. 
- At the bottom, it will show you the difference between the repositories. 
- You can then add a more detailed comment and create the request. 
- Now, in the original repo, the original user will then go to Pull Requests and you will see the active requests. 
- The user can then confirm the request, which will pull those changes into the repository. 
- But in some cases if there are multiple edits to the same file, it will create a conflict. 
- You can use the "resolve conflicts" button to resolve them. You can then delete anything you don't want, or edit the file and then commit the merge. 
- It is called a pull request because the owner is pulling from your repo into theirs. 
- Users can then pull these changes down into their fork by using "Sync Fork" button. You will then update the branch and it will sync all changes.

## Day 3 & 4: R 

### Basics of using R 

Lesson plans for today are [here](http://swcarpentry.github.io/r-novice-inflammation/) and [here](https://preview.carpentries.org/R-ecology-lesson/01-intro-to-r.html). 

- String = anything enclosed in quotes
- `print("string of things")` prints the string into the output
- To glue strings together, you can use `paste("string1", "string2")`
- Use hashtags to make notes as in Terminal 
- In `paste`, you can also set the `sep` argument to paste things together with anything other than a space. E.g., `paste("banana", "apple", "orange", sep="/")` returns `banana/apple/orange`
- A variable is anything assigned a value with a `<-`. You can also use `=`, but this can be tricky because we often use `=` for other tasks. 
- You can use the `args` function to display arguments that the function requires. For example `args(round)`. The `round` function is used to round numbers (`round(1.2345, digits=2` returns 1.23). 
- Create a vector or list of multiple things with `vector <- c(1,2,3,4)`. You can add new things to the front (`vector<-c(5, vector)`) or end (vector<-c(vector, 10)`) of a vector as well. 
- Boolean logical values are `TRUE` and `FALSE`. These can be used to filter data later on. 
- Summarize variables and view structure with `str()`. 
- `vector[2]` pulls out the second value of the vectorand `vector[c(3,2)]` pulls out the third and second values. 
- Square values = retrieve something, parentheses = do something. 
- All values are coerced to strings if there is at least one string in the vector. 
- You can also subset in a vector format using `vector_new<-vector[vector<5]`. A `==` will evaluate if two things are equal.
- `!=` means "is not equal to" and you can combine <= or >= to mean less than/greater than or equal to. Combine these with `|` to mean "or". So you can filter by `vector <=  5 | vector ==10`. And `&` will look for both conditions to be met. 
- To subset by things in a list, you would filter by `vector %in% list`. Where you create a list to be `list<-c(1,2,3)`. 
- NA's can make things tricky. For example when calculating a mean, you have to use `mean(data, na.rm=TRUE)` if you have NA's. We can use `!is.na(vector)` to remove them. 

### Working with data 

Below, I'll paste in the code from our practice session.  **You can copy and paste this code into an R script and follow along!**   

To download data from a URL and store it on your computer, use this function: 

```
download.file(url = "https://ndownloader.figshare.com/files/2292169", destfile = "data_raw/portal_data_joined.csv")
``` 

This is a helpful function! 

We are going to use `tidyverse` for the data today using `library(tidyverse)` after you `install.packages("tidyverse")`. 

To read in a file after downloading, use `read_csv(path/file.csv)`. 

Set up

```
library(tidyverse)
```

Download data

```
download.file(url = "https://ndownloader.figshare.com/files/2292169", destfile = "data_raw/portal_data_joined.csv")
```

Read data 

```
surveys<-read_csv("data_raw/portal_data_joined.csv")
```

Examine and practice subsetting data
 
```
str(surveys) #view structure of dataframe
names(surveys) #column names
summary(surveys) #summary stats
head(surveys,3) #view first 3 rows
surveys[1,] #gives first row
surveys[,2] #gives second column
surveys[1:3, c("month", "day", "year")] #how month/day/year columns for rows 1 to 3
surveys[,-1] #select all rows but not the first column
surveys$species #view all species column entries
surveys[surveys["sex"]=="F",] #select all rows where sex is female
surveys[!is.na(surveys["sex"]) & surveys["sex"]=="F",] #select rows where there isn't an NA and sex is female 
na.omit(surveys[surveys["sex"]=="F",]) #omit NAs and return all rows where sex is female
```

Summarizing data and using factors 

```
summary(surveys) #this doesn't work to summarize any non-numeric/integer variable 

#factors are categorical and R has knowledge that data points belong to these categories 
surveys$sex #as character strings
factor(surveys$sex) #as factors

surveys$sex<-factor(surveys$sex) #save as factor
levels(surveys$sex) #shows the levels/order of the factor 
nlevels(surveys$sex) #shows number of categories 
surveys$sex<-factor(surveys$sex, levels=c("M", "F")) # order levels manually
levels(surveys$sex) #M is now first 
```

Plot data 

```
plot(surveys$sex) #shows counts of observations in each level of factor 
surveys$sex<-addNA(surveys$sex) #include NA's in factor levels
plot(surveys$sex) #shows the count of NA's
levels(surveys$sex)
levels(surveys$sex)[3]<-"undetermined" #rename NA's as "undetermined" - a good way to rename a level of a factor 
levels(surveys$sex) #shows this change
plot(surveys$sex) #now the plot shows M and F and undetermined 
```

Formatting dates

We will use the lubridate package for formatting dates. 

```
library(lubridate)
```

Working with dates

```
my_date<-ymd("2015-01-01") #specifies string as date from YYYY-MM-DD format
str(my_date)

my_date<-ymd(paste("2015", "1", "1", sep="-")) #this is another way we could format a date string into an R date object
str(my_date) #same result as above

surveys$date <- ymd(paste(surveys$year, surveys$month, surveys$day, sep = "-")) #we can format a new column as dates by pasting together separate year month day columns
str(surveys)
summary(surveys$date) #we can see there are 129 NA's from the date we made above

#investigate the missing dates
missing_dates <- surveys[is.na(surveys$date), c("year", "month", "day")] #extract NAs and only show year month day columns
head(missing_dates)

#the reason these dates are NA is because September and April do not have 31 days, only 30. This is an example of QC'ing your data
```

### Manipulating data 

Cheat sheets for tidyverse can be found [here](https://raw.githubusercontent.com/rstudio/cheatsheets/main/data-transformation.pdf) and [here](https://raw.githubusercontent.com/rstudio/cheatsheets/main/data-import.pdf). 

Manipulating data  

```
select(surveys, plot_id, species_id, weight) #keeps only plot, species, and weight columns
select(surveys, -record_id, -species_id) #keeps all variables except record and species
filter(surveys, year == 1995) #this will filter the dataset by evaluating a condition
```

Using pipes to do multiple actions at the same time  

```
#here is an example of nesting functions
surveys_sml <- select(filter(surveys, weight < 5), species_id, sex, weight) #keeps species sex and weight columns from the dataframe for rows with weight less than 5

#the same way to do this with a pipe is: 
surveys_sml <- surveys %>%
  filter(weight < 5) %>%
  select(species_id, sex, weight); surveys_sml
```

Hint: A shortcut for putting in a pipe is to use command + shift + m! 

Create new columns with `mutate` 

```
# filter out na's and create a new column and display the top 5 rows
surveys %>%
  filter(!is.na(weight)) %>%
  mutate(weight_kg = weight / 1000) %>%
  head()

#you can create multiple columns in one mutate command like this 
surveys %>%
  mutate(weight_kg = weight / 1000,
         weight_lb = weight_kg * 2.2)
```

Summarise data with `summarize`

```
#summarize mean weights for each sex and display
surveys %>%
  group_by(sex) %>%
  summarize(mean_weight = mean(weight, na.rm = TRUE))

#summarize by two groups and remove na's before the calculation
surveys %>%
  filter(!is.na(weight)) %>%
  group_by(sex, species_id) %>%
  summarize(mean_weight = mean(weight))

#summarize and then print 15 lines, use this to customize the output
surveys %>%
  filter(!is.na(weight)) %>%
  group_by(sex, species_id) %>%
  summarize(mean_weight = mean(weight)) %>%
  print(n = 15)

#summarize multiple metrics and arrange them in ascending order by weight
surveys %>%
  filter(!is.na(weight)) %>%
  group_by(sex, species_id) %>%
  summarize(mean_weight = mean(weight),
            min_weight = min(weight)) %>%
  arrange(min_weight)

#now arrange in descending order 
surveys %>%
  filter(!is.na(weight)) %>%
  group_by(sex, species_id) %>%
  summarize(mean_weight = mean(weight),
            min_weight = min(weight)) %>%
  arrange(desc(mean_weight))
```

Counting data 

```
#we can use the count function to calculate the number of observations in our data 
surveys %>%
    count(sex)

#this is also equivalent (but more streamlined) to: 
surveys %>%
    group_by(sex) %>%
    summarise(count = n())

#you can also sort within the count function
surveys %>%
    count(sex, sort = TRUE)

#now calculate number of observations for a combination of two factors and arrange in descending order (alphabetical) by species
surveys %>%
  count(sex, species)%>%
  arrange(species, desc(n))

#example: what was the heaviest animal measured in each year? 
surveys %>%
    filter(!is.na(weight)) %>%
    group_by(year) %>%
    filter(weight == max(weight)) %>%
    select(year, genus, species, weight) %>%
    arrange(year)
```

Reshaping data with `pivot_longer` and `pivot_wider` (previously `gather` and `spread`). 

What are the 4 rules of a "tidy" dataframe? 
1. Each variable has its own column
2. Each observatin has its own row
3. Each value has its own cell
4. Each type of observational unit forms a table 

Pivoting tables helps with rule #4!

Pivoting from wide to long

```
#pivot_wider needs three things:
# data
# names from (values will become new columns)
# values from (values will fill the new columns)

surveys_gw <- surveys %>%
  filter(!is.na(weight)) %>%
  group_by(plot_id, genus) %>%
  summarize(mean_weight = mean(weight))
str(surveys_gw)
#this data frame has observations for each plot across multiple rows

surveys_wide <- surveys_gw %>%
  pivot_wider(names_from = genus, values_from = mean_weight) #convert from long to wide format

head(surveys_wide)
str(surveys_wide)
#now each row is a plot and each column is a species populated by mean weight

#for missing values, add a weight of 0 rather than NA
surveys_gw %>%
  pivot_wider(names_from = genus, values_from = mean_weight, values_fill = 0) %>%
  head()

```

Pivoting from wide to long format 

```
#we would do this if we want to have the values in columns as a variable (e.g., genus rather a column for each genus)

#pivot_longer needs 4 things
# the data
# the names_to column variable we wish to create from column names.
# the values_to column variable we wish to create and fill with values.
# cols are the name of the columns we use to make this pivot (or to drop).

surveys_long <- surveys_wide %>%
  pivot_longer(names_to = "genus", values_to = "mean_weight", cols = -plot_id) #column names now go into a new column called genus, the data goes into a new column called "mean weight", and the function ignores the plot column

str(surveys_long)
head(surveys_long) #now we are back to the wide version

#another example: 
surveys_wide_genera <- surveys %>%
  group_by(plot_id, year) %>%
  summarize(n_genera = n_distinct(genus)) %>%
  pivot_wider(names_from = year, values_from = n_genera) #summarize the number of unique genera for each plot and year, then make into a wide format
head(surveys_wide_genera)

#now pivot so that each row is a unique plot and year combo in long format
surveys_wide_genera %>%
  pivot_longer(names_to = "year", values_to = "n_genera", cols = -plot_id)

#here is an example of how you would pivot with multiple data columns and put the measurement type as a new column
surveys_long2 <- surveys %>%
  pivot_longer(names_to = "measurement", values_to = "value", cols = c(hindfoot_length, weight))
head(surveys_long2)
```

Exporting data 

```
#obtain dataset with no NAs
surveys_complete <- surveys %>%
  filter(!is.na(weight),           # remove missing weight
         !is.na(hindfoot_length),  # remove missing hindfoot_length
         !is.na(sex))                # remove missing sex

# Extract the most common species_id
species_counts <- surveys_complete %>%
    count(species_id) %>%
    filter(n >= 50)

# Only keep the most common species
surveys_complete <- surveys_complete %>%
  filter(species_id %in% species_counts$species_id)

#write to csv 
write_csv(surveys_complete, file = "data/surveys_complete.csv")
```

### Visualizing data  

ggplot2 is the main package we will use for visualizations. It uses the following structure:  

`ggplot(data = <DATA>, mapping = aes(<MAPPINGS>)) +  <GEOM_FUNCTION>()`

Build a plot with points and try a hex plot 

```
#full plot
ggplot(data = surveys_complete, aes(x = weight, y = hindfoot_length)) +
  geom_point()

# the + sign add layers to the plots and anything you add to the main plot header (data and aes) applies to all other functions in the plot unless otherwise specified  

#you can also assign to a variable and then draw additional objects on that variable 
# Assign plot to a variable
surveys_plot <- ggplot(data = surveys_complete,
                       mapping = aes(x = weight, y = hindfoot_length))

# Draw the plot with points
#geom_point, geom_line, geom_boxplot are common geoms to add 
surveys_plot +
    geom_point()

# Draw with hex 
install.packages("hexbin")
library(hexbin)

#draw a hex plot
surveys_plot +
 geom_hex()
#the advantage of a hex plot is that you can see count frequency in addition to distribution, cool!

#add point transparency
ggplot(data = surveys_complete, aes(x = weight, y = hindfoot_length)) +
    geom_point(alpha = 0.1)

#add color
ggplot(data = surveys_complete, mapping = aes(x = weight, y = hindfoot_length)) +
    geom_point(alpha = 0.1, color = "blue")

#vary color by species
ggplot(data = surveys_complete, mapping = aes(x = weight, y = hindfoot_length)) +
    geom_point(alpha = 0.1, aes(color = species_id))
```

Build a boxplot

```
#basic boxplot
ggplot(data = surveys_complete, mapping = aes(x = species_id, y = weight)) +
    geom_boxplot()

#add data points over the box plot
ggplot(data = surveys_complete, mapping = aes(x = species_id, y = weight)) +
    geom_boxplot(alpha = 0) +
    geom_jitter(alpha = 0.3, color = "tomato")

#try a violin plot
ggplot(data = surveys_complete, mapping = aes(x = species_id, y = weight)) +
    geom_violin()

#scale by log 10
ggplot(data = surveys_complete, mapping = aes(x = species_id, y = weight)) +
    geom_violin()+
    scale_y_log10()

```

Plotting time series data 

```
#calculate the number of counts per year per genus
yearly_counts <- surveys_complete %>%
  count(year, genus)

#plot this data
ggplot(data = yearly_counts, aes(x = year, y = n)) +
     geom_line()

#tell the plot to group by genus
ggplot(data = yearly_counts, aes(x = year, y = n, group = genus)) +
    geom_line()

#now add color
ggplot(data = yearly_counts, aes(x = year, y = n, color = genus)) +
    geom_line()
```

Piping in from other functions

```
#for simplicity, you can pipe in from other tidyverse functions
yearly_counts %>%
    ggplot(mapping = aes(x = year, y = n, color = genus)) +
    geom_line()

#link data manipulation functions
yearly_counts_graph <- surveys_complete %>%
    count(year, genus) %>%
    ggplot(mapping = aes(x = year, y = n, color = genus)) +
    geom_line();yearly_counts_graph

```

Faceting

```
#facet to show multiple plots side by side split by categories of choosing
ggplot(data = yearly_counts, aes(x = year, y = n)) +
    geom_line() +
    facet_wrap(facets = vars(genus))

#now split line in each plot by sex
yearly_sex_counts <- surveys_complete %>%
                      count(year, genus, sex)

ggplot(data = yearly_sex_counts, mapping = aes(x = year, y = n, color = sex)) +
  geom_line() +
  facet_wrap(facets =  vars(genus))

#set facets by rows and columns manually 
ggplot(data = yearly_sex_counts,
       mapping = aes(x = year, y = n, color = sex)) +
  geom_line() +
  facet_grid(rows = vars(sex), cols =  vars(genus))

# One column, facet by rows
ggplot(data = yearly_sex_counts,
       mapping = aes(x = year, y = n, color = sex)) +
  geom_line() +
  facet_grid(rows = vars(genus))

# One row, facet by column
ggplot(data = yearly_sex_counts,
       mapping = aes(x = year, y = n, color = sex)) +
  geom_line() +
  facet_grid(cols = vars(genus))
```

Plot themes

```
#black and white theme 
 ggplot(data = yearly_sex_counts,
        mapping = aes(x = year, y = n, color = sex)) +
     geom_line() +
     facet_wrap(vars(genus)) +
     theme_bw()

#minimal
 ggplot(data = yearly_sex_counts,
        mapping = aes(x = year, y = n, color = sex)) +
     geom_line() +
     facet_wrap(vars(genus)) +
     theme_minimal()
 
 #light
 ggplot(data = yearly_sex_counts,
        mapping = aes(x = year, y = n, color = sex)) +
     geom_line() +
     facet_wrap(vars(genus)) +
     theme_light()
 
 #void
 ggplot(data = yearly_sex_counts,
        mapping = aes(x = year, y = n, color = sex)) +
     geom_line() +
     facet_wrap(vars(genus)) +
     theme_void()
 
#classic
 ggplot(data = yearly_sex_counts,
        mapping = aes(x = year, y = n, color = sex)) +
     geom_line() +
     facet_wrap(vars(genus)) +
     theme_classic()
```

Customizing plots

```
#change axis titles
ggplot(data = yearly_sex_counts, aes(x = year, y = n, color = sex)) +
    geom_line() +
    facet_wrap(vars(genus)) +
    labs(title = "Observed genera through time",
         x = "Year of observation",
         y = "Number of individuals") +
    theme_bw()

#change text size
ggplot(data = yearly_sex_counts, mapping = aes(x = year, y = n, color = sex)) +
    geom_line() +
    facet_wrap(vars(genus)) +
    labs(title = "Observed genera through time",
        x = "Year of observation",
        y = "Number of individuals") +
    theme_bw() +
    theme(text=element_text(size = 16))

#change text orientation and face
ggplot(data = yearly_sex_counts, mapping = aes(x = year, y = n, color = sex)) +
    geom_line() +
    facet_wrap(vars(genus)) +
    labs(title = "Observed genera through time",
        x = "Year of observation",
        y = "Number of individuals") +
    theme_bw() +
    theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
                        axis.text.y = element_text(colour = "grey20", size = 12),
                        strip.text = element_text(face = "italic"),
                        text = element_text(size = 16))

#create your own custom theme to reuse! 
grey_theme <- theme(axis.text.x = element_text(colour="grey20", size = 12,
                                               angle = 90, hjust = 0.5,
                                               vjust = 0.5),
                    axis.text.y = element_text(colour = "grey20", size = 12),
                    text=element_text(size = 16))

ggplot(surveys_complete, aes(x = species_id, y = hindfoot_length)) +
    geom_boxplot() +
    grey_theme
```

Arranging plots with `patchwork` package. In patchwork, you use `+` to put plots next to each other, `/` to arrange them vertically, and `plot_layout()` to determine the total space. More examples are here: https://patchwork.data-imaginist.com/. In patchwork you can add annotations and customize layouts or add labels. 

```
install.packages("patchwork")
library(patchwork)

plot_weight <- ggplot(data = surveys_complete, aes(x = species_id, y = weight)) +
  geom_boxplot() +
  labs(x = "Species", y = expression(log[10](Weight))) +
  scale_y_log10()

plot_count <- ggplot(data = yearly_counts, aes(x = year, y = n, color = genus)) +
  geom_line() +
  labs(x = "Year", y = "Abundance")

plot_weight / plot_count + plot_layout(heights = c(3, 2))

#examples of other layouts
#(p1 | p2 | p3) /
#      p4
#this would create a plot with three plots on top and one plot on the bottom - super cool and easier than cowplot! 
```

Note: check out patchwork for arranging plots, it seems really useful! 

Exporting plots

```
my_plot <- ggplot(data = yearly_sex_counts,
                  aes(x = year, y = n, color = sex)) +
    geom_line() +
    facet_wrap(vars(genus)) +
    labs(title = "Observed genera through time",
        x = "Year of observation",
        y = "Number of individuals") +
    theme_bw() +
    theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90,
                                     hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(colour = "grey20", size = 12),
          text = element_text(size = 16))

ggsave("fig/plot.png", my_plot, width = 15, height = 10) #save to folder

## This also works for plots combined with patchwork
plot_combined <- plot_weight / plot_count + plot_layout(heights = c(3, 2))
ggsave("fig/plot_combined.png", plot_combined, width = 10, dpi = 300)
```