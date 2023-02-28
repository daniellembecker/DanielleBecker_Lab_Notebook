---
layout: post
title: 2023 UW Software Carpentry Workshop in Terminal Git and R
date: '2023-02-27'
categories: Analysis Tutorial
tags: R GitHub
---

# Workshop: Software Carpentry in R - Terminal, GitHub, and R 

## Day 1: Terminal and Git 

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



