# How to Fast Forward in Git

If you have been following class using Git, you likely have been making notes in the files in your cloned copy of the repo. Git has a Fast Forward function that will help you update the files in the repo each week, while not overwriting your notes. This how-to describes how, however many functions must be preformed using the Git command line. You can find a command line in RStudio next to the "console" window

## 1. Clone into the Repository

Clone into the repository, either by using RStudio's cloning GUI or through the terminal. Put the repository in a good location for storing files and note the repo path on your computer.

Using git on the command line, the following command clones the repo into where the current working directory of the terminal happens to be:

`git clone git@github.com:jonathantemplin/BayesianPsychometricModeling.git`

Once you have created your cloned repo, change the current working directory of the terminal to the highest-level folder of the repo:

`cd [your path to repo]/BayesianPsychometricModeling`

## 2. Create a new branch in the repository

Next, to make things easy to maintain, create a local branch of the master branch of the repo. A local branch will allow you to make changes to the files (these will be your class notes). You can create the local branch with the following command:

`git checkout -b [name of your branch]`

Now, when you start RStudio by opening the .Rproj file, you can tell you are in your branch by looking at the Git tab and seeing the name of it in the upper right-hand corner.

## 3. Create a Notes Folder

Next, create a folder in your repo just for notes. Here, put all files you modify from lectures or homework to be sure to easily merge your branch with an updated master branch:

On Linux or macOS:

`mkdir Notes`

Here, I copy Lecture 8 to notes:

`cp -R Lectures/08-MCMC-Algorithms/ notes/Lecture08/`

Now, you can take notes in the files of Lecture 08, and then commit them to your local branch when you are done.

## 4. Updating from the Master Branch

When lecture starts again, now you will need a new set of notes, so you want the latest version of the master branch. You can merge it with your local branch so that you don't overwrite any of your notes in your notes folder. To do this, use the following commands in the terminal:

`git fetch`


`git rebase origin/master`

Note: The rebase command is very powerful. Any changes you made to files that also changed on the master branch will be overwritten. That's the reason for creating a separate notes folder.
