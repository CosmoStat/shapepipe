[Home](./shapepipe.md)

# Contribution Guidelines

Please read the following guidelines if you would like to contribute to the development of ShapePipe.

## Contents

1. [ShapePipe Development Plan](#ShapePipe-Development-Plan)
   1. [Issues](#Issues)
   1. [Milestones](#Milestones)
   1. [Branching](#Branching)
   1. [Merging](#Merging)
   1. [Documenting](#Documenting)
1. [Git Help](#Git-Help)
   1. [Cloning](#Cloning)
   1. [Remote vs Local](#Remote-vs-Local)
   1. [Pulling Updates](#Pulling-Updates)
   1. [Branching](#Branching)
   1. [Staging](#Staging)
   1. [Committing](#Committing)
   1. [Pushing to the Remote Repository](#Pushing-to-the-Remote-Repository)
   1. [Merging](#Merging)

## ShapePipe Development Plan

This document outlines the procedures everyone should follow for the continuous development of the pipeline.

### Issues

For every development task in the pipeline an issue, with a clear description, should be created and assigned to an individual and a milestone.

When an issue is created it is by default labelled as `To Do`. When development begins on this issue it should be updated to `In Progress`. Finalised issues will be labelled as `Closed`.

### Milestones

Milestones will set specific deadlines for a set of issues to be accomplished.

### Branching

New branches should be created for specific issues with the naming convention `package_developer` to specify which package is being modified and who is working on this branch.

After merging branches should be deleted.

### Merging

Merge requests should be made regularly and contain a small number of stable modifications. This reduces potential conflicts between branches.  

### Documenting

Docstrings should adhere to [Sphinx](http://www.sphinx-doc.org/en/stable/) [Numpydoc](http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html) conventions.

**[Back to Wiki](home)**

## Git Help

This page contains some tips and tricks for interacting with a Git repository.

### Cloning

To clone a repository (*i.e.* make a local copy of the remote repository) use the `git clone` command. For example to clone the pipeline:

```bash
   $ git@github.com:CosmoStat/shapepipe.git
```

This will copy the contents of the `master` branch to your computer.

### Remote vs Local

After cloning the repository you will be able to make modifications on your local copy and track updates on the "remote" (*i.e.* the online githab repository). You can display the current remote settings using the `git remote` command. For the pipeline you should see the following output:

```bash
   $ git remote -v
   origin	git@github.com:CosmoStat/shapepipe.git (fetch)
   origin	git@github.com:CosmoStat/shapepipe.git (push)
```

`origin` is the remote repository name and has the same address for both `fetch` and `push` commands.

### Pulling Updates

Before making any modifications to the pipeline you should make sure your local copy is up to date with the remote. You can update your copy with the `git pull` command.

```bash
   $ git pull
```

### Branching

Once you have an up-to-date copy of the repository you should create a branch, which is linked to an open issue (see [Development Plan](Development-plan)). You should avoid making direct commits to `master`.

You can view existing branches with the `git branch` command. For example to list all branches on the local and remote repositories:

```bash
   $ git branch -a
```

To create a new branch simply specify a branch name. In general, you should choose a name linked to an issue and specify who is working on this branch, *e.g.*

```bash
   $ git branch updateTemplate_sam
```

To switch to the new branch use the `git checkout` command:

```bash
   $ git checkout updateTemplate_sam
```

These two commands can be combined into a shortcut:

```bash
   $ git checkout -b updateTemplate_sam
```

### Staging

After modifying a file on your branch on your local repository, the file needs to be staged for a commit. Files can be added to the staging area with the `git add` command. You can add individual files:

```bash
   $ git add FILE
```

or all files within a given directory:

```bash
   $ git add .
```

Files can be removed from the staging area with the command `git reset`, *e.g.*

```bash
   $ git reset FILE
```

Finally, you can view the status of the staging area with the command `git status`. This command will list all files currently in the staging area. In addition, it will show untracked files (*i.e.* files that have been modified but not staged).

### Committing

When you are ready to submit changes to files in the staging area to your local repository you can use the `git commit` command as follows:

```bash
   $ git commit -m "Short description of the changes made"
```

Make sure to provide a clear and concise description of the changes your are committing to the repository. You should make regular commits each of which constitutes a small number of changes as the repository can always be reset to a previous commit.

### Pushing to the Remote Repository

When you have committed the changes you want you can update the remote repository with the changes made to your local repository with the `git push` command as follows:

```bash
   $ git push origin BRANCH_NAME
```

where `origin` is the name of the remote repository and `BRANCH_NAME` is the name of your branch on your local repository. If this is the first push, a new branch will be uploaded to the remote, otherwise your commits will simply be synced with the remote branch.

### Merging

Once all of the changes related to a given issue have been pushed to the remote repository a merge request can be made. To do this simply select your branch on gitlab and press the merge request button. You will need to submit a short description of the changes that will be made to `master` and assign an individual and milestone to the request. After this, the merge request will be reviewed by one of the repository administrators. If accepted, your changed will be merged into `master` and your remote branch will be deleted. You should also delete your local branch as follows:

```bash
   $ git branch -d BRANCH_NAME
```
