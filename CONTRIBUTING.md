# Contributing to ShapePipe

ShapePipe is an open source galaxy shape measurement pipeline developed within the [CosmoStat](http://www.cosmostat.org/) lab at CEA Paris-Saclay.

These guidelines are intended to help you contribute to improving this package by fixing bugs, adding new features or simply making suggestions. Further details about ShapePipe can be found in the [documentation](https://cosmostat.github.io/shapepipe).

## Contents

1. [Introduction](#introduction)  
2. [Issues](#issues)   
3. [Pull Requests](#pull-requests)  
   1. [Before opening a PR](#before-opening-a-pr)  
   2. [Opening a PR](#opening-a-pr)  
   3. [Revising a PR](#revising-a-pr)
   3. [After a PR is merged](#after-a-pr-is-merged)  
   4. [Content](#content)  
   5. [CI tests](#ci-tests)   
   6. [Style guide](#style-guide)  

## Introduction

If you wish to contribute to the development of ShapePipe we kindly requested that you to adhere to the following guidelines and the [code of conduct](./CODE_OF_CONDUCT.md).

## Issues

The easiest way to contribute to ShapePipe is by raising a [new issue](https://github.com/CosmoStat/shapepipe/issues/new/choose). This will give you the opportunity to ask questions, report bugs or even request new features.

Remember to use clear and descriptive titles for issues. This will help other users that encounter similar problems find quick solutions. We also ask that you read the available [documentation](https://cosmostat.github.io/shapepipe) and browse [existing issues](https://github.com/CosmoStat/shapepipe/issues) on similar topics before raising a new issue in order to avoid repetition.  

We provide templates for various different types of issues that should help you clearly communicate the problems and/or ideas you have. If you are uncertain which template to choose we suggest you opt for [Help!](https://github.com/CosmoStat/shapepipe/issues/new?assignees=&labels=help+wanted&template=help-.md&title=%5BHELP%5D) and we will do out best to help you out.

## Pull Requests

If you would like to take a more active roll in the development of ShapePipe you can do so by submitting a [pull request](https://github.com/CosmoStat/shapepipe/pulls). A Pull Request (PR) is a way by which a user can submit modifications or additions to the ShapePipe package directly. PRs need to be reviewed by the package moderators and if accepted are merged into the [develop](https://github.com/CosmoStat/shapepipe/tree/develop) branch of the repository.

Before opening a PR, be sure to carefully read the following guidelines.

### Before opening a PR

The following steps should be followed before opening a pull request.

1. Log into your [GitHub](https://github.com/) account or create an account if you do not already have one.

2. Go to the main [ShapePipe repository](https://github.com/CosmoStat/shapepipe) page.

3. Check if an [issue](#issues) already exists for the changes you plan to make. If not, [open one](https://github.com/CosmoStat/shapepipe/issues/new/choose). In either case, specify in the issue that you plan to open a PR resolve the points raised.

3. Fork the repository, *i.e.* press the button on the top right with this symbol <img src="https://upload.wikimedia.org/wikipedia/commons/d/dd/Octicons-repo-forked.svg" height="20">. This will create an independent but linked copy of the repository on your account.

4. Clone your fork of ShapePipe.  

```bash
git clone git@github.com:<YOUR_GITHUB_ACCOUNT>/shapepipe.git
```

5. Add the original repository (*upstream*) to remote.

```bash
git remote add upstream git@github.com:CosmoStat/shapepipe
```

### Opening a PR

The following steps should be followed to make a pull request:

1. Pull the latest updates to the original repository.

```bash
git pull upstream develop
```

2. Create a new feature branch for your modifications.

```bash
git checkout -b <BRANCH NAME>
```

3. Make the desired modifications/additions to the relevant modules adhering to the [style guide](#style-guide).

4. Add the modified files to the staging area.

```bash
git add .
```

5. Make sure all of the appropriate files have been staged. Note that all files listed in green will be included in the following commit.

```bash
git status
```

6. Commit the changes with an appropriate description.

```bash
git commit -m "Description of commit"
```

7. Push the commits to a branch on your fork of ShapePipe.

```bash
git push origin <BRANCH NAME>
```

8. Make a [pull request](https://github.com/CosmoStat/shapepipe/compare) for your branch with a clear description of what has been done, why and what issues this relates to.

### Revising a PR

Once you have opened a PR, one of the ShapePipe maintainers will assign a reviewer to check the changes you have made and if they appropriately address the corresponding issue.

The reviewer may ask you questions about your implementation and may request changes. This process should be seen as a discussion with the objective of arriving a the best possible solution to the corresponding issue and ensuring that the standards of the code are retained.

To implement these changes simply repeat steps 3-7 from the [previous section](#opening-a-pr).

Please be patient as this process may take some time and may depend on the availability of reviewers. You can track the progress of this processes as the reviewer ticks off each bullet point in the *Reviewer Checklist*.

Once your PR has been approved the reviewer will merge this into the [develop](https://github.com/CosmoStat/shapepipe/tree/develop) branch.

### After a PR is merged

If your PR is accepted and merged it is recommended that you follow these steps to keep your fork up to date.

1. Make sure you switch back to your local `develop` branch.

```bash
git checkout develop
```

2. Delete the local branch you used for the PR.

```bash
git branch -d <BRANCH NAME>
```

3. Pull the latest updates to the original repository, which include your PR changes.

```bash
git pull upstream develop
```

4. Push the commits to your fork.

```bash
git push origin develop
```

### Content

Every PR should correspond to a bug fix or new feature issue that has already been raised. When you make a PR be sure to tag the issue that it resolves (*e.g.* this PR closes issue #1). This way the issue can be closed once the PR has been merged.

The content of a given PR should be as concise as possible. To that end, aim to restrict modifications to those needed to resolve a single issue. Additional bug fixes or features should be made as separate PRs.

### CI tests

Continuous Integration (CI) tests are implemented via [GitHub Actions](https://github.com/features/actions). All PRs must pass the CI tests before being merged. Your PR may not be fully reviewed until all CI test are passed. Therefore, try to resolve any issues in your PR that may cause the tests to fail.

In some cases it may be necessary to modify the unit tests, but this should be clearly justified in the PR description.

### Style guide

All contributions should adhere to the following style guides currently implemented in ShapePipe.

1. All code should be compatible with the Python package versions listed in the [Conda environment](https://github.com/CosmoStat/shapepipe/blob/develop/environment.yml).

1. All code should adhere to [PEP8](https://www.python.org/dev/peps/pep-0008/) standards.

1. Docstrings need to be provided for all new modules, methods and classes. These should adhere to [numpydoc](https://numpydoc.readthedocs.io/en/latest/format.html) standards.

1. When in doubt look at the existing code for inspiration.

Some particular style conventions to note are:

- We prefer single quotes `''` to double quotes `""` for strings.

```python
good_string = 'my good string'
bad_string = "my bad string"
```

- We prefer explicit floats to implicit floats.

```python
my_good_float = 1.0
my_bad_float = 1.
```

- Long lines (i.e. >79 characters) should be split using `()` as follows.

```python
my_short_line = ...
my_long_line = (
   ...
)
def my_short_func(arg1, arg2):
  ...
def my_long_func(
   first_arg=arg1,
   second_arg=arg2,
   third_arg=arg3,
   ...
):
  ...
```

- Strings should be formatting using [f-strings](https://realpython.com/python-f-strings/) and concatenated explicitly with a `+` at beginning of each line.

```python
my_short_str = f'This is a short string with param value {param_val}.'
my_long_str = (
   f'This is my long string with several parameters values {param_val1}, '
   + f'{param_val2}, {param_val3}, etc.'
)
```
