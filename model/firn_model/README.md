# Markov Chain Monte Carlo Example for Glaciology reading seminar

This is a worked example inverse problem using a simple Monte Carlo Method.

## Dependencies

The example is written in Python (tested with 3.7) but should be fine with >=2.7. You can easily get python from [anaconda](https://www.anaconda.com/download/). You also need the scientific python stack (numpy, scipy, and matplotlib) and jupyter. If you use the bare miniconda installation, you will then need to `conda install numpy scipy matplotlib jupyter`. If you are running a different python than anaconda, you will want to do `pip install numpy scipy matplotlib jupyter`.

## Getting and running the example

You need [git](https://git-scm.com/) to work with github. Git is the version control software, github is in some sense clever storage that lets you share with others and can be a backup. You don't need to use github to get the benefits of version control (get the old non-broken version of your code back)!
You don't need to learn much git--you want to know clone, commit, pull, and push. You can probably figure out the rest later as needs arise. For now, just get git and from a folder that you want to contain the project clone this repository with `git clone https://github.com/dlilien/mcmc_seminar.git`. Go into the directory (`cd mcmc_seminar`) and then (assuming you let conda modify your path) `jupyter-notebook`.

![alt text](https://imgs.xkcd.com/comics/git.png "Git according to xkcd").
