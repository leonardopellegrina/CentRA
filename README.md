# CentRA: Efficient Centrality Maximization with Rademacher Averages

This folder contains the source code of CentRA, a randomized algorithm for
approximating the centrality of sets of nodes in a network based on progressive sampling, presented at KDD 2023: https://doi.org/10.1145/3580305.3599325


An extended version of the KDD 2023 paper is available in arXiv: https://arxiv.org/abs/2306.03651


For any questions or bugs please contact me at [leonardo.pellegrina@unipd.it](mailto:leonardo.pellegrina@unipd.it)


## Installation

The software requires the [OpenMP API](http://openmp.org/wp/).
Build the software with the `make` command, obtaining the `centra` executable.

## Running CentRA

The Python script `run_centra.py` allows to run CentRA by taking several arguments:

```
usage: run_centra.py [-h] [-db DB] [-e EPSILON] [-d DELTA] [-k K]
                     [-c MCTRIALS] [-p PROGRESSIVE] [-m NUMSAMPLES] [-t TYPE]
                     [-u UNIONBOUND] [-o OUTPUTPATH] [-r RESULTSPATH]
optional arguments:
  -h, --help            show this help message and exit
  -db DB                path to graph input file
  -e EPSILON, --epsilon EPSILON
                        approximation accuracy parameter (in (0,1))
  -d DELTA, --delta DELTA
                        approximation confidence parameter (in (0,1), def. 0.1)
  -k K, --k K           number of nodes to select in the greedy algorithm
  -c MCTRIALS, --mctrials MCTRIALS
                        number of mc trials to do to estimate the mcera
  -p PROGRESSIVE, --progressive PROGRESSIVE
                        1 means progressive sampling, 0 fixed number of samples (def. progressive)
  -m NUMSAMPLES, --numsamples NUMSAMPLES
                        number of samples to take when progressuve = 0
  -t TYPE, --type TYPE  type of graph. 1 means directed, 0 undirected (def. undirected)
  -u UNIONBOUND, --unionbound UNIONBOUND
                        1 use the union bound to check the stopping condition
  -o OUTPUTPATH, --outputpath OUTPUTPATH
                        output file path to use
  -r RESULTSPATH, --resultspath RESULTSPATH
                        output file path to use for results
```

## Reproducing the experiments

Download the graphs from the following link: https://tinyurl.com/CENTRAgraphs .
Extract the zip within the `graphs` folder.
Run the Python scripts `run_all_experiments_fixed.py`, `run_all_experiments_prog.py`, and `run_all_experiments_trials.py` to reproduce all experiments described in the paper.
All results are found in `.csv` files produced by the scripts.

## Licensing

This work is licensed under the Apache License 2.0.
