# spatial-hypergraph-epidemics

This repo is the code and data accompanying the papers [A spatial hypergraph model where epidemic spread demonstrates clear higher-order effects](https://arxiv.org/abs/2410.12688) for ComplexNetworks2024 and it's follow up paper which is currently under review at PLOS Complex Systems. 

The main code for running epidemic experiments is `code/hysteresis/epidemic-experiments.jl` and `code/hysteresis/quick-epi-experiments.jl`.

The old epidemic code (for the conference paper experiments) is located at code/`parallel-epidemics.jl`.

## Blurb on Aggregated Data
The bulk of the data is in the keys "ntransmissions_hyperedge" and "ninfected_hyperedge". These are vectors of dictionaries. The maximum possible dimensions in a file are 
tmax X max_hsize X total_simulations but they are stored as dictionaries as we expect quite a bit of sparsity in those quantities.

Infection data is stored in 2 locations.

## TODO
- clean up epidemic plotting code [X]
  - map figures in paper to code and split into separate files [X]
    - figure10 - total infections demo [X]
    - figure11 - no bistability [X]
    - figure12 - total infections [X]
    - figure13 - hyperedge effects p1 [X]
    - figure14 - hyperedge effects p2 [X]
    - figure17 - transmissions [X]
  - centralize data loading for epidemic figures [X]
    - pushed into epidemic utils
  - clustering coefficient figure [X]
    - split data generation and figure
    - cache data
    - retest plot
  - hypergraph stats [X]
    - split data generation and figure
    - cache data
    - retest plot
  - alpha figure [X]
    - split data generation and figure
    - cache data
    - retest plot
  - ppr figure [X]
    - split data generation and figure
    - cache data
    - retest plot
- unify code for _custom_heatmap (skipping - not super important)
  - hypergraph_stats
  - epidemic figures
  - retest plotting code 
- get MWE for experiments 
  - figure out parameters for experiments 
  - wrap into a function
  - have 4 function calls or wrap into a single one. 
- rename files("..hysteresis.."->"..epidemics...") [X]
- build out aggregated data [X]
  - keep all infs but aggregate transmissions [X]
  - retest plotting code [X]
- compress raw data and send it to spectral for long term storage [waiting on data caching]
- compress aggregated data
  - if small enough, put on github otherwise link to zenodo
- end-to-end example for reproducing plots from github/info
  - git clone into new directory
  - try to regenerate plots
  - code to setup environment?

