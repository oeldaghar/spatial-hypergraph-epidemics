# spatial-hypergraph-epidemics

This repo is the code and data accompanying the papers [A spatial hypergraph model where epidemic spread demonstrates clear higher-order effects](https://arxiv.org/abs/2410.12688) for ComplexNetworks2024 and it's follow up paper which is currently under review.

The main code for running epidemic experiments is code/hysteresis/epidemic-experiments.jl and code/hysteresis/quick-epi-experiments.jl.

The code for the conference version is located in the branch [ComplexNetworks2024](https://github.com/oeldaghar/spatial-hypergraph-epidemics/tree/ComplexNetworks2024).


## Data

The compressed data is publically available at https://zenodo.org/records/15701837. Note the the compressed file is ~1GB while the uncompressed version is ~3GB. 

The data can be uncompressed by downloading the compressed file into the repo and running 

```tar -xzvf data.tar.gz ```

Once uncompressed, the epidemic data is stored in **data/epidemics/sirs/**.
Each file is stored as a json file. Loading in all epidemic data files can take a few minutes. 

Graphs used for epidemics can be found in **data/hypergraphs** where **\*.xy** files give the spatial encoding for each node. **\*.txt** files give the list of hyperedges.