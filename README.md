This repository contains the configuration files, simulation routines and plotting scripts for
the paper `A PDE-Constrained Generalized Nash Equilibrium Approach
for Modeling Gas Markets with Transport` by Grimm et al.

## Setup

To install the required dependencies for the julia code, from the toplevel directory run
```
julia --project=./julia
] instantiate
```
Follow the instruction at https://github.com/jump-dev/Gurobi.jl to setup `Gurobi.jl`

## Running the simulations

In each directory under `cfgs`, run

```
julia --project=../../julia ../../julia/src/run.jl ./*
```

to generate the simulation results (stored as hdf5 files).
Each subdirectory contains the following files:
- `cfg.json` contains the physical characteristics of the pipe, the market data and the operational steady state.
- `solver_opts.json` contains the solver options
- `model_opts.json` contains the penalization parameter values as well as some discretization choice.
- `solver` contains the name of the solver
- `nb_x_pts` contains the number of (discretization) nodes along the pipe.

## Generating the figures

The required python packages include `matplotlib`, `numpy`, `colorama`.

To generate the figures displaying the various quantities of interest, run
```
make -f $(git rev-parse --show-toplevel)/python/figures/Makefile
```
from the directory containing the hdf5 files of interest.

To compare the evolution of the economical quantities with different number of players, run
```
python $(git rev-parse --show-toplevel)/python/figures/fig-econ-vals.py --dirs cfg-price-3-firms,cfg-price-perf-comp
```
from the `cfgs` directory.


