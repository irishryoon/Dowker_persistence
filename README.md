# Dowker persistence

This repository contains code for computing the [Dowker persistence diagrams](https://arxiv.org/abs/1608.05432) of two point clouds in dimensions 0 or 1. Please note that we use the term "Dowker" and "Witness" interchangeably.  

## Install
Download [Julia](https://julialang.org/downloads/).

## Quick start

```
import Pkg
Pkg.activate(".")
include("dowker_persistence.jl")
using .Dowker
using Eirene
using DelimitedFiles

# load points
P1 = readdlm("toy_example/R1_coordinates.csv", ',', Float64)
P2 = readdlm("toy_example/R2_coordinates.csv", ',', Float64)

# compute cross-system distance
D1, D2, D_P1_P2, D_P2_P1 = compute_distance(P1, P2)

# compute Dowker persistence
W_P1_P2 = compute_Witness_persistence(D_P1_P2, maxdim = 1)
W_barcode = Eirene.barcode(W_P1_P2["eirene_output"], dim = 1)

# plot Dowker persistence diagram
plot_PD(W_barcode, title = "Dowker persistence diagram")
```

## Examples
* See `1_compute_dowker_persistence.ipynb` for example computation of Dowker persistence diagram, cycle representatives, and persistence images.
* See `2_more_examples.ipynb` for Dowker persistence diagrams on a few toy examples. 