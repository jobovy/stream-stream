# stream-stream

This small repository contains the code to do the analysis in [Bovy
(2015)](http://arxiv.org/abs/1512.XXXXX) of the interaction between a
stellar stream and a disrupting dark-matter halo.

## AUTHOR

Jo Bovy - bovy - at - astro - dot - utoronto -  dot - ca

## Requirements

* [galpy](https://github.com/jobovy/galpy)
* [NEMO](bima.astro.umd.edu/nemo/)

and the usual scientific Python packages (Numpy, Scipy, matplotlib,
seaborn).

## Code overview

### (StreamKicks.ipynb)[py/StreamKicks.ipynb]

Brief notebook that computes the kicks due to the interaction of a
stellar stream with a dark-matter stream in the impulse approximation.

### (SnapshotAnalysis.ipynb)[py/SnapshotAnalysis.ipynb]

Notebook analyzing *N*-body simulations of the interaction between a
stellar stream and various dark-matter streams. The initial conditions
for the *N*-body simulation are computed in (this
notebook)[py/Orbits-for-Nbody.ipynb].

The *N*-body simulations are run using gyrfalcON and
[NEMO](bima.astro.umd.edu/nemo/) and use a variety of NEMO tools. The
initial conditions and the necessary commands to run the *N*-body
simulations are given in (the sim directory)[sim/].

