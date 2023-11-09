# GitHub Repository: Towards Scalable Digital Modeling of Networks of Biorealistic Silicon Neurons

Welcome to the GitHub repository for the paper titled "Towards Scalable Digital Modeling of Networks of Biorealistic Silicon Neurons." This repository contains all the code used to produce the results presented in the paper, including implementations in C (floating and fixed point), MATLAB, and Python. The neuron model used in this platform is derived from the FH neuron model as described in the paper "E. Farquhar and P. Hasler, A bio-physically inspired silicon neuron, March 2005, doi: 10.1109/TCSI.2004.842871."

## Table of Contents

1. [Introduction](#introduction)
2. [Paper Overview](#paper-overview)
3. [Code Structure](#code-structure)
4. [Cite](#license)

## Introduction

This repository houses the codebase for the digital modeling of networks of biorealistic silicon neurons, as detailed in the research paper "Towards Scalable Digital Modeling of Networks of Biorealistic Silicon Neurons." Our work draws inspiration from the FH neuron model, as presented in the paper "A Bio-Physically Inspired Silicon Neuron."

## Paper Overview

In our paper, we explore the development of scalable digital models for networks of biorealistic silicon neurons. We present detailed biomimetic Spiking Neural Network implementations in multiple programming languages, including C (floating and fixed point), MATLAB, and Python. These models allow for the simulation and study of the behavior of biorealistic silicon neurons in various contexts.

## Code Structure

The repository is organized as follows:

- `C_fixedpoint/`: This directory contains the most efficient implementations of the neuron model, along with various network simulations featuring as many as 100,000 neurons, all utilizing fixed-point operations in the C programming language.
- `C_floatingpoint/`: This directory has the C code  to simulate both neurons and basic Spiking Neural Network (SNN) operations in floating-point format.
- `MATLAB/`: Inside this folder, you'll find MATLAB code used for running SNN simulations.
- `Python/`: In this directory, you can access Python code for the neuron simulations.

## Cite

If you use this simulation framework in your research or educational material, please cite the work as follows: 
Bibtex:
```
@ARTICLE{10305556,
  author={Bhattacharyya, Swagat and Ayyappan, Praveen Raj and Hasler, Jennifer O.},
  journal={IEEE Journal on Emerging and Selected Topics in Circuits and Systems}, 
  title={Towards Scalable Digital Modeling of Networks of Biorealistic Silicon Neurons}, 
  year={2023},
  volume={},
  number={},
  pages={1-1},
  doi={10.1109/JETCAS.2023.3330069}}

}
```

Formatted:
```
S. Bhattacharyya, P. R. Ayyappan and J. O. Hasler, "Towards Scalable Digital Modeling of Networks of Biorealistic Silicon Neurons," in IEEE Journal on Emerging and Selected Topics in Circuits and Systems, doi: 10.1109/JETCAS.2023.3330069.
``` 
The research paper corresponding to the above citation can be found [(https://ieeexplore.ieee.org/document/10305556].
# Small comment