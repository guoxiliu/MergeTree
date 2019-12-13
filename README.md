# Merge Forest

## Introduction

The program is based on *Toward Localized Topological Data Structures: Querying the Forest for the Tree* by Pavol Klacansky et al. In the paper, the author introduced a localized topological data structure for merge tree, merge forest, to represent topological features.

Our implementation is slightly different from what is proposed in the paper. We make the merge forest as a collection of local structures which contain the merge tree and local bridge set. The merge tree is represented as a list of nodes instead of a set of arcs.

## How to Run

Before running the program, you need to make sure `VTK` and `OpenMP` package are installed and configured correctly in the system. The detailed installation and configuration guide of VTK can be found on the [offical website](https://vtk.org/Wiki/VTK/Configure_and_Build).

The program supports Mac OS X, Linux and Windows. 
- For Unix-based system, please use the provided `Makefile`. 
  - To generate both the serial and parallel program, please use the command `make` or `make all` in the terminal; 
  - To generate the serial program only, please use the command `make serial` in the terminal;
  - To generate the parallel program only, please use the command `make parallel` in the terminal.
- For Windows system, a Visual Studio project file is provided. The project only contains the solution for parallel program, but it is quite straightforward to make another solution for serial program.

