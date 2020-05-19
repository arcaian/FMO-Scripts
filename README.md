# FMO-Scripts

This repository contains a series of (primarily Python) scripts for use in the running of GAMESS-US FMO jobs, or jobs derived from these. A short summary of the scripts follows:

## extractFragmentCoord.py
This file is a Python script that takes a GAMESS-US FMO job as an input, with a file extension of .inp, and provides the .xyz coordinates of each fragment. Additionally, a selection of fragments are written out as individual Psi4 input files, and a joint Psi4 input files with the coordinates of each of these fragments together. This enables the running of specific dimers or trimers in a non-FMO context, if issues are occurring with these individual FMO fragment combinations.
