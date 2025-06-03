# Streptomycin Resistance Landscape
This repository contains the data, code, and analysis for a project investigating streptomycin (STR) resistance across the bacterial tree of life, focusing on resistance mutations in the rpsL and rrs genes.

## Overview
This project aims to assess both intrinsic resistance—where streptomycin resistance mutations are already present in bacterial genomes—and the evolutionary potential for bacteria to acquire resistance through single-nucleotide mutations. Using a computational framework adapted from previous work on rifampicin [^1], the analysis screens a curated panel of resistance-associated mutations in the rpsL and rrs genes to map resistance patterns across the bacterial tree of life and quantify species-specific evolvability.


## Contents

- R/: Contains all analysis scripts, including:

  - master.R: Main pipeline for analyzing rpsL mutations.

  - master_nt.R: Main pipeline for analyzing rrs mutations.

- output/: Intermediate analysis outputs, such as filtered sequences and mutation screens.

- plots/: Final figures generated for the report or thesis.

- results/: Summary tables and final processed results.


## Citation

[^1] Bolourchi, N., Brown, C. R. P., Letten, A. D., & Engelstädter, J. (2024). Evolution and evolvability of rifampicin resistance across the bacterial tree of life. bioRxiv. 
[![DOI:10.1101/2024.11.05.622190](https://zenodo.org/badge/DOI/10.1101/2024.11.05.622190.svg)](https://doi.org/10.1101/2024.11.05.622190)


