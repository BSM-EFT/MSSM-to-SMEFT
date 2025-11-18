# MSSM-to-SMEFT

[![DOI](https://zenodo.org/badge/946582313.svg)](https://doi.org/10.5281/zenodo.15602469) 
[![arXiv](https://img.shields.io/badge/arXiv-2506.05201-00aa00.svg)](https://arxiv.org/abs/2506.05201)

This repository contains the supplementary material for the paper "SUSY meets SMEFT: Complete one-loop matching of the general MSSM", [arXiv:2506.05201](https://arxiv.org/abs/2506.05201), by Sabine Kraml, Andre Lessa, Suraj Prakash and Felix Wilsch. 
Concretely we provide:

- the implementation of the MSSM in Matchete (see folder [matching](matching)); 
- the full matching results in the form of Mathematica and PDF files  (folder [matching-results](matching-results));
- example Mathematica notebooks showing how to use these results (folder [examples](examples));
- setup for matching the MSSM onto the 2 Higgs doublet model EFT (see folder [2HDM-EFT](2HDM-EFT)), which contains:
    * the Matchete implementation of the MSSM with two light Higg;
    * the Matchete implementation of a [2HDM-EFT basis](https://arxiv.org/abs/2405.20511);
    * full matching results for this scenario in the form of Mathematica and PDF files.

The matching is performed up to dimension six at a single scale $\bar\mu$, with the masses of all BSM states kept generic and non-degenerate.  
*Note that the code provided here requires [Matchete v0.4.0](https://gitlab.com/matchete/matchete/tree/v0.4.0) or higher!*  
