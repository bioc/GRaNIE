<img src="man/figures/logo.png" align="right" width="250"/>

## GRaNIE: Reconstruction and evaluation of data-driven, cell type specific gene regulatory networks including enhancers using chromatin accessibility and RNAseq data


### Repository structure

This repository is for our *GRaNIE* package (**G**ene **R**egul**a**tory **N**etwork **I**nference including **E**nhancers). The accompanying [website](https://grp-zaugg.embl-community.io/GRaNIE) is **not** part of Bioconductor and independently build and maintained based on the development version of *GRaNIE*, which is hosted [in this repository](https://git.embl.de/grp-zaugg/GRaNIE). It may therefore differ from the *release version* of the package.

For the *GRaNPA* package, see [here](https://grp-zaugg.embl-community.io/GRaNPA).

### Installation and Documentation

*GRaNIE* is now included on *Bioconductor*. The package and installation instructions can be found here and is suitable for the most users, which installs the package from the release branch of Bioconductor: [https://bioconductor.org/packages/GRaNIE](https://bioconductor.org/packages/GRaNIE)

**The package is still under very active development. We recommend to either use the development version of Bioconductor for being able to isntall the latest versions OR to install the newest version of the package outside of the official Bioconductor release cycle (i.e., a devel version that however does not necessarily require even the latest Bioconductor version to be installed)**

For the latter option, you may use the following line that will install the package directly from the Gitlab instance that also hosts the package website. Due to potential package conflicts with existing (Bioconductor) packages, we cannot guarantee this always works, though:
`devtools::install_gitlab("grp-zaugg/GRaNIE", host = "git.embl.de", subdir = "src/GRaNIE", force = TRUE)`. If you run into any problems, let us know!


In addition, **the full documentation for the development version is available at [https://grp-zaugg.embl-community.io/GRaNIE](https://grp-zaugg.embl-community.io/GRaNIE) and regularly updated and extended**.


### Citation
**If you use our packages, please use the following citation (we will update it once the paper has been officially published):**

GRaNIE and GRaNPA: Inference and evaluation of enhancer-mediated gene regulatory networks applied to study macrophages. Aryan Kamal, Christian Arnold, Annique Claringbould, Rim Moussa, Neha Daga, Daria Nogina, Maksim Kholmatov, Nila Servaas, Sophia Mueller-Dott, Armando Reyes-Palomares, Giovanni Palla, Olga Sigalova, Daria Bunina, Caroline Pabst, Judith B. Zaugg. bioRxiv 2021.12.18.473290; doi: https://doi.org/10.1101/2021.12.18.473290

### Bug Reports, Feature Requests and Contact Information

For issues, bugs, and feature request, please see the [Issue Tracker](https://git.embl.de/grp-zaugg/GRaNIE/issues). 

**We are actively working on the package and regularly improve upon features, add features, or change features for increased clarity. This sometimes results in minor changes to the workflow, changed argument names or other small incompatibilities that may result in errors when running a version of the package that differs from the version this vignette has been run for.**
**Thus, make sure to run a version of `GRaNIE` that is compatible with this vignette. If in doubt or when you receive errors, check the R help, which always contains the most up-to-date documentation.**

If you have other questions or comments, feel free to contact us. We will be happy to answer any questions related to this project as well as questions related to the software implementation. For method-related questions, contact Judith B. Zaugg (judith.zaugg@embl.de). For technical questions, contact Christian Arnold (christian.arnold@embl.de). We will aim to respond in a timely manner.

 


