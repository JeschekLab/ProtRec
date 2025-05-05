# uASPIre_TEV

This repository contains the data of _**u**ltradeep **A**cquisition of **S**equence-**P**henotype **I**nter**re**lations_ (uASPIre) and the code for processing raw NGS data as presented in the publication "Data-driven Protease Engineering by DNA-Recording and Epistasis-aware Machine Learning".


## Organization of the repository

1. [**Illumina**](Illumina) contains scripts for processing raw Illumina Novaseq NGS data.
2. [**PacBio**](PacBio) contains scripts for processing raw PacBio SMRT-seq data.
3. [**Illumina_PacBio_combined**](Illumina_PacBio_combined) contains scripts for combining Illumina and PacBio data.


## Hardware requirements
This code package was run on a Unix system with >120 GB RAM, > 1TB of storage and 12 cores.


## Software requirements
The code runs on Unix systems only, and was developed and tested on R version 4.2.1 running on a Red Hat Enterprise Linux Server (release 7.9).


## Software dependencies
This code requires the following UNIX tools and packages:

+ AGREP 3.41.5 (e.g. available via https://github.com/Wikinaut/agrep)
+ TRE agrep 0.8.0 (e.g. available via https://wiki.ubuntuusers.de/tre-agrep/)


## Date
April 2023


## Contact
simon.hoellerer@gmail.com

Simon Höllerer, ETH Zürich, D-BSSE, Basel, Switzerland
