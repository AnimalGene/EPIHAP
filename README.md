[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)<br />
# EPIHAP
This is a computing tool for genomic prediction of quantitative traits using epistasis and haplotype effects.<br />
# Overview
**EPIHAP** is an open-source software program written in C++ for genomic best linear unbiased prediction (GBLUP) methods, which are used for genomic prediction and heritability estimation. By using genome-wide SNPs, this tool can perform multifactorial models to estimate genetic values and heritability for the trait of interest. These models can include any or all types of the following effects: single SNP additive (A) and dominance (D) effects, epistatic effects (including AA,AD,DD,AAA,AAD,ADD,DDD,AA-intra,AA-inter,AD-intra,AD-inter,DD-intra, and DD-inter), and haplotype additive (HA) effects. Two methods can be implemented in EPIHAP to calculate genetic relationship matrices (GRMs) for epistatic effects. One is approximate genomic epistasis relationship matrices (AGERM) and the other exact genomic epistasis relationship matrices (EGERM). By providing property model parameters, EPIHAP also allows estimation of the partitioned pairwise epistatic effects — that are divided into intra-chromosomes and inter-chromosomes. Lastly, the variance components for these types of effects are optimized using genome-based restricted maximum likelihood (GREML) algorithms including EM-REML and AI-REML.
# Installation
The binary executable of EPIHAP is available [here](https://github.com/AnimalGene/EPIHAP/tree/master/bin). This program can only run on Linux system.
# Input Data
To run this tool, you need to prepare the following files:
* Parameter file
* SNP genotype file
* SNP map file
* Haplotype genotype file
* Phenotype file

The example data in the folder `example` is available for demonstration. The `example` contains eight files: `example_parameter_step1.txt`, `example_parameter_step2.txt`,`example_parameter_step2_for_A+D+AA-intra+AA-inter.txt`,`example_parameter_for_marker_effects.txt`,`example.dat`, `example.map`, `example.hap`, and `example.phen`. The `example_parameter_step1.txt`, `example_parameter_step2.txt`,`example_parameter_step2_for_A+D+AA-intra+AA-inter.txt`, and `example_parameter_for_marker_effects.txt` are the parameter files. The parameter file stores all the parameters that can be read by EPIHAP. The file is separated into sections, each containing default values as provided. The genotypic data is stored in the `example.dat` file, which includes a header line followed by 100 lines, each representing an individual with data for 300 SNPs. The `example.map` file contains SNP information, consisting of a header line followed by 300 lines, each corresponding to a specific SNP. Haplotype genotypes are stored in the `example.hap` file, which includes a header line followed by 100 lines, each representing an individual with data for 110 haplotype blocks (with two columns per block). Finally, the phenotypic data is stored in the `example.phen` file, which contains a header line followed by 100 lines, each corresponding to an individual’s phenotypic information.
# Quick Start
EPIHAP can be run on a terminal using command line with the parameter file name as an argument:
```console
./EPIHAP parameter.txt
```
EPIHAP can run without any input if a parameter file named `parameter.txt` exists in the current working directory. Thus, the command mentioned above can also be typed as follows:
```console
./EPIHAP
```
Below is an example of how to run EPIHAP for genomic prediction in two steps:
The first step of genomic prediction is making genomic relationship matrices using SNP genotypes or haplotypes data. To execute this program, run for example with the following command:
```console
./EPIHAP example_parameter_step1.txt
```
The second step is to estimate the variance components or heritability, genetic values and reliability for all types of genetic effects. An example run for the model A+D+AA is as follows:
```console
./EPIHAP example_parameter_step2.txt
```
EPIHAP also can estimate the genetic values and reliability for three partitioned pairwise epistatic effects. Before proceeding with this task, the **make_partitioned_egrms** parameter in `example_parameter_step2.txt` should be set to **Y**, and the starting values of the **var_snp_aa**, **var_snp_ad**, and **var_snp_dd** parameters, which must be skipped by EPIHAP, are set to be less than or equal to 0. For example, if we decide to run model A+D+AA-intra+AA-inter,we should use the following command:
```console
./EPIHAP example_parameter_step2_for_A+D+AA-intra+AA-inter.txt
```
In addition, EPIHAP also can estimate the genetic effects and heritability for a SNP, a pair of SNPs or a haplotype block after estimating the variance components using GREML method. We here give an example for the model A+D+AA+AD+DD+AH to show how to calculate these marker effects. To do this, the **marker_effects**, **pairwise_effects**, and **cin_var** parameters in `example_parameter_for_marker_effects.txt` should be set to **Y**. An example run is as follows:
```console
./EPIHAP example_parameter_for_marker_effects.txt
```
**Detailed instructions for running the program are available in the [user manual](https://github.com/AnimalGene/EPIHAP/blob/master/doc/EPIHAP_user_manual6.pdf).**
## License
It is a free and open source software, licensed under [GPLv3 LICENSE](https://github.com/Leon-Liang591/EPIHAP/blob/main/LICENSE).
