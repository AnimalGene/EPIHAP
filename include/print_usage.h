#include <iostream>
#include <iomanip>
using namespace std;
#ifndef __PRINT_USAGE_H
#define __PRINT_USAGE_H
void print_usage(){
    cerr <<"   EPIHAP Program Usage:\n\n";
    cerr <<"             --help or -h	Print help statement.\n\n";
    cerr << "# The lines starting with '#' sign are comments for the parameter definitions, which cannot be read by EPIHAP. Only the lines starting with the parameter name can be delivered to EPIHAP. The parameter name and parameter values in each non-comment line are delimited by a whitespace ' '.\n";
    cerr << "# The capital letters 'Y' or 'N' observed in some non-comment lines are used to specify whether some certain parameters should be passed to EPIHAP or not (Y=Yes, N=No).\n";
    cerr << "## ###################################################################################\n";
    cerr << R"(
    # Specify the full path to the filename of the genotype file.
    geno_snp /home/path/to/example.dat
    # Specify the full path to the filename of the SNP map file.
    geno_map /home/path/to/example.map
    # Set Y/N to specify whether the parameter geno_hap should be passed to EPIHAP or not.
    use_geno_hap Y
    # Specify the full path to the filename of the haplotype genotypes file.
    geno_hap /home/path/to/example.hap
    # Specify the full path to the filename of the phenotype file.
    phenotype /home/path/to/example.phen
    # Set the missing phenotype values [default=-9999]. This value must be same as the missing values occurred in phenotype file.
    missing_phen_val -9999111
    # Set the missing haplotype values [default=-9999]. This value must be same as the missing values occurred in haplotype file.
    missing_hap_val -9999111
    # Set the position of the desired trait of interest in the phenotype file.
    trait_col 7
    # Set the number of the fixed non-genetic factors only for discrete variables in the phenotype file. The parameter factors_pos will be skipped if the integer number is set to less than 1.
    factors_counts 2
    # Set the positions of the fixed non-genetic factors only for discrete variables in the phenotype file if the parameter factors_counts ≥ 1.
    factors_pos 2 3
    # Set the number of the covariables in the phenotype file. The parameter covar_pos will be skipped if the integer number is set to less than 1.
    covar_counts 2
    # Set the positions of the covariables in the phenotype file.
    covar_pos 4 5 6
    # Set Y/N to turn on/off running the construction of GRMs.
    make_grms Y
    # Set Y/N to turn on/off running the construction of GRMs for the pairwise epistatic effects (AA, AD and DD) that are partitioned into intra-chromosomes and inter-chromosomes.
    make_partitioned_egrms N
    # Set the method to construct epistatic GRMs via integers 1 or 2. 1: Approximate Genomic Epistasis Relationship Matrices (AGERM); 2: Exact Genomic Epistasis Relationship Matrices (EGERM), [default=1].
    egrms_method 1
    # Specify the full path to the prefix of GRM file names.
    grm_prefix /home/path/to/grmfiles/example
    # Set Y/N to turn on/off loading GRMs and executing variance components estimation using GREML method and calculating heritability estimates and genetic values.
    load_grms N
    # Set the starting value of the additive variance component. The user can set a starting value less than or equal to 0 (var_snp_a ≤ 0) to skip this parameter.
    var_snp_a 3
    # Set the starting value of the dominance variance component.
    var_snp_d 1
    # Set the starting value of the additive × additive variance component.
    var_snp_aa 6
    # Set the starting value of the additive × additive variance component only for inter-chromosomes.
    var_snp_aa-inter 0
    # Set the starting value of the additive × additive variance component only for intra-chromosomes.
    var_snp_aa-intra 0
    # Set the starting value of the additive × dominance variance component.
    var_snp_ad 4
    # Set the starting value of the additive × dominance variance component only for inter-chromosomes.
    var_snp_ad-inter 0
    # Set the starting value of the additive × dominance variance component only for intra-chromosomes.
    var_snp_ad-intra 0
    # Set the starting value of the dominance × dominance variance component.
    var_snp_dd 2
    # Set the starting value of the dominance × dominance variance component only for inter-chromosomes.
    var_snp_dd-inter 0
    # Set the starting value of the dominance × dominance variance component only for intra-chromosomes.
    var_snp_dd-intra 0
    # Set the starting value of the additive × additive × additive variance component.
    var_snp_aaa 9
    # Set the starting value of the additive × additive × dominance variance component.
    var_snp_aad 7
    # Set the starting value of the additive × dominance × dominance variance component.
    var_snp_add 5
    # Set the starting value of the dominance × dominance × dominance variance component.
    var_snp_ddd 3
    # Set the starting value of the haplotype additive variance component.
    var_hap_a 3
    # Set the starting value of the residual variance.
    var_e 1
    # Set the maximum number of iterations that are allowed in the GREML_CE method [default=1000].
    num_iter 1000
    # Set the starting iteration number for converting the EM-REML to the AI-REML to estimate the variance components [default=3].
    ai-reml-iter-start 3
    # Set the tolerance threshold as a convergence criterion to stop estimating the variance components [default=1E-8].
    tolerance 1.0E-08
    # Set the tolerance threshold as a convergence criterion to stop estimating the heritability [default=1E-6].
    tolerance_her 1.0E-06
    # Set Y/N to turn on/off calculation of the reliability of genomic breeding values.
    reml-ce-rel N
    # Set Y/N to turn on/off computing the genetic effects and heritability estimates partitioned by SNPs or haplotype blocks.
    marker_effects N
    # Set Y/N to turn on/off computing the pairwise epistatic effects and heritability estimates that are partitioned by SNP pairs.
    pairwise_effects N
    # Set the number of top-ranked SNP pairs with highest pairwise epistatic effects and heritability estimates [default = 30].
    num_pairwise_out 30
    # Set Y/N to turn on/off calculating the genetic values with the specified variance components.
    cin_var N
    # Specify the full path to the prefix filenames of main GBLUP output files.
    output_gblup_prefix /home/path/to/out/example
    # Set the number of threads for parallel computing [default=16].
    numThreads 16
    # Set the full path to the prefix filename of log file.
     log_prefix /home/path/to/log/example
 )";
}
#endif