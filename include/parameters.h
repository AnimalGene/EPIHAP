#ifndef __PARAMETERS_H
#define __PARAMETERS_H
#include <iostream>
#include <cstring>
#include <string>
#include <vector>
#include <list>
using namespace std;
/// command line information and global parameters
class parameters{
public:
    int iter_n,iter_ai_start,num_pairwise_out,grm_method,numThreads,factors_counts,covar_counts,\
    trait_pos,additive_variances_f,dominance_variances_f,\
    additive_additive_variances_f,additive_dominance_variances_f,dominance_dominance_variances_f,\
    additive_additive_intra_variances_f,additive_dominance_intra_variances_f,dominance_dominance_intra_variances_f,\
    additive_additive_inter_variances_f,additive_dominance_inter_variances_f,dominance_dominance_inter_variances_f,\
    additive_additive_additive_variances_f,additive_additive_dominance_variances_f,additive_dominance_dominance_variances_f,\
    dominance_dominance_dominance_variances_f,additive_haplotype_variances_f,permanent_env_variances_f,residual_variances_f;
    double missing_phen_val,missing_hap_val,additive_variances,dominance_variances,\
    additive_additive_variances,additive_dominance_variances,dominance_dominance_variances,\
    additive_additive_intra_variances,additive_additive_inter_variances,additive_dominance_intra_variances,\
    additive_dominance_inter_variances,dominance_dominance_intra_variances,dominance_dominance_inter_variances,\
    additive_additive_additive_variances,additive_additive_dominance_variances,\
    additive_dominance_dominance_variances,dominance_dominance_dominance_variances,\
    additive_haplotype_variances,permanent_env_variances,residual_variances,iter_tolerance,iter_tolerance_her;
    /// target trait column in phenotype file
    string fixed_factor_column,covariate_factor_column,genotype,parameter_file,genomap,phenotype,haplotype,load_file_prefix,grm_file_prefix,log_prefix,outprefix;
    bool mk_grm_f,factors_pos_f,covar_pos_f,cin_var_f,loading_grms_f,num_pairwise_out_f,\
    output_mrk_effect,reml_ce_rel_f,rawgeno_f,rawmap_f,hap_f,modelpe_f,chrsplit_f,iter_tolerance_f,iter_tolerance_her_f;
    /// read relevant information
    void start ();
    void read_cmd_line (int argc,char **argv);
};
#endif