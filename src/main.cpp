#include <iostream>
#include <vector>
#include <map>
#include <time.h>
#include <string>
#include <fstream>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "read_geno.h"
#include <cstdlib>
using namespace std;
using namespace Eigen;
int main (int argc, char **argv){
	clock_t start_t,finish_t;
	time_t now = time(0);
	char* current_time = ctime(&now);
    parameters params;
    params.start();
	params.read_cmd_line(argc,argv);
	string log_filename;
    params.mk_grm_f==true ? log_filename=params.log_prefix+"_for_making_grms.log":log_filename=params.log_prefix+"_for_gblup.log";
	ofstream flog(log_filename.c_str());
	flog<<"***********************EPIHAP log file****************************\n"<<endl;
	flog<<"Version                                            :\t1.0.0"<<endl;
	flog<<"The local date and time is                         :\t"<<current_time<<endl;
	flog<<"Reading parameter file " << params.parameter_file << " ... \n"<<endl;
	flog<<"geno_snp "<<params.genotype<<endl;
	flog<<"geno_map "<<params.genomap<<endl;
	flog<<"use_geno_hap "<<(params.hap_f==true ? "Y ":"N ")<<endl;
	flog<<"geno_hap "<<params.haplotype<<endl;
	flog<<"phenotype "<<params.phenotype<<endl;
	flog<<"missing_phen_val "<<params.missing_phen_val<<endl;
	flog<<"trait_col "<<params.trait_pos<<endl;
	flog<<"factors_counts "<<params.factors_counts<<endl;
	flog<<"factors_pos "<<params.fixed_factor_column<<endl;
	flog<<"covar_counts "<<params.covar_counts<<endl;
	flog<<"covar_pos "<<params.covariate_factor_column<<endl;
	flog<<"make_grms "<<(params.mk_grm_f==true ? "Y":"N")<<endl;
	flog<<"make_partitioned_egrms "<<(params.chrsplit_f==true ? "Y":"N")<<endl;
	flog<<"egrms_method "<<params.grm_method<<endl;
	flog<<"grm_prefix "<<params.grm_file_prefix<<endl;
	flog<<"load_grms "<<(params.loading_grms_f==true ? "Y":"N")<<endl;
	flog<<"var_snp_a "<<params.additive_variances<<endl;
	flog<<"var_snp_d "<<params.dominance_variances<<endl;
	flog<<"var_snp_aa "<<params.additive_additive_variances<<endl;
	flog<<"var_snp_aa-inter "<<params.additive_additive_inter_variances<<endl;
	flog<<"var_snp_aa-intra "<<params.additive_additive_intra_variances<<endl;
	flog<<"var_snp_ad "<<params.additive_dominance_variances<<endl;
	flog<<"var_snp_ad-inter "<<params.additive_dominance_inter_variances<<endl;
	flog<<"var_snp_ad-intra "<<params.additive_dominance_intra_variances<<endl;
	flog<<"var_snp_dd "<<params.dominance_dominance_variances<<endl;
	flog<<"var_snp_dd-inter "<<params.dominance_dominance_inter_variances<<endl;
	flog<<"var_snp_dd-intra "<<params.dominance_dominance_intra_variances<<endl;
	flog<<"var_snp_aaa "<<params.additive_additive_additive_variances<<endl;
	flog<<"var_snp_aad "<<params.additive_additive_dominance_variances<<endl;
	flog<<"var_snp_add "<<params.additive_dominance_dominance_variances<<endl;
	flog<<"var_snp_ddd "<<params.dominance_dominance_dominance_variances<<endl;
	flog<<"var_hap_a "<<params.additive_haplotype_variances<<endl;
	flog<<"var_e "<<params.residual_variances<<endl;
	flog<<"num_iter "<<params.iter_n<<endl;
	flog<<"ai-reml-iter-start "<<params.iter_ai_start<<endl;
    flog<<"tolerance "<<params.iter_tolerance<<endl;
    flog<<"tolerance_her "<<params.iter_tolerance_her<<endl;
    flog<<"reml-ce-rel "<<(params.reml_ce_rel_f==true ? "Y":"N")<<endl;
    flog<<"marker_effects "<<(params.output_mrk_effect==true ? "Y":"N")<<endl;
    flog<<"pairwise_effects "<<(params.num_pairwise_out_f==true ? "Y":"N")<<endl;
    flog<<"num_pairwise_out "<<params.num_pairwise_out<<endl;
    flog<<"cin_var "<<(params.cin_var_f==true ? "Y":"N")<<endl;
    flog<<"output_gblup_prefix "<<params.outprefix<<endl;
	flog<<"numThreads "<<params.numThreads<<endl;
	flog<<"log_prefix "<<params.log_prefix<<endl<<endl;
	Greml_ce greml_ce;
	greml_ce.init(params);
	start_t = clock();
	if (params.mk_grm_f) {
	    if(params.grm_method==1){
	        if(params.rawgeno_f==true) greml_ce.read_genotype_method_1(params,flog);
	        if(params.hap_f==true) greml_ce.read_hap(params,flog);
	        fprintf(stderr,"\nApproximate Genomic Epistasis Relationship Matrices (AGERM) method is used for GRMs inference!\n");
	    }
	    else if(params.grm_method==2){
	        if(params.rawgeno_f==true) greml_ce.read_genotype_method_2(params,flog);
	        if(params.hap_f==true) greml_ce.read_hap(params,flog);
	        fprintf(stderr,"\nExact Genomic Epistasis Relationship Matrices (EGERM) method is used for GRMs inference!\n");
	    }
	    else{
	        ;
	    }
	}
	else if (params.loading_grms_f) {
	    greml_ce.load_genotype(params,flog);
	    greml_ce.read_phen_pre(params,flog);
	    greml_ce.read_phen(params);
	    greml_ce.reml(params,flog);
	}
	else{
	;
	}
	finish_t = clock();
    if(params.mk_grm_f){
        flog<<"\nThe time (seconds) cost for all GRMs inference is  :\t"<<(int)((double)(finish_t - start_t)/CLOCKS_PER_SEC)<<endl;
    }
    else{
        flog<<"\nThe time (seconds) cost for GREML_CE is            :\t"<<(int)((double)(finish_t - start_t)/CLOCKS_PER_SEC)<<endl;
    }
    flog.close();
    log_filename.clear();
	return 0;
}