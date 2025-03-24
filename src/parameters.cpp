#include "parameters.h"
#include "print_usage.h"
#include "tools.h"
void parameters::start() {
    additive_variances=3,dominance_variances=1,\
    additive_additive_variances=6,additive_dominance_variances=4,dominance_dominance_variances=2,\
    additive_additive_intra_variances=5,additive_additive_inter_variances=5,additive_dominance_intra_variances=4,\
    additive_dominance_inter_variances=4,dominance_dominance_intra_variances=2,dominance_dominance_inter_variances=2,\
    additive_additive_additive_variances=9,additive_additive_dominance_variances=7,additive_dominance_dominance_variances=5,\
    dominance_dominance_dominance_variances=3,additive_haplotype_variances=2,permanent_env_variances=1,residual_variances=1;

    additive_variances_f=0,dominance_variances_f=0,\
    additive_additive_variances_f=0,additive_dominance_variances_f=0,dominance_dominance_variances_f=0,\
    additive_additive_intra_variances_f=0,additive_dominance_intra_variances_f=0,dominance_dominance_intra_variances_f=0,\
    additive_additive_inter_variances_f=0,additive_dominance_inter_variances_f=0,dominance_dominance_inter_variances_f=0,\
    additive_additive_additive_variances_f=0,additive_additive_dominance_variances_f=0,additive_dominance_dominance_variances_f=0,\
    dominance_dominance_dominance_variances_f=0,additive_haplotype_variances_f=0,permanent_env_variances_f=0,residual_variances_f=0;
    iter_tolerance=1E-8,iter_tolerance_her=1E-6;
    genotype="null",genomap="null",phenotype="null",haplotype="null",outprefix="null",parameter_file="null.txt",grm_file_prefix="null_grm",log_prefix="null",\
    fixed_factor_column="-1",covariate_factor_column="-1";
//    missing_phen_val=“-999999999999”
    missing_phen_val=-9999,missing_hap_val=-9999;
    grm_method=1;
    factors_counts=0,covar_counts=0;
    iter_ai_start=3,iter_n=2000,numThreads=16,num_pairwise_out=30;
    mk_grm_f=false,cin_var_f=false,output_mrk_effect=false,\
    chrsplit_f=false,reml_ce_rel_f=false,rawgeno_f=false,rawmap_f=false,hap_f=false,modelpe_f=false;
    factors_pos_f=false,covar_pos_f=false,loading_grms_f=false,num_pairwise_out_f=false,\
    iter_tolerance_f=false,iter_tolerance_her_f=false;
    return;
}
void parameters::read_cmd_line(int argc,char **argv) {
    if(argc==1){
        parameter_file="parameter.txt";
        cerr<<string(argv[0])<<" only accepts a parameter file named \"parameter.txt\" if you don't specify any file as parameter file."<<endl;
    }
    else if(argc==2){
        if ((strcmp(argv[1],"-h")==0)||(strcmp(argv[1],"--help")==0)){
            print_usage();
            exit(1);
        }
        else{
            parameter_file=string(argv[1]);
            cerr<<string(argv[0])<<" "<<parameter_file<<endl;
        }
    }
    else{
        for (int i=1;i<argc;i++){
            if ((strcmp(argv[i],"-h")==0)||(strcmp(argv[i],"--help")==0)){
                print_usage();
                exit(1);
            }
            if (strcmp(argv[i],"-p")==0){
                parameter_file=string(argv[++i]);
                cerr<<string(argv[0])<<" -p "<<parameter_file<<endl;
            }
            else{
                fprintf(stderr,"\nERROR: please specify parameter_file using '-p'\n");
            }
        }
    }
    ifstream input_par(parameter_file.c_str());
    if(!input_par){
      cerr << "ERROR: can not open the file "<< parameter_file << endl;
      throw;
    }
    else{
        cout<<"reading parameter file " << parameter_file << " ... "<<endl;
    }
    string line;
    string line_org;
    while(getline(input_par,line_org).good()){
        unsigned found = line_org.find('#');
	    if(found==0){
		    line_org.clear();
			continue;
		}
		line = line_org.substr(0,found);
		line_org.clear();
		if(!line.empty()){
            vector<string>words=str_split(" \r\t\n\x0b",line);
            if(words[0] == "geno_snp"){
                if(words.size()<2){
					cerr<<"Error: # geno_snp fields were wrong!"<<endl;
					throw;
				}
				genotype=words[1];
				rawgeno_f=true;
            }
            else if(words[0] == "geno_map"){
                if(words.size()<2){
					cerr<<"Error: # geno_map fields is wrong!"<<endl;
					throw;
				}
				genomap=words[1];
				rawmap_f=true;
            }
            else if(words[0] == "use_geno_hap"){
                if(words.size()<2){
					cerr<<"Error: # geno_hap fields is wrong!"<<endl;
					throw;
				}
                if(words[1] == "Y"){
				    hap_f=true;
				}
				else if(words[1] == "N"){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts 'Y' and 'N' "<<endl;
					throw;
				}
            }
            else if(words[0] == "geno_hap"){
                if(hap_f == true){
                    if(words.size()<2){
					cerr<<"Error: # geno_hap fields is wrong!"<<endl;
					throw;
				    }
				    haplotype=words[1];
                }
                else{
				    continue;
                }
            }
            else if(words[0] == "phenotype"){
                if(words.size()<2){
					cerr<<"Error: # phenotype fields is wrong!"<<endl;
					throw;
				}
				phenotype=words[1];
            }
            else if(words[0] == "missing_phen_val"){
                if(words.size()<2){
					cerr<<"Error: # missing_phen_val fields were wrong!"<<endl;
					throw;
				}
				missing_phen_val=atof(words[1].c_str());
            }
            else if(words[0] == "missing_hap_val"){
                if(words.size()<2){
					cerr<<"Error: # missing_hap_val fields were wrong!"<<endl;
					throw;
				}
				missing_hap_val=atof(words[1].c_str());
            }
            else if(words[0] == "trait_col"){
                if(words.size()<2){
					cerr<<"Error: # trait_col fields were wrong!"<<endl;
					throw;
				}
				trait_pos=atoi(words[1].c_str());
            }
            else if(words[0] == "factors_counts"){
                if(words.size()<2){
					cerr<<"Error: # factors_counts fields were wrong!"<<endl;
					throw;
				}
				factors_counts=atoi(words[1].c_str());
            }
            else if(words[0] == "factors_pos"){
                if (factors_counts>=1){
                    if(words.size()<2){
					    cerr<<"Error: # factors_pos fields were wrong!"<<endl;
					    throw;
				    }
				    stringstream ss;
                    for(int j=1;j<(words.size()-1);j++) ss<<words[j]<<" ";
                    ss<<words[words.size()-1]<<"";
                    fixed_factor_column=ss.str();
				    ss.clear();
				    factors_pos_f=true;
                }
                else{
                    continue;
                }

            }
            else if(words[0] == "covar_counts"){
                if(words.size()<2){
					cerr<<"Error: # covar_counts fields were wrong!"<<endl;
					throw;
				}
				covar_counts=atoi(words[1].c_str());
            }
            else if(words[0] == "covar_pos"){
                if (covar_counts>=1){
                    if(words.size()<2){
					    cerr<<"Error: # covar_pos fields were wrong!"<<endl;
					    throw;
				    }
				    stringstream ss;
                    for(int j=1;j<(words.size()-1);j++) ss<<words[j]<<" ";
                    ss<<words[words.size()-1]<<"";
                    covariate_factor_column=ss.str();
				    ss.clear();
				    covar_pos_f=true;
				}
				else{
				    continue;
				}
            }
            else if(words[0] == "make_grms"){
                if(words.size()<2){
					cerr<<"Error: # make_grms fields were wrong!"<<endl;
					throw;
				}
				if(words[1] == "Y"){
				    mk_grm_f=true;
				}
				else if(words[1] == "N"){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts 'Y' and 'N' "<<endl;
					throw;
				}
            }
            else if(words[0] == "make_partitioned_egrms"){
                if(words.size()<2){
					cerr<<"Error: # make_partitioned_egrms fields were wrong!"<<endl;
					throw;
				}
				if(words[1] == "Y"){
				    chrsplit_f=true;
				}
				else if(words[1] == "N"){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts 'Y' and 'N' "<<endl;
					throw;
				}
            }
            else if(words[0] == "egrms_method"){
                if(words.size()<2){
					cerr<<"Error: # egrms_method fields were wrong!"<<endl;
					throw;
				}
				grm_method=atoi(words[1].c_str());
				if ((grm_method<1||grm_method>2)){
				    cerr<<"Error: # egrms_method only accepts 1 or 2"<<endl;
					throw;
				}
            }
            else if(words[0] == "grm_prefix"){
                if(words.size()<2){
					cerr<<"Error: # grm_prefix fields were wrong!"<<endl;
					throw;
				}
				grm_file_prefix=words[1];
            }
            else if(words[0] == "load_grms"){
                if(words.size()<2){
					cerr<<"Error: # load_grms fields were wrong!"<<endl;
					throw;
				}
				if(words[1] == "Y"){
				    loading_grms_f=true;
				}
				else if(words[1] == "N"){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts 'Y' and 'N' "<<endl;
					throw;
				}
            }
            else if(words[0] == "var_snp_a"){
                if(words.size()<2){
					cerr<<"Error: # var_snp_a fields were wrong!"<<endl;
					throw;
				}
				additive_variances=atof(words[1].c_str());
				if(additive_variances>0){
				    additive_variances_f=1;
				}
				else if(additive_variances<=0){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts double value "<<endl;
					throw;
				}
            }
            else if(words[0] == "var_snp_d"){
                if(words.size()<2){
					cerr<<"Error: # var_snp_d fields were wrong!"<<endl;
					throw;
				}
				dominance_variances=atof(words[1].c_str());
				if(dominance_variances>0){
				    dominance_variances_f=1;
				}
				else if(dominance_variances<=0){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts double value "<<endl;
					throw;
				}
            }
            else if(words[0] == "var_snp_aa"){
                if(words.size()<2){
					cerr<<"Error: # var_snp_aa fields were wrong!"<<endl;
					throw;
				}
				additive_additive_variances=atof(words[1].c_str());
				if(additive_additive_variances>0){
				    additive_additive_variances_f=1;
				}
				else if(additive_additive_variances<=0){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts double value "<<endl;
					throw;
				}
            }
            else if(words[0] == "var_snp_aa-inter"){
                if(words.size()<2){
					cerr<<"Error: # var_snp_aa-inter fields were wrong!"<<endl;
					throw;
				}
				additive_additive_inter_variances=atof(words[1].c_str());
				if(additive_additive_inter_variances>0){
				    additive_additive_inter_variances_f=1;
				}
				else if(additive_additive_inter_variances<=0){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts double value "<<endl;
					throw;
				}
            }
            else if(words[0] == "var_snp_aa-intra"){
                if(words.size()<2){
					cerr<<"Error: # var_snp_aa-intra fields were wrong!"<<endl;
					throw;
				}
				additive_additive_intra_variances=atof(words[1].c_str());
				if(additive_additive_intra_variances>0){
				    additive_additive_intra_variances_f=1;
				}
				else if(additive_additive_intra_variances<=0){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts double value "<<endl;
					throw;
				}
            }
            else if(words[0] == "var_snp_ad"){
                if(words.size()<2){
					cerr<<"Error: # var_snp_ad fields were wrong!"<<endl;
					throw;
				}
				additive_dominance_variances=atof(words[1].c_str());
				if(additive_dominance_variances>0){
				    additive_dominance_variances_f=1;
				}
				else if(additive_dominance_variances<=0){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts double value "<<endl;
					throw;
				}
            }
            else if(words[0] == "var_snp_ad-inter"){
                if(words.size()<2){
					cerr<<"Error: # var_snp_ad-inter fields were wrong!"<<endl;
					throw;
				}
				additive_dominance_inter_variances=atof(words[1].c_str());
				if(additive_dominance_inter_variances>0){
				    additive_dominance_inter_variances_f=1;
				}
				else if(additive_dominance_inter_variances<=0){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts double value "<<endl;
					throw;
				}
            }
            else if(words[0] == "var_snp_ad-intra"){
                if(words.size()<2){
					cerr<<"Error: # var_snp_ad-intra fields were wrong!"<<endl;
					throw;
				}
				additive_dominance_intra_variances=atof(words[1].c_str());
				if(additive_dominance_intra_variances>0){
				    additive_dominance_intra_variances_f=1;
				}
				else if(additive_dominance_intra_variances<=0){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts double value "<<endl;
					throw;
				}
            }
            else if(words[0] == "var_snp_dd"){
                if(words.size()<2){
					cerr<<"Error: # var_snp_dd fields were wrong!"<<endl;
					throw;
				}
				dominance_dominance_variances=atof(words[1].c_str());
				if(dominance_dominance_variances>0){
				    dominance_dominance_variances_f=1;
				}
				else if(dominance_dominance_variances<=0){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts double value "<<endl;
					throw;
				}
            }
            else if(words[0] == "var_snp_dd-inter"){
                if(words.size()<2){
					cerr<<"Error: # var_snp_dd-inter fields were wrong!"<<endl;
					throw;
				}
				dominance_dominance_inter_variances=atof(words[1].c_str());
				if(dominance_dominance_inter_variances>0){
				    dominance_dominance_inter_variances_f=1;
				}
				else if(dominance_dominance_inter_variances<=0){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts double value "<<endl;
					throw;
				}
            }
            else if(words[0] == "var_snp_dd-intra"){
                if(words.size()<2){
					cerr<<"Error: # var_snp_dd-intra fields were wrong!"<<endl;
					throw;
				}
				dominance_dominance_intra_variances=atof(words[1].c_str());
				if(dominance_dominance_intra_variances>0){
				    dominance_dominance_intra_variances_f=1;
				}
				else if(dominance_dominance_intra_variances<=0){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts double value "<<endl;
					throw;
				}
            }
            else if(words[0] == "var_snp_aaa"){
                if(words.size()<2){
					cerr<<"Error: # var_snp_aaa fields were wrong!"<<endl;
					throw;
				}
				additive_additive_additive_variances=atof(words[1].c_str());
				if(additive_additive_additive_variances>0){
				    additive_additive_additive_variances_f=1;
				}
				else if(additive_additive_additive_variances<=0){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts double value "<<endl;
					throw;
				}
            }
            else if(words[0] == "var_snp_aad"){
                if(words.size()<2){
					cerr<<"Error: # var_snp_aad fields were wrong!"<<endl;
					throw;
				}
				additive_additive_dominance_variances=atof(words[1].c_str());
				if(additive_additive_dominance_variances>0){
				    additive_additive_dominance_variances_f=1;
				}
				else if(additive_additive_dominance_variances<=0){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts double value "<<endl;
					throw;
				}
            }
            else if(words[0] == "var_snp_add"){
                if(words.size()<2){
					cerr<<"Error: # var_snp_add fields were wrong!"<<endl;
					throw;
				}
				additive_dominance_dominance_variances=atof(words[1].c_str());
				if(additive_dominance_dominance_variances>0){
				    additive_dominance_dominance_variances_f=1;
				}
				else if(additive_dominance_dominance_variances<=0){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts double value "<<endl;
					throw;
				}
            }
            else if(words[0] == "var_snp_ddd"){
                if(words.size()<2){
					cerr<<"Error: # var_snp_ddd fields were wrong!"<<endl;
					throw;
				}
				dominance_dominance_dominance_variances=atof(words[1].c_str());
				if(dominance_dominance_dominance_variances>0){
				    dominance_dominance_dominance_variances_f=1;
				}
				else if(dominance_dominance_dominance_variances<=0){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts double value "<<endl;
					throw;
				}
            }
            else if(words[0] == "var_hap_a"){
                if(words.size()<2){
					cerr<<"Error: # var_hap_a fields were wrong!"<<endl;
					throw;
				}
				additive_haplotype_variances=atof(words[1].c_str());
				if(additive_haplotype_variances>0){
				    additive_haplotype_variances_f=1;
				}
				else if(additive_haplotype_variances<=0){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts double value "<<endl;
					throw;
				}
            }
            else if(words[0] == "var_e"){
                if(words.size()<2){
					cerr<<"Error: # var_e fields were wrong!"<<endl;
					throw;
				}
				residual_variances=atof(words[1].c_str());
				if(residual_variances<=0){
				    cerr<<"Error: # var_e cannot be less than or equal to 0!"<<endl;
					throw;
				}
				else{
				    residual_variances_f=1;
				}

            }
            else if(words[0] == "num_iter"){
                if(words.size()<2){
					cerr<<"Error: # num_iter fields were wrong!"<<endl;
					throw;
				}
				iter_n=atoi(words[1].c_str());
            }
            else if(words[0] == "ai-reml-iter-start"){
                if(words.size()<2){
					cerr<<"Error: # ai-reml-iter-start fields were wrong!"<<endl;
					throw;
				}
				iter_ai_start=atoi(words[1].c_str());
            }
            else if(words[0] == "tolerance"){
                if(words.size()<2){
					cerr<<"Error: # tolerance fields were wrong!"<<endl;
					throw;
				}
				iter_tolerance=atof(words[1].c_str());
				if(iter_tolerance>0){
				    iter_tolerance_f=1;
				}
				else if(iter_tolerance<=0){
                    cerr<<"Error: # iter_tolerance must be great than 0 "<<endl;
					throw;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts double value "<<endl;
					throw;
				}
            }
            else if(words[0] == "tolerance_her"){
                if(words.size()<2){
					cerr<<"Error: # tolerance_her fields were wrong!"<<endl;
					throw;
				}
				iter_tolerance_her=atof(words[1].c_str());
				if(iter_tolerance_her>0){
				    iter_tolerance_her_f=1;
				}
				else if(iter_tolerance_her<=0){
                    cerr<<"Error: # iter_tolerance_her must be great than 0 "<<endl;
					throw;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts double value "<<endl;
					throw;
				}
            }
            else if(words[0] == "reml-ce-rel"){
                if(words.size()<2){
					cerr<<"Error: # reml-ce-rel fields were wrong!"<<endl;
					throw;
				}
				if(words[1] == "Y"){
				    reml_ce_rel_f=true;
				}
				else if(words[1] == "N"){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts 'Y' and 'N' "<<endl;
					throw;
				}
            }
            else if(words[0] == "marker_effects"){
                if(words.size()<2){
					cerr<<"Error: # marker_effects fields were wrong!"<<endl;
					throw;
				}
				if(words[1] == "Y"){
				    output_mrk_effect=true;
				}
				else if(words[1] == "N"){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts 'Y' and 'N' "<<endl;
					throw;
				}
            }
            else if(words[0] == "cin_var"){
                if(words.size()<2){
					cerr<<"Error: # cin_var fields were wrong!"<<endl;
					throw;
				}
				if(words[1] == "Y"){
				    cin_var_f=true;
				}
				else if(words[1] == "N"){
                    continue;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts 'Y' and 'N' "<<endl;
					throw;
				}
            }
            else if(words[0] == "pairwise_effects"){
                if(words.size()<2){
					cerr<<"Error: # geno_hap fields is wrong!"<<endl;
					throw;
				}
                if(words[1] == "Y"){
				    num_pairwise_out_f=true;
				}
				else if(words[1] == "N"){
                    num_pairwise_out_f=false;
				}
				else{
				    cerr<<"Error: # EPIHAP only accepts 'Y' and 'N' "<<endl;
					throw;
				}
            }
            else if(words[0] == "num_pairwise_out"){
                if(num_pairwise_out_f==true){
                    if(words.size()<2){
					    cerr<<"Error: # num_pairwise_out fields were wrong!"<<endl;
					    throw;
				    }
				    num_pairwise_out=atoi(words[1].c_str());
				}
				else{
				    continue;
                }
            }
            else if(words[0] == "numThreads"){
                if(words.size()<2){
					cerr<<"Error: # numThreads fields were wrong!"<<endl;
					throw;
				}
				numThreads=atoi(words[1].c_str());
            }
            else if(words[0] == "output_gblup_prefix"){
                if(words.size()<2){
					cerr<<"Error: # output_gblup_prefix fields were wrong!"<<endl;
					throw;
				}
				outprefix=words[1];
            }
            else if(words[0] == "log_prefix"){
                if(words.size()<2){
					cerr<<"Error: # log_prefix fields were wrong!"<<endl;
					throw;
				}
				log_prefix=words[1];
            }
            else{
                continue;
            }
		}
    }
    input_par.close();
    if (additive_additive_variances_f==1){
        additive_additive_inter_variances_f=0;
        additive_additive_intra_variances_f=0;
    }
    if (additive_dominance_variances_f==1){
        additive_dominance_inter_variances_f=0;
        additive_dominance_intra_variances_f=0;
    }
    if (dominance_dominance_variances_f==1){
        dominance_dominance_inter_variances_f=0;
        dominance_dominance_intra_variances_f=0;
    }
    if (mk_grm_f==true) loading_grms_f=false;
    return;
}