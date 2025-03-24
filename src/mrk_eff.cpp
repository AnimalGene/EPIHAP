#include "read_geno.h"
void Greml_ce::insert_epistasis_effect(const parameters &p,Pairwise_result_one *epistasis_effects,int &n_inserted,Pairwise_result_one epistasis_effect){
    int i,n,imin, imax, imid, key;
    if(n_inserted<(p.num_pairwise_out-1)) n_inserted+=1;
    imin=0;
    imax=n_inserted;
    key=-1;
    if (epistasis_effect.effect<epistasis_effects[0].effect){
        key=0;
    }
    else{
        while(key==-1){
            imid=imin+((imax-imin)/2);
            if(epistasis_effects[imid].effect == epistasis_effect.effect){
                key=imid;
                break;
            }
            else if(epistasis_effects[imid].effect < epistasis_effect.effect){
                imin=imid;
            }
            else{
                imax=imid;
            }
            if(imax-imin<=1) key=imax;
        }
    }
    n=n_inserted-key+1;
    vector<Pairwise_result_one>epistasis_effects_vec(n);
    for(i=0;i<n;i++) epistasis_effects_vec[i]=epistasis_effects[key+i];
    rotate(epistasis_effects_vec.rbegin(),epistasis_effects_vec.rbegin()+1,epistasis_effects_vec.rend());
    epistasis_effects_vec[0]=epistasis_effect;
    for(i=0;i<n;i++) epistasis_effects[key+i]=epistasis_effects_vec[i];
    epistasis_effects_vec.clear();
}
void Greml_ce::calc_each_epi_eff(const parameters &p,MatrixXd &mrk_epi_eff,string eff_type,Pairwise_result_one *epistasis_effects){
    int j1,j2,n_inserted;
    double neg_abs_pairwise_value;
    n_inserted=-1;
    for(j1=0;j1<mrk_n;j1++){
        for(j2=1;j2<mrk_n;j2++){
            if(j2>j1){
                neg_abs_pairwise_value=(-1)*abs(mrk_epi_eff(j1,j2));
                if(neg_abs_pairwise_value<epistasis_effects[p.num_pairwise_out-1].effect){
                    Pairwise_result_one only_one_result;
                    only_one_result.effect_type=eff_type;
                    only_one_result.snp1=j1;
                    only_one_result.snp2=j2;
                    only_one_result.real_value=mrk_epi_eff(j1,j2);
                    only_one_result.effect=neg_abs_pairwise_value;
                    insert_epistasis_effect(p,epistasis_effects,n_inserted,only_one_result);
                }
            }
        }
    }
    cout<<"The top pair snps were been calculated!"<<endl;
}
void Greml_ce::get_SNP_eff(const parameters &p,VectorXd &r){
    time_t time_a,time_b,time_a_geno,time_b_geno,time_a_w;
    time(&time_a);
    int i,j,have_single_var,have_epi_var;
    if (p.rawmap_f==true){
        fprintf(stderr,"Map file was found!\n");
        string rawmapFile=p.genomap;
        mrk_n=-1;
        read_map(rawmapFile);
    }
    else{
        fprintf(stderr,"ERROR: The input genotype file!\n");
    }
    MatrixXd WT;
    VectorXd mrk_eff_tmp;
    vector<VectorXd> mrk_eff_tot;
    have_single_var=0;
    have_epi_var=0;
    for (i=0;i<var_name.size()-1;i++){
        if((var_name[i]=="A")||(var_name[i]=="D")){
            have_single_var++;
        }
        else if((var_name[i]=="AA")||(var_name[i]=="AD")||(var_name[i]=="DD")){
            have_epi_var++;
        }
        else{
            ;
        }
    }
    if(have_single_var>0){
        for (i=0;i<var_name.size()-1;i++){
            if((var_name[i]=="A")||(var_name[i]=="D")){
//                string gmat_filename =p.load_file_prefix+".g.W"+var_name[i]+"T";
//                read_binary(gmat_filename.c_str(),WT);
                string gmat_filename =p.grm_file_prefix+".g.W"+var_name[i]+"T.txt";
                read_plain_txt(gmat_filename.c_str(),WT);
                WT/=G_diag_mean_sqrt(i);
                mrk_eff_tmp=WT*r*var(i);
                mrk_eff_tot.push_back(mrk_eff_tmp);
            }
        }
        mrk_eff_tmp.resize(0);
        get_single_mrk_eff(p,mrk_eff_tot);
        for (i=0;i<mrk_eff_tot.size();i++) mrk_eff_tot[i].resize(0);
        WT.resize(0,0);
    }
    if(have_epi_var>0 && p.num_pairwise_out_f==true){
        vector<string>var_epi_name={"AA","AD","DA","DD"};
        int j1,j2_size,j2,k;
        mrk_epi_n=(mrk_n*(mrk_n-1))/2;
        string wat_filename,wdt_filename;
        wat_filename=p.grm_file_prefix+".g.WAT.txt";
        wdt_filename=p.grm_file_prefix+".g.WDT.txt";
        MatrixXd mrk_2_order_epi_eff;
        Pairwise_result_one myinition={var_epi_name[i],0,0,0.0000,1.00000};
        for (i=0;i<var_name.size()-1;i++){
            if(var_name[i]=="AD"){
                cout<<"******************AD effect***********"<<endl;
                VectorXd her_mean_std_all=VectorXd::Zero(3);
                mrk_2_order_epi_eff=get_epi_eff_each_type(wat_filename,wdt_filename,var_name[i],r);
                cout<<mrk_2_order_epi_eff.block(0,0,5,5)<<endl;
                cout<<"Starting to compute the mean and std of heritabilities!"<<endl;
                get_epi_mrk_her_m_std(p,mrk_2_order_epi_eff,her_mean_std_all,i);
                Pairwise_result_one *pairwise_results_ad=new Pairwise_result_one[p.num_pairwise_out];
                for(j=0;j<p.num_pairwise_out;j++) pairwise_results_ad[j]=myinition;
                calc_each_epi_eff(p,mrk_2_order_epi_eff,var_name[i],pairwise_results_ad);
                cout<<"******************DA effect***********"<<endl;
                MatrixXd mrk_2_order_epi_eff_da;
                VectorXd her_mean_std_all_da=VectorXd::Zero(3);
                mrk_2_order_epi_eff_da=get_epi_eff_each_type(wdt_filename,wat_filename,var_name[i],r);
                cout<<mrk_2_order_epi_eff_da.block(0,0,5,5)<<endl;
                cout<<"Starting to compute the mean and std of heritabilities!"<<endl;
                get_epi_mrk_her_m_std(p,mrk_2_order_epi_eff_da,her_mean_std_all_da,i);
                cout<<"******************AD and DA effect (sorted by AD)***********"<<endl;
                get_epi_mrk_eff_AD(p,var_name[i],mrk_2_order_epi_eff,mrk_2_order_epi_eff_da,her_mean_std_all,her_mean_std_all_da,pairwise_results_ad);
                delete [] pairwise_results_ad;
                her_mean_std_all.resize(0);
                her_mean_std_all_da.resize(0);
                mrk_2_order_epi_eff_da.resize(0,0);
                mrk_2_order_epi_eff.resize(0,0);
            }
            else if(var_name[i]=="AA"){
                cout<<"******************AA effect***********"<<endl;
                VectorXd her_mean_std_all=VectorXd::Zero(3);
                mrk_2_order_epi_eff=get_epi_eff_each_type(wat_filename,wat_filename,var_name[i],r);
                cout<<mrk_2_order_epi_eff.block(0,0,5,5)<<endl;
                cout<<"Starting to compute the mean and std of heritabilities!"<<endl;
                get_epi_mrk_her_m_std(p,mrk_2_order_epi_eff,her_mean_std_all,i);
                Pairwise_result_one *pairwise_results_aa=new Pairwise_result_one[p.num_pairwise_out];
                for(j=0;j<p.num_pairwise_out;j++) pairwise_results_aa[j]=myinition;
                calc_each_epi_eff(p,mrk_2_order_epi_eff,var_name[i],pairwise_results_aa);
                get_epi_mrk_eff_AA_or_DD(p,var_name[i],mrk_2_order_epi_eff,her_mean_std_all,pairwise_results_aa);
                delete [] pairwise_results_aa;
                her_mean_std_all.resize(0);
                mrk_2_order_epi_eff.resize(0,0);
            }
            else if(var_name[i]=="DD"){
                cout<<"******************DD effect***********"<<endl;
                VectorXd her_mean_std_all=VectorXd::Zero(3);
                mrk_2_order_epi_eff=get_epi_eff_each_type(wdt_filename,wdt_filename,var_name[i],r);
                cout<<mrk_2_order_epi_eff.block(0,0,5,5)<<endl;
                cout<<"Starting to compute the mean and std of heritabilities!"<<endl;
                get_epi_mrk_her_m_std(p,mrk_2_order_epi_eff,her_mean_std_all,i);
                Pairwise_result_one *pairwise_results_dd=new Pairwise_result_one[p.num_pairwise_out];
                for(j=0;j<p.num_pairwise_out;j++) pairwise_results_dd[j]=myinition;
                calc_each_epi_eff(p,mrk_2_order_epi_eff,var_name[i],pairwise_results_dd);
                get_epi_mrk_eff_AA_or_DD(p,var_name[i],mrk_2_order_epi_eff,her_mean_std_all,pairwise_results_dd);
                delete [] pairwise_results_dd;
                her_mean_std_all.resize(0);
                mrk_2_order_epi_eff.resize(0,0);
            }
            else{
                ;
            }
        }
    }
    if(p.haplotype!="null") get_mrk_eff_hap(p,r);
}
MatrixXd Greml_ce::get_epi_eff_each_type(string &W1Tfile,string &W2Tfile,string var_epi_type,VectorXd &r){
    MatrixXd mrk_epi_eff;
    double alpha = 1.0,beta  = 0.0;
    int i,j,k,column_w;
    double var_val,gdiag_mean_sqrt;
    for (i=0;i<var_name.size()-1;i++){
        if(var_name[i]==var_epi_type){
            var_val=var(i);
            gdiag_mean_sqrt=G_diag_mean_sqrt(i);
        }
    }
    time_t time_epi_s,time_epi_e;
    char* w1file;
    char* w2file;
    w1file=&W1Tfile[0];
    w2file=&W2Tfile[0];
    matrix_struct *m_1 = get_matrix_struct(w1file);
    matrix_struct *m_2 = get_matrix_struct(w2file);
    matrix_struct *result_matrix = (matrix_struct*)malloc(sizeof(matrix_struct));
    result_matrix->rows = m_1->rows;
    result_matrix->cols = m_2->rows;
    result_matrix->mat_data = (double**)calloc(result_matrix->rows, sizeof(double*));
    for(i=0; i < result_matrix->rows; ++i) result_matrix->mat_data[i]=(double*)calloc(result_matrix->cols, sizeof(double));
    time(&time_epi_s);
#pragma omp parallel for schedule(dynamic,50) collapse(2) private(k,i,j) shared(m_1,m_2)
    for (k=0;k<ind_n;k++) {
        for (i =0;i<mrk_n; i++) {
            for (j=0;j<mrk_n;j++) {
                result_matrix->mat_data[i][j] += m_1->mat_data[i][k]*m_2->mat_data[j][k]*r(k)*var_val;
            }
        }
    }
    time(&time_epi_e);
    fprintf(stderr,"Reading matrix use %d seconds.\n",(int)(time_epi_e-time_epi_s));
    free_matrix(m_1);
    free_matrix(m_2);
    MatrixXd mat_tmp(mrk_n,mrk_n);
    for (i =0;i<mrk_n; i++) {
        for (j=0;j<mrk_n;j++){
            mat_tmp(i,j)=result_matrix->mat_data[i][j];
        }
    }
    free_matrix(result_matrix);
    mrk_epi_eff=MatrixXd(mat_tmp.triangularView<Upper>());
    mat_tmp.resize(0,0);
    VectorXd setdiag_zero=VectorXd::Zero(mrk_n);
    mrk_epi_eff.diagonal()=setdiag_zero;
    mrk_epi_eff/=gdiag_mean_sqrt;
    return mrk_epi_eff;
}
void Greml_ce::get_epi_mrk_her_m_std(const parameters &params,MatrixXd &epi_mrk_eff,VectorXd &her_all,int indx){
    double h2,one_u2mrk_tmp_sum_h2;
    MatrixXd u2mrk_tmp;
    MatrixXd h2mrk_tmp;
    vector<double>mean_stddev;
    u2mrk_tmp=epi_mrk_eff.cwiseProduct(epi_mrk_eff);
    h2 = var(indx)/var.sum();
    one_u2mrk_tmp_sum_h2=h2/u2mrk_tmp.sum();
    h2mrk_tmp = u2mrk_tmp*one_u2mrk_tmp_sum_h2;
    mean_stddev=mean_std_epi_snp(h2mrk_tmp);
    her_all(0)=mean_stddev[0];
    her_all(1)=mean_stddev[1];
    her_all(2)=one_u2mrk_tmp_sum_h2;
    u2mrk_tmp.resize(0,0);
    h2mrk_tmp.resize(0,0);
}
void Greml_ce::get_epi_mrk_eff_AA_or_DD(const parameters &p,string eff_type,MatrixXd &epi_mrk_eff,VectorXd &her_all,Pairwise_result_one *epistasis_effects){
    int i,j,mrk_index1,mrk_index2;
    double h_mean,h_stdd,u2mrk_tmp_sum_h2;
    h_mean=her_all(0);
    h_stdd=her_all(1);
    u2mrk_tmp_sum_h2=her_all(2);
    string mrk_eff_file=p.outprefix+"_"+eff_type+"_epi_effect.snpe";
    for(i=0;i<(p.num_pairwise_out-1);i++){
        for(j=1;j<(p.num_pairwise_out-i);j++){
            if((epistasis_effects[j-1].snp1>epistasis_effects[j].snp1)||\
            ((epistasis_effects[j-1].snp1==epistasis_effects[j].snp1)&&\
            ((epistasis_effects[j-1].snp2>epistasis_effects[j].snp2)||\
            (epistasis_effects[j-1].snp2==epistasis_effects[j].snp2)))){
                swap(epistasis_effects[j-1],epistasis_effects[j]);
            }
        }
    }
    double h2_eff_k,h2_eff_k_all,h2_eff_k_norm,h2_eff_k_norm_all,mrk_eff_k;
    ofstream output_mrk(mrk_eff_file.c_str());
    output_mrk<<setw(5)<<right<<"Chr1"<<setw(17)<<right<<"SNP1"<<setw(17)<<right<<"Pos1"<<setw(17)<<right<<"Chr2"<<setw(17)<<right<<"SNP2"<<setw(17)<<right<<"Pos2";
    output_mrk<<setw(17)<<right<<"Effect_"+eff_type<<setw(17)<<right<<"m_effect_"+eff_type\
    <<setw(17)<<right<<"Effect_abs_"+eff_type<<setw(17)<<right<<"m_effect_abs_"+eff_type<<setw(17)<<right<<"h2_mrk_"+eff_type\
    <<setw(17)<<right<<"m_h2_mrk_"+eff_type<<setw(17)<<right<<"h2_mrk_norm_"+eff_type<<setw(17)<<right<<"m_h2_mrk_norm_"+eff_type<<endl;
    for(j=0;j<p.num_pairwise_out;j++){
        mrk_index1=epistasis_effects[j].snp1;
        mrk_index2=epistasis_effects[j].snp2;
        output_mrk<<setw(5)<<right<<mysnpmap[mrk_index1].chr<<setw(17)<<right<<mysnpmap[mrk_index1].mrk<<setw(17)<<right<<mysnpmap[mrk_index1].pos<<setw(17)<<right\
        <<mysnpmap[mrk_index2].chr<<setw(17)<<right<<mysnpmap[mrk_index2].mrk<<setw(17)<<right<<mysnpmap[mrk_index2].pos<<setw(17)<<right;
        mrk_eff_k=epi_mrk_eff(mrk_index1,mrk_index2);
        h2_eff_k=pow(mrk_eff_k,2)*u2mrk_tmp_sum_h2;
        h2_eff_k_all=h2_eff_k*mrk_epi_n;
        h2_eff_k_norm=(h2_eff_k-h_mean)/h_stdd;
        h2_eff_k_norm_all=h2_eff_k_norm*mrk_epi_n;
        output_mrk<<setw(17)<<right<<setprecision(6)<<scientific<<mrk_eff_k<<setw(17)<<right<<(mrk_eff_k*mrk_epi_n)<<setw(17)<<right\
        <<abs(mrk_eff_k)<<setw(17)<<right<<(abs(mrk_eff_k)*mrk_epi_n)<<setw(17)<<right\
        <<h2_eff_k<<setw(17)<<right<<h2_eff_k_all<<setw(17)<<right<<h2_eff_k_norm<<setw(17)<<right<<h2_eff_k_norm_all<<endl;
    }
    output_mrk.close();
}
void Greml_ce::get_epi_mrk_eff_AD(const parameters &p,string eff_type,MatrixXd &epi_mrk_eff1,MatrixXd &epi_mrk_eff2,VectorXd &her_all1,VectorXd &her_all2,Pairwise_result_one *epistasis_effects){
    int i,j,mrk_index1,mrk_index2;
    double h_mean,h_stdd,u2mrk_tmp_sum_h2;
    double h_mean2,h_stdd2,u2mrk_tmp_sum_h22;
    h_mean=her_all1(0);
    h_stdd=her_all1(1);
    u2mrk_tmp_sum_h2=her_all1(2);
    h_mean2=her_all2(0);
    h_stdd2=her_all2(1);
    u2mrk_tmp_sum_h22=her_all2(2);
    string mrk_eff_file=p.outprefix+"_"+eff_type+"_effect.snpe";
    for(i=0;i<(p.num_pairwise_out-1);i++){
        for(j=1;j<(p.num_pairwise_out-i);j++){
            if((epistasis_effects[j-1].snp1>epistasis_effects[j].snp1)||\
            ((epistasis_effects[j-1].snp1==epistasis_effects[j].snp1)&&\
            ((epistasis_effects[j-1].snp2>epistasis_effects[j].snp2)||\
            (epistasis_effects[j-1].snp2==epistasis_effects[j].snp2)))){
                swap(epistasis_effects[j-1],epistasis_effects[j]);
            }
        }
    }
    double h2_eff_k,h2_eff_k_all,h2_eff_k_norm,h2_eff_k_norm_all,mrk_eff_k;
    double h2_eff_k2,h2_eff_k_all2,h2_eff_k_norm2,h2_eff_k_norm_all2,mrk_eff_k2;
    ofstream output_mrk(mrk_eff_file.c_str());
    output_mrk<<setw(5)<<right<<"Chr1"<<setw(17)<<right<<"SNP1"<<setw(17)<<right<<"Pos1"<<setw(17)<<right<<"Chr2"<<setw(17)<<right<<"SNP2"<<setw(17)<<right<<"Pos2";
    output_mrk<<setw(17)<<right<<"Effect_AD"<<setw(17)<<right<<"m_effect_AD"\
    <<setw(17)<<right<<"Effect_abs_AD"<<setw(17)<<right<<"m_effect_abs_AD"<<setw(17)<<right<<"h2_mrk_AD"\
    <<setw(17)<<right<<"m_h2_mrk_AD"<<setw(17)<<right<<"h2_mrk_norm_AD"<<setw(17)<<right<<"m_h2_mrk_norm_AD"\
    <<setw(17)<<right<<"Effect_DA"<<setw(17)<<right<<"m_effect_DA"<<setw(17)<<right<<"Effect_abs_DA"\
    <<setw(17)<<right<<"m_effect_abs_DA"<<setw(17)<<right<<"h2_mrk_DA"<<setw(17)<<right<<"m_h2_mrk_DA"<<setw(17)<<right<<"h2_mrk_norm_DA"\
    <<setw(17)<<right<<"m_h2_mrk_norm_DA"<<endl;
    for(j=0;j<p.num_pairwise_out;j++){
        mrk_index1=epistasis_effects[j].snp1;
        mrk_index2=epistasis_effects[j].snp2;
        output_mrk<<setw(5)<<right<<mysnpmap[mrk_index1].chr<<setw(17)<<right<<mysnpmap[mrk_index1].mrk<<setw(17)<<right<<mysnpmap[mrk_index1].pos<<setw(17)<<right\
        <<mysnpmap[mrk_index2].chr<<setw(17)<<right<<mysnpmap[mrk_index2].mrk<<setw(17)<<right<<mysnpmap[mrk_index2].pos<<setw(17)<<right;
        mrk_eff_k=epi_mrk_eff1(mrk_index1,mrk_index2);
        h2_eff_k=pow(mrk_eff_k,2)*u2mrk_tmp_sum_h2;
        h2_eff_k_all=h2_eff_k*mrk_epi_n;
        h2_eff_k_norm=(h2_eff_k-h_mean)/h_stdd;
        h2_eff_k_norm_all=h2_eff_k_norm*mrk_epi_n;
        output_mrk<<setw(17)<<right<<setprecision(6)<<scientific<<mrk_eff_k<<setw(17)<<right<<(mrk_eff_k*mrk_epi_n)<<setw(17)<<right\
        <<abs(mrk_eff_k)<<setw(17)<<right<<(abs(mrk_eff_k)*mrk_epi_n)<<setw(17)<<right\
        <<h2_eff_k<<setw(17)<<right<<h2_eff_k_all<<setw(17)<<right<<h2_eff_k_norm<<setw(17)<<right<<h2_eff_k_norm_all;
        mrk_eff_k2=epi_mrk_eff2(mrk_index1,mrk_index2);
        h2_eff_k2=pow(mrk_eff_k2,2)*u2mrk_tmp_sum_h22;
        h2_eff_k_all2=h2_eff_k2*mrk_epi_n;
        h2_eff_k_norm2=(h2_eff_k2-h_mean2)/h_stdd2;
        h2_eff_k_norm_all2=h2_eff_k_norm2*mrk_epi_n;
        output_mrk<<setw(17)<<right<<setprecision(6)<<scientific<<mrk_eff_k2<<setw(17)<<right<<(mrk_eff_k2*mrk_epi_n)<<setw(17)<<right\
        <<abs(mrk_eff_k2)<<setw(17)<<right<<(abs(mrk_eff_k2)*mrk_epi_n)<<setw(17)<<right\
        <<h2_eff_k2<<setw(17)<<right<<h2_eff_k_all2<<setw(17)<<right<<h2_eff_k_norm2<<setw(17)<<right<<h2_eff_k_norm_all2<<endl;
    }
    output_mrk.close();
}
void Greml_ce::get_single_mrk_eff(const parameters &p,vector<VectorXd> &single_mrk_eff_all){
    int i,j;
    double h2,tmp,tmpH,H_mean,H_stdd,one_u2mrk_tmp_sum_h2,onesite_val;
    VectorXd H2mrk = VectorXd::Zero(mrk_n);
    VectorXd u2mrk_tmp;
    VectorXd h2mrk_tmp;
    VectorXd h_mean;
    VectorXd h_stdd;
    VectorXd u2mrk_tmp_sum_h2;
    h_mean=VectorXd::Zero(var_name.size()-1);
    h_stdd=VectorXd::Zero(var_name.size()-1);
    u2mrk_tmp_sum_h2=VectorXd::Zero(var_name.size()-1);
	vector<double>mean_stddev;
    string mrk_eff_file=p.outprefix+"_snp_effect.snpe";
    for (i=0;i<var_name.size()-1;i++){
        if((var_name[i]=="A")||(var_name[i]=="D")){
			u2mrk_tmp = single_mrk_eff_all[i].cwiseProduct(single_mrk_eff_all[i]);
			h2 = var(i)/var.sum();
			one_u2mrk_tmp_sum_h2=h2/u2mrk_tmp.sum();
			h2mrk_tmp = u2mrk_tmp*one_u2mrk_tmp_sum_h2;
			mean_stddev=mean_std_single_snp(h2mrk_tmp);
            h_mean(i)=mean_stddev[0];
            h_stdd(i)=mean_stddev[1];
			H2mrk += h2mrk_tmp;
			u2mrk_tmp_sum_h2(i)=one_u2mrk_tmp_sum_h2;
        }
    }
    u2mrk_tmp.resize(0);
    h2mrk_tmp.resize(0);
    mean_stddev=mean_std_single_snp(H2mrk);
    H_mean=mean_stddev[0];
    H_stdd=mean_stddev[1];
    ofstream output_mrk(mrk_eff_file.c_str());
    output_mrk<<setw(17)<<right<<"Chr"<<setw(17)<<right<<"SNP"<<setw(17)<<right<<"Pos";
    for(i=0;i<var_name.size()-1;i++){
        if((var_name[i]=="A")||(var_name[i]=="D")){
            output_mrk << setw(17) << right << "Effect_" + var_name[i];
            output_mrk << setw(17) << right << "m_effect_" + var_name[i];
        }
    }
    for(i=0;i<var_name.size()-1;i++){
        if((var_name[i]=="A")||(var_name[i]=="D")){
            output_mrk << setw(17) << right << "Effect_" + var_name[i] + "2";
	        output_mrk << setw(17) << right << "m_effect_" + var_name[i] + "2";
	    }
    }
	for(i=0;i<var_name.size()-1;i++){
	    if((var_name[i]=="A")||(var_name[i]=="D")){
	        output_mrk << setw(17) << right << "h2_mrk_" + var_name[i];
	        output_mrk << setw(17) << right << "m_h2_mrk_" + var_name[i];
	    }
	}
	output_mrk << setw(17) << right << "H2_mrk";
	for(i=0;i<var_name.size()-1;i++){
	    if((var_name[i]=="A")||(var_name[i]=="D")){
	        output_mrk << setw(17) << right << "h2_mrk_norm_" + var_name[i];
	        output_mrk << setw(17) << right << "m_h2_mrk_norm_" + var_name[i];
	    }
	}
	output_mrk<<setw(17)<<right<<"H2_mrk_norm"<<setw(17)<<right<<"m_H2_mrk_norm";
	output_mrk<<endl;
	for(j=0;j<mrk_n;j++){
        output_mrk<<setw(17)<<right<<mysnpmap[j].chr<<setw(17)<<right<<mysnpmap[j].mrk<<setw(17)<<right<<setprecision(0)<<mysnpmap[j].pos;
        for(int jj=0;jj<single_mrk_eff_all.size();jj++){
            output_mrk<<setw(17)<<right<<setprecision(6)<<scientific<<(single_mrk_eff_all[jj])(j);
            output_mrk<<setw(17)<<right<<setprecision(6)<<scientific<<((single_mrk_eff_all[jj])(j)*mrk_n);
        }
        for(int jj=0;jj<single_mrk_eff_all.size();jj++){
	        output_mrk<<setw(17)<<right<<setprecision(6)<<scientific<<abs((single_mrk_eff_all[jj])(j));
			output_mrk<<setw(17)<<right<<setprecision(6)<<scientific<<(abs((single_mrk_eff_all[jj])(j))*mrk_n);
		}
		for (i=0;i<var_name.size()-1;i++){
            if((var_name[i]=="A")||(var_name[i]=="D")){
                onesite_val=pow(single_mrk_eff_all[i](j),2)*u2mrk_tmp_sum_h2(i);
                output_mrk<<setw(17)<<right<<setprecision(6)<<scientific<<onesite_val;
                output_mrk<<setw(17)<<right<<setprecision(6)<<scientific<<(onesite_val*mrk_n);
            }
        }
        output_mrk<<setw(17)<<right<<setprecision(6)<<scientific<<H2mrk(j);
        for (i=0;i<var_name.size()-1;i++){
            if((var_name[i]=="A")||(var_name[i]=="D")){
                onesite_val=pow(single_mrk_eff_all[i](j),2)*u2mrk_tmp_sum_h2(i);
                tmp = (onesite_val-h_mean(i))/h_stdd(i);
                output_mrk<<setw(17)<<right<<setprecision(6)<<scientific<<tmp;
                output_mrk<<setw(17)<<right<<setprecision(6)<<scientific<<(tmp*mrk_n);
            }
        }
        tmpH=(H2mrk(j)-H_mean)/H_stdd;
        output_mrk<<setw(17)<<right<<setprecision(6)<<scientific<<tmpH\
        <<setw(17)<<right<<setprecision(6)<<scientific<<(tmpH*mrk_n);
        output_mrk<<endl;
	}
	output_mrk.close();
}
void Greml_ce::get_mrk_eff_hap(const parameters &p,VectorXd &r){
    time_t time_a;
    time(&time_a);
    int i,j,ii,hapflag,n_2_hapcol,n_hapcol,b,start_b,count1,b1;
    double hetero_total  = 0.0,hetero_total2 = 0.0,alpha = 1.0,beta  = 0.0;
    string hap_input_file=p.haplotype,hap_eff_file=p.outprefix+"_haplotype_effect.snpe";
    if((p.genotype=="null")&&(p.haplotype!="null")){
        ind_n=0;
        int linetmp_num=0;
        string   haplinetmp;
        ifstream input_hap(hap_input_file);
        while(getline(input_hap,haplinetmp).good()){
            linetmp_num++;
        }
        input_hap.close();
        haplinetmp.clear();
        ind_n=linetmp_num-1;
    }
    ifstream input_hap(hap_input_file);
    if(!input_hap){
	    fprintf(stderr,"ERROR: Can not open the file %s\n",hap_input_file.c_str());
	    throw;
    }
    else{
	    fprintf(stderr,"Reading haplotype file %s ...\n",hap_input_file.c_str());
    }
    string   hapline;
    getline(input_hap,hapline);
    vector<string> hapwords = str_split(" \r\t\n\x0b", hapline);
    hapwords.erase(hapwords.begin());
    n_2_hapcol=hapwords.size();
    n_hapcol=n_2_hapcol/2;
    hapwords.clear();
    hapline.clear();
    MatrixXi geno_hap  = MatrixXi::Zero(ind_n,n_2_hapcol);
    hb_allele.resize(n_hapcol);
    hb_het_geno.resize(n_hapcol);
    hapflag=0;
    while(getline(input_hap,hapline).good()){
        vector<string> hapwords = str_split(" \r\t\n\x0b", hapline);
        if(ind_geno==NULL){
            get_ind_list_geno(hapwords[0],hapflag,&ind_geno);
        }
        else if(HASH_COUNT(ind_geno)!=ind_n){
            get_ind_list_geno(hapwords[0],hapflag,&ind_geno);
        }
        else{
            char name[STR_LEN];
            strcpy(name,hapwords[0].c_str());
            IDX *ind_geno_tmp;
            HASH_FIND_STR(ind_geno,name,ind_geno_tmp);
            if(ind_geno_tmp==NULL){
                fprintf(stderr,"ERROR: Individual is not in haplotype file\n");
                throw;
            }
            else{
                if(ind_geno_tmp->index!=hapflag){
                    fprintf(stderr,"ERROR: The postion of Individual\n");
                    throw;
                }
            }
        }
        for(i=0;i<n_hapcol;i++){
            geno_hap(hapflag,i*2)    = atoi(hapwords[i*2+1].c_str());
            geno_hap(hapflag,(i*2+1))= atoi(hapwords[i*2+2].c_str());
        }
        hapline.clear();
        hapflag++;
    }
    input_hap.close();
 /////////////////////////////////////////////////////////////////////////////////
    int geno_tmp_0,geno_tmp_1,max_tmp;
    double sum;
    for(i=0;i<n_hapcol;i++){
        for(j=0;j<ind_n;j++){
            geno_tmp_0 = geno_hap(j,i*2);
            geno_tmp_1 = geno_hap(j,i*2+1);
            allele_list::const_iterator got;
            if(geno_tmp_0!=0){
                got = hb_allele[i].find(geno_tmp_0);
                if(got==hb_allele[i].end()){
                    hb_allele[i][geno_tmp_0].count=1.0;
                }
                else{
                    hb_allele[i][geno_tmp_0].count+=1.0;
                }
            }
            if(geno_tmp_1!=0){
                got = hb_allele[i].find(geno_tmp_1);
                if(got==hb_allele[i].end()){
                    hb_allele[i][geno_tmp_1].count=1.0;
                }
                else{
                    hb_allele[i][geno_tmp_1].count+=1.0;
                }
            }
            array<int,2> het_geno;
            het_geno_list::const_iterator got_het_geno;
            if(geno_tmp_0!=0 && geno_tmp_1!=0 && geno_tmp_0!=geno_tmp_1){
                if(geno_tmp_0<geno_tmp_1){
                    het_geno[0] = geno_tmp_0;
                    het_geno[1] = geno_tmp_1;
                }
                else if(geno_tmp_0>geno_tmp_1){
                    het_geno[0] = geno_tmp_1;
                    het_geno[1] = geno_tmp_0;
                }
                got_het_geno = hb_het_geno[i].find(het_geno);
                if(got_het_geno==hb_het_geno[i].end()){
                    hb_het_geno[i][het_geno].count=1.0;
                }
                else{
                    hb_het_geno[i][het_geno].count+=1.0;
                }
            }
        }
    }
    hap_max_key.resize(n_hapcol);
/////////////////////////////////////////////////////////////////////////////////
    for(i=0;i<n_hapcol;i++){
        hb_allele[i].erase(0);
        sum = 0.0;
        max_tmp = hb_allele[i].begin()->first;
        for(auto it = hb_allele[i].begin();it!=hb_allele[i].end();++it){
            if((hb_allele[i][max_tmp]).count<(it->second).count){
                max_tmp = it->first;
            }
            sum += (it->second).count;
        }
        hap_max_key(i) = max_tmp;
        count1 = 0;
        for(auto it = hb_allele[i].begin();it!=hb_allele[i].end();++it){
            (it->second).count /= sum;
            if(it->first!=hap_max_key(i)){
                (it->second).idx = count1;
                count1++;
            }
            else{
                (it->second).idx = -1;
            }
        }
        count1 = 0;
        for(auto it = hb_het_geno[i].begin();it!=hb_het_geno[i].end();++it){
            (it->second).idx = count1;
            count1++;
        }
    }
////////////////////////////////////////////////////////////////////////////////////
    cout<<endl<<"Calculating haplotype heritabilities"<<endl;
    ofstream output_mrk(hap_eff_file.c_str());
    output_mrk<<setw(5)<<right<<"HAPID"<<setw(15)<<right<<"h2_hap_AH"<<setw(15)<<right<<"h2_hap_std_AH"<<endl;
    int indx_hap=var_name.size()-2;
	VectorXd u2;
	u2.resize(n_hapcol);
    MatrixXd W_hap;
    int hap_count,hap_max_key_j,vdegree;
    VectorXd u;
    double hap_val,u2_sum,h2,hap_mean,hap_stdd,one_block_value;
    h2=var(indx_hap)/var.sum();
    for(b1=0;b1<n_hapcol;b1++){
        hap_count = hb_allele[b1].size()-1;
        W_hap.resize(ind_n,hap_count);
        hap_max_key_j=hap_max_key(b1);
        each_W_hap_generate(b1,geno_hap,W_hap,hb_allele,ind_n,hap_max_key_j);
        u = W_hap.transpose()*r;
        hap_val=u.dot(u);
        u2(b1)=hap_val;
    }
    W_hap.resize(0,0);
    u.resize(0);
    u2_sum=u2.sum();
    VectorXd hap_ah_her;
    hap_ah_her=u2*(h2/u2_sum);
    u2.resize(0);
    vdegree=n_hapcol-1;
    hap_mean=hap_ah_her.mean();
    hap_stdd=sqrt((hap_ah_her.array()-hap_mean).square().sum()/vdegree);
    for(i=0;i<n_hapcol;i++){
        one_block_value=(hap_ah_her(i)-hap_mean)/hap_stdd;
        output_mrk<<setw(5)<<right<<i<<setw(15)<<right<<setprecision(6)<<scientific<<hap_ah_her(i)<<setw(15)<<right<<one_block_value<<endl;
    }
    output_mrk.close();
    time_t time_b;
    time(&time_b);
    fprintf(stderr,"\nhaplotype_heritibility processing took: %d seconds!\n",(unsigned int)(time_b-time_a));
}
vector<double> Greml_ce::mean_std_single_snp(VectorXd &vec){
    vector<double>out;
    double hmean;
    double stddev;
    int vdegree;
    vdegree=vec.size()-1;
    hmean=vec.mean();
    stddev=sqrt((vec.array()-hmean).square().sum()/vdegree);
    out.push_back(hmean);
    out.push_back(stddev);
    return out;
}
vector<double> Greml_ce::mean_std_epi_snp(MatrixXd &vec){
    vector<double>out;
    double hmean;
    double stddev;
    int vdegree;
    vdegree=mrk_epi_n-1;
    hmean=vec.sum()/mrk_epi_n;
    MatrixXd vec1;
    vec1=vec.array().pow(2).matrix()-2*hmean*vec;
    stddev=sqrt((vec1.sum()+mrk_epi_n*pow(hmean,2))/vdegree);
    out.push_back(hmean);
    out.push_back(stddev);
    return out;
}
