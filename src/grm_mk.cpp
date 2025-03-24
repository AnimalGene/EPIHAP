#include "read_geno.h"
void Greml_ce::init(const parameters &p){
    ind_idx  = NULL;
    ind_geno = NULL;
    round_geno = 0;
    var_name_values["A"]=p.additive_variances;
    var_name_values["D"]=p.dominance_variances;
    var_name_values["AA"]=p.additive_additive_variances;
    var_name_values["AA-inter"]=p.additive_additive_inter_variances;
    var_name_values["AA-intra"]=p.additive_additive_intra_variances;
    var_name_values["AD"]=p.additive_dominance_variances;
    var_name_values["AD-inter"]=p.additive_dominance_inter_variances;
    var_name_values["AD-intra"]=p.additive_dominance_intra_variances;
    var_name_values["DD"]=p.dominance_dominance_variances;
    var_name_values["DD-inter"]=p.dominance_dominance_inter_variances;
    var_name_values["DD-intra"]=p.dominance_dominance_intra_variances;
    var_name_values["AAA"]=p.additive_additive_additive_variances;
    var_name_values["AAD"]=p.additive_additive_dominance_variances;
    var_name_values["ADD"]=p.additive_dominance_dominance_variances;
    var_name_values["DDD"]=p.dominance_dominance_dominance_variances;
    var_name_values["AH"]=p.additive_haplotype_variances;
    var_name_values["E"] =p.residual_variances;
    var_name_flags["A"]=p.additive_variances_f;
    var_name_flags["D"]=p.dominance_variances_f;
    var_name_flags["AA"]=p.additive_additive_variances_f;
    var_name_flags["AA-inter"]=p.additive_additive_inter_variances_f;
    var_name_flags["AA-intra"]=p.additive_additive_intra_variances_f;
    var_name_flags["AD"]=p.additive_dominance_variances_f;
    var_name_flags["AD-inter"]=p.additive_dominance_inter_variances_f;
    var_name_flags["AD-intra"]=p.additive_dominance_intra_variances_f;
    var_name_flags["DD"]=p.dominance_dominance_variances_f;
    var_name_flags["DD-inter"]=p.dominance_dominance_inter_variances_f;
    var_name_flags["DD-intra"]=p.dominance_dominance_intra_variances_f;
    var_name_flags["AAA"]=p.additive_additive_additive_variances_f;
    var_name_flags["AAD"]=p.additive_additive_dominance_variances_f;
    var_name_flags["ADD"]=p.additive_dominance_dominance_variances_f;
    var_name_flags["DDD"]=p.dominance_dominance_dominance_variances_f;
    var_name_flags["AH"]=p.additive_haplotype_variances_f;
    var_name_flags["E"] =p.residual_variances_f;
    if(p.chrsplit_f==true){
        G_diag_mean.resize(10);
        var_name={"A","D","AA-inter","AA-intra","AD-inter","AD-intra","DD-inter","DD-intra","AH","E"};
    }
    else{
        G_diag_mean.resize(11);
        var_name={"A","D","AA","AD","DD","AAA","AAD","ADD","DDD","AH","E"};
    }
    cerr<<"A total of "<<omp_get_num_procs()<<" threads are available.\n\n";
    omp_set_num_threads(p.numThreads);
    Eigen::initParallel();
    Eigen::setNbThreads(p.numThreads);
    cerr<<Eigen::nbThreads()<<" threads are used.\n";
    if(p.factors_pos_f==true){
        vector<string> fixed_tmp=str_split(" ",p.fixed_factor_column);
        for (int j=1;j<fixed_tmp.size();j++) cross.push_back(atoi(fixed_tmp[j].c_str()));
        fixed_tmp.clear();
    }
    if(p.covar_pos_f==true){
        vector<string> covariate_tmp=str_split(" ",p.covariate_factor_column);
        for (int j=1;j<covariate_tmp.size();j++) covar.push_back(atoi(covariate_tmp[j].c_str()));
        covariate_tmp.clear();
    }
    string ma_aa,ma_ad,ma_da,ma_dd;
    ma_aa="AA",ma_ad="AD",ma_da="DA",ma_dd="DD";
    my_ma_type[ma_aa]=0;
    my_ma_type[ma_ad]=1;
    my_ma_type[ma_da]=2;
    my_ma_type[ma_dd]=3;
    for(int i=0;i<4;i++){
        MatrixXi epi_indx(3,3);
        epi_indx(0,0)=0+i*9;
        epi_indx(0,1)=1+i*9;
        epi_indx(0,2)=2+i*9;
        epi_indx(1,0)=3+i*9;
        epi_indx(1,1)=4+i*9;
        epi_indx(1,2)=5+i*9;
        epi_indx(2,0)=6+i*9;
        epi_indx(2,1)=7+i*9;
        epi_indx(2,2)=8+i*9;
        my_geno_stat_epi_idx.push_back(epi_indx);
        epi_indx.resize(0,0);
    }
}
void Greml_ce::read_genotype_method_1(const parameters &p,ofstream &flog){
//    time_t time_a,time_b,time_a_geno,time_b_geno,time_a_w;
    clock_t time_a,time_b,time_a_geno,time_b_geno,time_a_w;
//    time(&time_a);
    time_a = clock();
    int i,b,start_b,is_hwd=0;
    double hetero_total  = 0.0,hetero_total2 = 0.0,alpha = 1.0,beta  = 1.0;
    string SA_filename =p.grm_file_prefix+".g.A",SD_filename =p.grm_file_prefix+".g.D",\
    SAA_filename =p.grm_file_prefix+".g.AA",SAD_filename =p.grm_file_prefix+".g.AD",SDD_filename =p.grm_file_prefix+".g.DD",\
    SAA_intra_filename =p.grm_file_prefix+".g.AA-intra",SAD_intra_filename =p.grm_file_prefix+".g.AD-intra",SDD_intra_filename =p.grm_file_prefix+".g.DD-intra",\
    SAA_inter_filename =p.grm_file_prefix+".g.AA-inter",SAD_inter_filename =p.grm_file_prefix+".g.AD-inter",SDD_inter_filename =p.grm_file_prefix+".g.DD-inter",\
    SAAA_filename =p.grm_file_prefix+".g.AAA",SAAD_filename =p.grm_file_prefix+".g.AAD",SADD_filename =p.grm_file_prefix+".g.ADD",\
    SDDD_filename =p.grm_file_prefix+".g.DDD";
    string g_diag_filename=p.grm_file_prefix+".gdiag",indgeno_filename=p.grm_file_prefix+".indgeno";
//    time(&time_a_geno);
    time_a_geno = clock();
    MatrixXi geno;
    MatrixXd geno_stat;
    read_all_genofiles(p,geno,geno_stat,flog);
//    if(p.hap_f==true) flog<<"Haplotype genotypes file                           :\t"<<p.haplotype<<endl;
//    flog<<"The prefix of GRM files                            :\t"<<p.grm_file_prefix<<endl;
//    flog<<"The number of threads is                           :\t"<<p.numThreads<<endl<<endl;
    flog<<"The number of SNPs is                              :\t"<<mrk_n<<endl;
//    time(&time_b_geno);
    time_b_geno = clock();
//    cerr<<"\nThe time (second) cost for reading genotypes is:\t"<<(int)(time_b_geno-time_a_geno)<<endl;
    flog<<"The time (seconds) cost for reading genotypes is   :\t"<<(int)((double)(time_b_geno-time_a_geno)/CLOCKS_PER_SEC)<<endl;
    genostat_update(geno_stat);
    get_het_sum(geno_stat,is_hwd,&hetero_total,&hetero_total2);
    int block_size=20000;
    save_wt(p,geno,geno_stat,block_size);
    VectorXd sgdiag;
    sgdiag=(p.chrsplit_f)?VectorXd::Zero(8):VectorXd::Zero(9);
    int remainders = mrk_n % block_size;
    vector<MatrixXd>WWT_AD;
    MatrixXd WWT_tmp,W;
    for(i=0;i<var_name.size();i++){
        if((var_name[i]=="A")||(var_name[i]=="D")){
            WWT_tmp.noalias() = MatrixXd::Zero(ind_n,ind_n);
            if(mrk_n>block_size){
                for (b=block_size;b<=mrk_n;b+=block_size){
                    W = MatrixXd::Zero(ind_n,block_size);
                    geno_transform(geno.block(0,(b-block_size),ind_n,block_size),geno_stat.block(0,(b-block_size),12,block_size),var_name[i],is_hwd,W);
                    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,ind_n,ind_n,block_size,alpha,W.data(),ind_n,W.data(),ind_n,beta,WWT_tmp.data(),ind_n);
                }
                if(remainders>0){
                    start_b=mrk_n-remainders;
                    W = MatrixXd::Zero(ind_n,remainders);
                    geno_transform(geno.block(0,start_b,ind_n,remainders),geno_stat.block(0,start_b,12,remainders),var_name[i],is_hwd,W);
                    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,ind_n,ind_n,remainders,alpha,W.data(),ind_n,W.data(),ind_n,beta,WWT_tmp.data(),ind_n);
                }
                else{
                    ;
                    }
            }
            else{
                W = MatrixXd::Zero(ind_n,mrk_n);
                geno_transform(geno,geno_stat,var_name[i],is_hwd,W);
                cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,ind_n,ind_n,mrk_n,alpha,W.data(),ind_n,W.data(),ind_n,beta,WWT_tmp.data(),ind_n);
                }
            WWT_AD.push_back(WWT_tmp);
         }
    }
    WWT_tmp.resize(0,0);
    MatrixXd WWTA,WWTD,WWTAA,WWTAD,WWTDD,SAA,SAD,SDD,SAAA,SAAD,SADD,SDDD,SA,SD;
    double ka,kd,kaa,kad,kdd;
    for(i=0;i<var_name.size();i++){
        if(var_name[i]=="A"){
            WWTA=WWT_AD[i];
            ka=WWTA.diagonal().mean();
            SA=WWTA/ka;
            write_binary(SA_filename.c_str(),SA);
            SA_filename.clear();
            sgdiag(0)=sqrt(ka);
            cerr<<"\n\n*******************SA************\n\n";
            cerr<<SA.block(0,0,5,5)<<endl;

        }
        else if (var_name[i]=="D"){
            WWTD=WWT_AD[i];
            kd=WWTD.diagonal().mean();
            SD=WWTD/kd;
            write_binary(SD_filename.c_str(),SD);
            SD_filename.clear();
            sgdiag(1)=sqrt(kd);
            cerr<<"\n\n*******************SD************\n\n";
            cerr<<SD.block(0,0,5,5)<<endl;
        }
        else{
            ;
        }
    }
    for(i=0;i<var_name.size();i++){
        if((var_name[i]=="A")||(var_name[i]=="D")) WWT_AD[i].resize(0,0);
    }
    MatrixXd SAA_intra,SAD_intra,SDD_intra,SAA_inter,SAD_inter,SDD_inter,WWTAA_intra,WWTAD_intra,WWTDD_intra,WWTAA_inter,WWTAD_inter,WWTDD_inter;
    string varname;
    double kaa_intra,kad_intra,kdd_intra,kaa_inter,kad_inter,kdd_inter;
    WWTAA=WWTA.cwiseProduct(WWTA);
    if(p.chrsplit_f){
        SA.resize(0,0);
        SD.resize(0,0);
        varname="A";
        henderson_chrsplit_WWTAA_WWTDD(WWTAA_intra,geno,geno_stat,varname,is_hwd);
        kaa_intra=WWTAA_intra.diagonal().mean();
        WWTAA_inter=WWTAA - WWTAA_intra;
        kaa_inter=WWTAA_inter.diagonal().mean();
        sgdiag(2)=sqrt(kaa_intra);
        sgdiag(3)=sqrt(kaa_inter);
        cerr<<"\n\n*******************SAA_intra************\n"<<endl;
        WWTAA_intra/=kaa_intra;
        write_binary(SAA_intra_filename.c_str(),WWTAA_intra);
        SAA_intra_filename.clear();
        cerr<<WWTAA_intra.block(0,0,5,5)<<endl;
        WWTAA_intra.resize(0,0);
        cerr<<"\n\n*******************SAA_inter************\n"<<endl;
        WWTAA_inter/=kaa_inter;
        write_binary(SAA_inter_filename.c_str(),WWTAA_inter);
        SAA_inter_filename.clear();
        cerr<<WWTAA_inter.block(0,0,5,5)<<endl;
        WWTAA_inter.resize(0,0);
    }
    else{
        geno.resize(0,0);
        geno_stat.resize(0,0);
        kaa=WWTAA.diagonal().mean();
        sgdiag(2)=sqrt(kaa);
        cerr<<"\n\n*******************SAA************\n"<<endl;
        WWTAA/=kaa;
        write_binary(SAA_filename.c_str(),WWTAA);
        SAA_filename.clear();
        cerr<<WWTAA.block(0,0,5,5)<<endl;

        SAAA=WWTAA.cwiseProduct(SA);
        write_binary(SAAA_filename.c_str(),SAAA);
        SAAA_filename.clear();
        sgdiag(5)=sqrt(ka*ka*ka);
        cerr<<"\n\n*******************SAAA************\n"<<endl;
        cerr<<SAAA.block(0,0,5,5)<<endl;
        SAAA.resize(0,0);

        SAAD=WWTAA.cwiseProduct(SD);
        write_binary(SAAD_filename.c_str(),SAAD);
        SAAD_filename.clear();
        sgdiag(6)=sqrt(ka*ka*kd);
        cerr<<"\n\n*******************SAAD************\n"<<endl;
        cerr<<SAAD.block(0,0,5,5)<<endl;
        SAAD.resize(0,0);
    }
    WWTAA.resize(0,0);
    WWTAD=WWTA.cwiseProduct(WWTD);
    WWTA.resize(0,0);
    if(p.chrsplit_f){
        henderson_chrsplit_WWTAD(WWTAD_intra,geno,geno_stat,is_hwd);
        kad_intra=WWTAD_intra.diagonal().mean();
        WWTAD_inter=WWTAD-WWTAD_intra;
        kad_inter=WWTAD_inter.diagonal().mean();
        sgdiag(4)=sqrt(kad_intra);
        sgdiag(5)=sqrt(kad_inter);
        cerr<<"\n\n*******************SAD_intra************\n"<<endl;
        WWTAD_intra/=kad_intra;
        write_binary(SAD_intra_filename.c_str(),WWTAD_intra);
        SAD_intra_filename.clear();
        cerr<<WWTAD_intra.block(0,0,5,5)<<endl;
        WWTAD_intra.resize(0,0);
        cerr<<"\n\n*******************SAD_inter************\n"<<endl;
        WWTAD_inter/=kad_inter;
        write_binary(SAD_inter_filename.c_str(),WWTAD_inter);
        SAD_inter_filename.clear();
        cerr<<WWTAD_inter.block(0,0,5,5)<<endl;
        WWTAD_inter.resize(0,0);
    }
    else{
        kad=WWTAD.diagonal().mean();
        sgdiag(3)=sqrt(kad);
        cerr<<"\n\n*******************SAD************\n"<<endl;
        WWTAD/=kad;
        write_binary(SAD_filename.c_str(),WWTAD);
        SAD_filename.clear();
        cerr<<WWTAD.block(0,0,5,5)<<endl;
    }
    WWTAD.resize(0,0);
    WWTDD=WWTD.cwiseProduct(WWTD);
    WWTD.resize(0,0);
    if(p.chrsplit_f){
        varname="D";
        henderson_chrsplit_WWTAA_WWTDD(WWTDD_intra,geno,geno_stat,varname,is_hwd);
        geno.resize(0,0);
        geno_stat.resize(0,0);
        kdd_intra=WWTDD_intra.diagonal().mean();
        WWTDD_inter=WWTDD-WWTDD_intra;
        kdd_inter=WWTDD_inter.diagonal().mean();
        sgdiag(6)=sqrt(kdd_intra);
        sgdiag(7)=sqrt(kdd_inter);
        cerr<<"\n\n*******************SDD_intra************\n"<<endl;
        WWTDD_intra/=kdd_intra;
        write_binary(SDD_intra_filename.c_str(),WWTDD_intra);
        SDD_intra_filename.clear();
        cerr<<WWTDD_intra.block(0,0,5,5)<<endl;
        WWTDD_intra.resize(0,0);
        cerr<<"\n\n*******************SDD_inter************\n"<<endl;
        WWTDD_inter/=kdd_inter;
        write_binary(SDD_inter_filename.c_str(),WWTDD_inter);
        SDD_inter_filename.clear();
        cerr<<WWTDD_inter.block(0,0,5,5)<<endl;
        WWTDD_inter.resize(0,0);
    }
    else{
        kdd=WWTDD.diagonal().mean();
        sgdiag(5)=sqrt(kdd);
        cerr<<"\n\n*******************SDD************\n"<<endl;
        WWTDD/=kdd;
        write_binary(SDD_filename.c_str(),WWTDD);
        SDD_filename.clear();
        cerr<<WWTDD.block(0,0,5,5)<<endl;

        SADD=SA.cwiseProduct(WWTDD);
        write_binary(SADD_filename.c_str(),SADD);
        SADD_filename.clear();
        sgdiag(7)=sqrt(ka*kd*kd);
        cerr<<"\n\n*******************SADD************\n"<<endl;
        cerr<<SADD.block(0,0,5,5)<<endl;
        SADD.resize(0,0);
        SA.resize(0,0);

        SDDD=WWTDD.cwiseProduct(SD);
        SD.resize(0,0);
        write_binary(SDDD_filename.c_str(),SDDD);
        SDDD_filename.clear();
        sgdiag(8)=sqrt(kd*kd*kd);
        cerr<<"\n\n*******************SDDD************\n"<<endl;
        cerr<<SDDD.block(0,0,5,5)<<endl;
        WWTDD.resize(0,0);
        SDDD.resize(0,0);
    }
    write_binary(g_diag_filename.c_str(),sgdiag);
    flog<<"The method used for epistatic GRMs inference is    :\tAGRM method"<<endl;
    for (i=0;i<sgdiag.size();i++) {
        flog<<setprecision(6)<<scientific<<"The mean of the diagonal elements for K"<<var_name[i];
        if(var_name[i].size()==1){
            flog<<" is        ";
        }
        else if (var_name[i].size()==2){
            flog<<" is       ";
        }
        else if (var_name[i].size()==3){
            flog<<" is      ";
        }
        else{
            flog<<" is ";
        }
        flog<<":\t"<<pow(sgdiag[i],2)<<endl;
    }
    sgdiag.resize(0);
//    time(&time_b);
    time_b = clock();
    flog<<"The time (seconds) cost for SNP GRMs inference is  :\t"<<(int)((double)(time_b-time_a)/CLOCKS_PER_SEC)<<endl;
    ifstream input_indgeno(indgeno_filename);
    if(!input_indgeno) save_ind_geno(indgeno_filename,&ind_geno);
    input_indgeno.close();
}
void Greml_ce::save_wt(const parameters &p,const MatrixXi &geno_org,const MatrixXd &geno_stat,int block_size){
    string wat_filename=p.grm_file_prefix+".g.WAT",wdt_filename=p.grm_file_prefix+".g.WDT";
    string wat_filename1=p.grm_file_prefix+".g.WAT.txt",wdt_filename1=p.grm_file_prefix+".g.WDT.txt";
    int i,b,start_b,remainders,is_hwd=0;
    remainders = mrk_n % block_size;
    MatrixXd Wsave;
    MatrixXd WsaveT;
    cerr<<"\n\n print WA' and WD'!\n\n";
    for (i=0;i<var_name.size()-1;i++){
        if((var_name[i]=="A")||(var_name[i]=="D")){
            WsaveT=MatrixXd::Zero(mrk_n,ind_n);
            if(mrk_n>block_size){
                for (b=block_size;b<=mrk_n;b+=block_size){
                    Wsave = MatrixXd::Zero(ind_n,block_size);
                    geno_transform(geno_org.block(0,(b-block_size),ind_n,block_size),geno_stat.block(0,(b-block_size),12,block_size),var_name[i],is_hwd,Wsave);
                    WsaveT.block((b-block_size),0,block_size,ind_n)=Wsave.transpose();
                }
                if(remainders>0){
                    start_b=mrk_n-remainders;
                    Wsave = MatrixXd::Zero(ind_n,remainders);
                    geno_transform(geno_org.block(0,start_b,ind_n,remainders),geno_stat.block(0,start_b,12,remainders),var_name[i],is_hwd,Wsave);
                    WsaveT.block(start_b,0,remainders,ind_n)=Wsave.transpose();
                }
                else{
                    ;
                }
            }
            else{
                Wsave = MatrixXd::Zero(ind_n,mrk_n);
                geno_transform(geno_org,geno_stat,var_name[i],is_hwd,Wsave);
                WsaveT=Wsave.transpose();
            }
            if(var_name[i]=="A"){
//                write_binary(wat_filename.c_str(),WsaveT);
                cerr<<"\n\n*****************WA'*************\n\n";
                cerr<<WsaveT.block(0,0,5,5)<<endl;
                saveData(wat_filename1.c_str(),WsaveT);
                WsaveT.resize(0,0);

            }
            else if (var_name[i]=="D"){
//                write_binary(wdt_filename.c_str(),WsaveT);
                cerr<<"\n\n*****************WD'*************\n\n";
                cerr<<WsaveT.block(0,0,5,5)<<endl;
                saveData(wdt_filename1.c_str(),WsaveT);
                WsaveT.resize(0,0);
            }
            else{
                cerr<<"\nerror in WsaveT\n"<<endl;
            }
        }
    }
    Wsave.resize(0,0);
}
void Greml_ce::read_genotype_method_2(const parameters &p,ofstream &flog){
    clock_t time_a,time_b,time_a_geno,time_b_geno,time_a_w;
    time_a = clock();
    int i,j,is_hwd=0,order;
    double hetero_total  = 0.0,hetero_total2 = 0.0,alpha = 1.0,beta  = 0.0;
    string SA_filename =p.grm_file_prefix+".g.A",SD_filename =p.grm_file_prefix+".g.D",\
    SAA_filename =p.grm_file_prefix+".g.AA",SAD_filename =p.grm_file_prefix+".g.AD",SDD_filename =p.grm_file_prefix+".g.DD",\
    SAA_intra_filename =p.grm_file_prefix+".g.AA-intra",SAD_intra_filename =p.grm_file_prefix+".g.AD-intra",SDD_intra_filename =p.grm_file_prefix+".g.DD-intra",\
    SAA_inter_filename =p.grm_file_prefix+".g.AA-inter",SAD_inter_filename =p.grm_file_prefix+".g.AD-inter",SDD_inter_filename =p.grm_file_prefix+".g.DD-inter",\
    SAAA_filename =p.grm_file_prefix+".g.AAA",SAAD_filename =p.grm_file_prefix+".g.AAD",SADD_filename =p.grm_file_prefix+".g.ADD",\
    SDDD_filename =p.grm_file_prefix+".g.DDD";
    string g_diag_filename=p.grm_file_prefix+".gdiag",indgeno_filename=p.grm_file_prefix+".indgeno";
    time_a_geno = clock();
    MatrixXi geno;
    MatrixXd geno_stat;
    read_all_genofiles(p,geno,geno_stat,flog);
//    if(p.hap_f==true) flog<<"Haplotype genotypes file                           :\t"<<p.haplotype<<endl;
//    flog<<"The prefix of GRM files                            :\t"<<p.grm_file_prefix<<endl;
//    flog<<"The number of used threads is                      :\t"<<p.numThreads<<endl<<endl;
    flog<<"The number of SNPs is                              :\t"<<mrk_n<<endl;
    time_b_geno = clock();
//    fprintf(stderr,"Reading genotypes use %d seconds.\n",(int)(time_b_geno-time_a_geno));
    flog<<"The time (seconds) cost for reading genotypes is   :\t"<<(int)((double)(time_b_geno-time_a_geno)/CLOCKS_PER_SEC)<<endl;
    genostat_update(geno_stat);
    get_het_sum(geno_stat,is_hwd,&hetero_total,&hetero_total2);
    int block_size=20000;
    save_wt(p,geno,geno_stat,block_size);
    VectorXd sgdiag;
    sgdiag=(p.chrsplit_f)?VectorXd::Zero(8):VectorXd::Zero(9);
    order=3;
    /////////////
    vector<MatrixXd>WW_A_all;
    vector<MatrixXd>WW_D_all;
    vector<vector<MatrixXd>>WW_AD_all((order+1));
    for(i=0;i<=order;i++) WW_AD_all[i]=vector<MatrixXd>((order+1));
    MatrixXd W_init;
    W_init=  MatrixXd::Constant(ind_n,ind_n,1.0);
    WW_A_all.push_back(W_init);
    WW_D_all.push_back(W_init);
    WW_AD_all[0][0]=W_init;
    W_init.resize(0,0);
    jiangyong_wwt_mk(WW_A_all,WW_D_all,WW_AD_all,geno,geno_stat,is_hwd,order,mrk_n,block_size);
    MatrixXd WAAWAAT_intra,WADWADT_intra,WDDWDDT_intra;
    if(p.chrsplit_f){
        MatrixXd WAAWAAT_intra_tmp,WADWADT_intra_tmp,WDDWDDT_intra_tmp;
        int chr,startp,genol;
        chr=1;
        for (auto chri:mychrindexmap){
            startp=chri.second[0];
            genol=chri.second[1]-startp+1;
            if(chr==1){
                jiangyong_chrsplit_wwt(WAAWAAT_intra,WADWADT_intra,WDDWDDT_intra,geno,geno_stat,startp,genol,block_size,is_hwd);
            }
            else{
                jiangyong_chrsplit_wwt(WAAWAAT_intra_tmp,WADWADT_intra_tmp,WDDWDDT_intra_tmp,geno,geno_stat,startp,genol,block_size,is_hwd);
                WAAWAAT_intra.noalias()+=WAAWAAT_intra_tmp;
                WADWADT_intra.noalias()+=WADWADT_intra_tmp;
                WDDWDDT_intra.noalias()+=WDDWDDT_intra_tmp;
            }
            chr+=1;
        }
        WAAWAAT_intra_tmp.resize(0,0);
        WADWADT_intra_tmp.resize(0,0);
        WDDWDDT_intra_tmp.resize(0,0);
    }
    geno.resize(0,0);
    geno_stat.resize(0,0);
    double ka,kd,kaa,kad,kdd,kaaa,kaad,kadd,kddd;
    MatrixXd WAWAT,WDWDT,WAAWAAT,WADWADT,WDDWDDT,WAAAWAAAT,WAADWAADT,WADDWADDT,WDDDWDDDT;
    WAWAT=WW_A_all[1];
    ka=WAWAT.diagonal().mean();
    sgdiag(0)=sqrt(ka);
    WAWAT/=ka;
    write_binary(SA_filename.c_str(),WAWAT);
//    saveData(SA_filename.c_str(),WAWAT);
    SA_filename.clear();
    cerr<<"\n\n*******************SA************\n\n";
    cerr<<setprecision(6)<<WAWAT.block(0,0,5,5)<<endl;
    WAWAT.resize(0,0);

    WDWDT=WW_D_all[1];
    kd=WDWDT.diagonal().mean();
    sgdiag(1)=sqrt(kd);
    WDWDT/=kd;
    write_binary(SD_filename.c_str(),WDWDT);
    SD_filename.clear();
    cerr<<"\n\n*******************SD************\n\n";
    cerr<<setprecision(6)<<WDWDT.block(0,0,5,5)<<endl;
    WDWDT.resize(0,0);
    int ordera,orderd;
    double kaa_intra,kad_intra,kdd_intra,kaa_inter,kad_inter,kdd_inter;
    MatrixXd WAAWAAT_inter,WADWADT_inter,WDDWDDT_inter;
    WAAWAAT=jiangyong_recurs(WW_A_all,2);
    if(p.chrsplit_f){
        WAAWAAT_inter=WAAWAAT-WAAWAAT_intra;
        kaa_inter=WAAWAAT_inter.diagonal().mean();
        sgdiag(3)=sqrt(kaa_inter);
        WAAWAAT_inter/=kaa_inter;
        write_binary(SAA_inter_filename.c_str(),WAAWAAT_inter);
//        saveData(SAA_inter_filename.c_str(),WAAWAAT_inter);
        SAA_inter_filename.clear();
        cerr<<"\n\n*******************SAA_inter***********\n"<<endl;
        cerr<<setprecision(6)<<WAAWAAT_inter.block(0,0,5,5)<<endl;
        WAAWAAT_inter.resize(0,0);
        kaa_intra=WAAWAAT_intra.diagonal().mean();
        sgdiag(2)=sqrt(kaa_intra);
        WAAWAAT_intra/=kaa_intra;
        write_binary(SAA_intra_filename.c_str(),WAAWAAT_intra);
//        saveData(SAA_intra_filename.c_str(),WAAWAAT_intra);
        SAA_intra_filename.clear();
        cerr<<"\n\n*******************SAA_intra************\n"<<endl;
        cerr<<setprecision(6)<<WAAWAAT_intra.block(0,0,5,5)<<endl;
        WAAWAAT_intra.resize(0,0);
    }
    else{
        kaa=WAAWAAT.diagonal().mean();
        sgdiag(2)=sqrt(kaa);
        WAAWAAT/=kaa;
        write_binary(SAA_filename.c_str(),WAAWAAT);
        SAA_filename.clear();
        cerr<<"\n\n*******************SAA************\n"<<endl;
        cerr<<setprecision(6)<<WAAWAAT.block(0,0,5,5)<<endl;
        WAAWAAT.resize(0,0);

        WAAAWAAAT=jiangyong_recurs(WW_A_all,order);
        kaaa=WAAAWAAAT.diagonal().mean();
        sgdiag(5)=sqrt(kaaa);
        WAAAWAAAT/=kaaa;
        write_binary(SAAA_filename.c_str(),WAAAWAAAT);
        SAAA_filename.clear();
        cerr<<"\n\n*******************SAAA************\n"<<endl;
        cerr<<setprecision(6)<<WAAAWAAAT.block(0,0,5,5)<<endl;
        WAAAWAAAT.resize(0,0);
    }
    ordera=1;
    orderd=1;
    WADWADT=jiangyong_recurs_2(WW_AD_all,ordera,orderd);
    if(p.chrsplit_f){
        WADWADT_inter=WADWADT-WADWADT_intra;
        kad_inter=WADWADT_inter.diagonal().mean();
        sgdiag(5)=sqrt(kad_inter);
        WADWADT_inter/=kad_inter;
        write_binary(SAD_inter_filename.c_str(),WADWADT_inter);
        SAD_inter_filename.clear();
        cerr<<"\n\n*******************SAD_inter************\n"<<endl;
        cerr<<setprecision(6)<<WADWADT_inter.block(0,0,5,5)<<endl;
        WADWADT_inter.resize(0,0);
        kad_intra=WADWADT_intra.diagonal().mean();
        sgdiag(4)=sqrt(kad_intra);
        WADWADT_intra/=kad_intra;
        write_binary(SAD_intra_filename.c_str(),WADWADT_intra);
//        saveData(SAD_intra_filename.c_str(),WADWADT_intra);
        SAD_intra_filename.clear();
        cerr<<"\n\n*******************SAD_intra************\n"<<endl;
        cerr<<setprecision(6)<<WADWADT_intra.block(0,0,5,5)<<endl;
        WADWADT_intra.resize(0,0);
    }
    else{
        kad=WADWADT.diagonal().mean();
        sgdiag(3)=sqrt(kad);
        WADWADT/=kad;
        write_binary(SAD_filename.c_str(),WADWADT);
        SAD_filename.clear();
        cerr<<"\n\n*******************SAD************\n"<<endl;
        cerr<<setprecision(6)<<WADWADT.block(0,0,5,5)<<endl;
        WADWADT.resize(0,0);

        ordera=2;
        orderd=1;
        WAADWAADT=jiangyong_recurs_2(WW_AD_all,ordera,orderd);
        kaad=WAADWAADT.diagonal().mean();
        sgdiag(6)=sqrt(kaad);
        WAADWAADT/=kaad;
        write_binary(SAAD_filename.c_str(),WAADWAADT);
        SAAD_filename.clear();
        cerr<<"\n\n*******************SAAD************\n"<<endl;
        cerr<<setprecision(6)<<WAADWAADT.block(0,0,5,5)<<endl;
        WAADWAADT.resize(0,0);

        ordera=1;
        orderd=2;
        WADDWADDT=jiangyong_recurs_2(WW_AD_all,ordera,orderd);
        kadd=WADDWADDT.diagonal().mean();
        sgdiag(7)=sqrt(kadd);
        WADDWADDT/=kadd;
        write_binary(SADD_filename.c_str(),WADDWADDT);
        SADD_filename.clear();
        cerr<<"\n\n*******************SADD************\n"<<endl;
        cerr<<setprecision(6)<<WADDWADDT.block(0,0,5,5)<<endl;
        WADDWADDT.resize(0,0);
    }
    WDDWDDT=jiangyong_recurs(WW_D_all,2);
    if(p.chrsplit_f){
        WDDWDDT_inter=WDDWDDT-WDDWDDT_intra;
        kdd_inter=WDDWDDT_inter.diagonal().mean();
        sgdiag(7)=sqrt(kdd_inter);
        WDDWDDT_inter/=kdd_inter;
        write_binary(SDD_inter_filename.c_str(),WDDWDDT_inter);
        SDD_inter_filename.clear();
        cerr<<"\n\n*******************SDD_inter************\n"<<endl;
        cerr<<setprecision(6)<<WDDWDDT_inter.block(0,0,5,5)<<endl;
        WDDWDDT_inter.resize(0,0);
        kdd_intra=WDDWDDT_intra.diagonal().mean();
        sgdiag(6)=sqrt(kdd_intra);
        WDDWDDT_intra/=kdd_intra;
        write_binary(SDD_intra_filename.c_str(),WDDWDDT_intra);
        SDD_intra_filename.clear();
        cerr<<"\n\n*******************SDD_intra************\n"<<endl;
        cerr<<setprecision(6)<<WDDWDDT_intra.block(0,0,5,5)<<endl;
        WDDWDDT_intra.resize(0,0);
    }
    else{
        kdd=WDDWDDT.diagonal().mean();
        sgdiag(4)=sqrt(kdd);
        WDDWDDT/=kdd;
        write_binary(SDD_filename.c_str(),WDDWDDT);
        SDD_filename.clear();
        cerr<<"\n\n*******************SDD************\n"<<endl;
        cerr<<setprecision(6)<<WDDWDDT.block(0,0,5,5)<<endl;
        WDDWDDT.resize(0,0);

        WDDDWDDDT=jiangyong_recurs(WW_D_all,order);
        kddd=WDDDWDDDT.diagonal().mean();
        sgdiag(8)=sqrt(kddd);
        WDDDWDDDT/=kddd;
        write_binary(SDDD_filename.c_str(),WDDDWDDDT);
        SDDD_filename.clear();
        cerr<<"\n\n*******************SDDD************\n"<<endl;
        cerr<<setprecision(6)<<WDDDWDDDT.block(0,0,5,5)<<endl;
        WDDDWDDDT.resize(0,0);
    }
    write_binary(g_diag_filename.c_str(),sgdiag);
    flog<<"The method used for epistatic GRMs inference is    :\tEGRM method"<<endl;
    for (i=0;i<sgdiag.size();i++) {
        flog<<setprecision(6)<<scientific<<"The mean of the diagonal elements for K"<<var_name[i];
        if(var_name[i].size()==1){
            flog<<" is        ";
        }
        else if (var_name[i].size()==2){
            flog<<" is       ";
        }
        else if (var_name[i].size()==3){
            flog<<" is      ";
        }
        else{
            flog<<" is ";
        }
        flog<<":\t"<<pow(sgdiag[i],2)<<endl;
    }
    sgdiag.resize(0);
    for(i=0;i<WW_A_all.size();i++)  WW_A_all[i].resize(0,0);
    for(i=0;i<WW_D_all.size();i++)  WW_D_all[i].resize(0,0);
    for(i=0;i<WW_AD_all.size();i++){
        for(j=0;j<WW_AD_all[i].size();j++) WW_AD_all[i][j].resize(0,0);
    }
    time_b = clock();
    flog<<"The time (seconds) cost for SNP GRMs inference is  :\t"<<(int)((double)(time_b-time_a)/CLOCKS_PER_SEC)<<endl;
    ifstream input_indgeno(indgeno_filename);
    if(!input_indgeno) save_ind_geno(indgeno_filename,&ind_geno);
    input_indgeno.close();
}
void Greml_ce::load_genotype(const parameters &p,ofstream &flog){
    int i,var_vec_size_calc,count=0;
    if(p.chrsplit_f){
        var_vec_size_calc=p.additive_variances_f+p.dominance_variances_f+\
        p.additive_additive_intra_variances_f+p.additive_additive_inter_variances_f+\
        p.additive_dominance_intra_variances_f+p.additive_dominance_inter_variances_f+\
        p.dominance_dominance_intra_variances_f+p.dominance_dominance_inter_variances_f+\
        p.additive_haplotype_variances_f+p.residual_variances_f;
    }
    else{
        var_vec_size_calc=p.additive_variances_f+p.dominance_variances_f+\
        p.additive_additive_variances_f+p.additive_dominance_variances_f+p.dominance_dominance_variances_f+\
        p.additive_additive_additive_variances_f+p.additive_additive_dominance_variances_f+\
        p.additive_dominance_dominance_variances_f+p.dominance_dominance_dominance_variances_f+\
        p.additive_haplotype_variances_f+p.residual_variances_f;
    }
    for (var_name_flags_iter = var_name_flags.begin(); var_name_flags_iter != var_name_flags.end(); var_name_flags_iter++){
        if(var_name_flags_iter->second==0) var_name_values.erase(var_name_flags_iter->first);
    }
    vector<string>var_name_new;
    vector<double>var_new;
    for (string elem : var_name){
        auto it = var_name_values.find(elem);
        if (it != var_name_values.end()){
            var_name_new.push_back(it->first);
            var_new.push_back(it->second);
        }
    }
    var.conservativeResize(var_vec_size_calc);
    var = Map<VectorXd>(var_new.data(), var_new.size());
//    for (i=0; i<var.size();i++) cout<<var(i)<<","<<var_name_new[i]<<endl;
	MatrixXd gdiag1;
    string gdiag_filename = p.grm_file_prefix+".gdiag";
	read_binary(gdiag_filename.c_str(),gdiag1);
	G_diag_mean_sqrt=VectorXd::Zero(var_vec_size_calc-1);
	count=0;
	for (i=0;i<var_name_new.size()-1;i++){
	    if(var_name_new[i]!="AH"){
	        G_diag_mean_sqrt(count)=gdiag1(i,0);
	        count++;
	    }
	}
	var_name=var_name_new;
	if (var_name[var_vec_size_calc-1]!="E"){
	    fprintf(stderr,"\nError: You forget to set the initial residual variance value!\n\n");
	    throw;
	}else{
	    fprintf(stderr,"\nLoading GRMs from files with prefix %s!\n",p.grm_file_prefix.c_str());
	    flog<<"Loading GRMs from files with prefix:\t"<<p.grm_file_prefix<<endl;
	    MatrixXd g_tmp;
	    for (i=0; i<var_name.size()-1;i++){
	        string gmat_filename =p.grm_file_prefix+".g."+var_name[i];
            read_binary(gmat_filename.c_str(),g_tmp);
	        G.push_back(g_tmp);
        }
        g_tmp.resize(0,0);
	    string indgeno_filename=p.grm_file_prefix+".indgeno";
	    load_ind_geno(indgeno_filename,&ind_geno);
	}
/////////Logout file///////////////////////////////////////
//    flog<<"The prefix of GRM files                            :\t"<<p.grm_file_prefix<<endl;
//    flog<<"Phenotypes file                                    :\t"<<p.phenotype<<endl;
//	flog<<"The column number for the trait to be analyzed     :\t"<<p.trait_pos<<endl;
//	flog<<"The missing value is                               :\t"<<p.missing_phen_val<<endl;
//	flog<<"The prefix of output files is                      :\t"<<p.outprefix<<endl;
//	flog<<"The number of threads is                           :\t"<<p.numThreads<<endl;
//	if(p.cin_var_f){
//	    for(i=0;i<var_name.size();i++) {
//            flog<<"The value of V"<<var_name[i];
//            if(var_name[i].size()==1){
//                flog<<" is                                ";
//            }
//            else if (var_name[i].size()==2){
//                flog<<" is                               ";
//            }
//            else if (var_name[i].size()==3){
//                flog<<" is                              ";
//            }
//            else{
//                flog<<" is                         ";
//            }
//            flog<<":\t"<<var_name_values[var_name[i]]<<endl;
//        }
//	}
//	else{
//	    flog<<"The number of iterations is                        :\t"<<p.iter_n<<endl;
//	    flog<<"The iteration number for starting AI-REML is       :\t"<<p.iter_ai_start<<endl;
//	    for(i=0;i<var_name.size();i++) {
//            flog<<"The initial value of V"<<var_name[i];
//            if(var_name[i].size()==1){
//                flog<<" is                         ";
//            }
//            else if (var_name[i].size()==2){
//                flog<<" is                        ";
//            }
//            else if (var_name[i].size()==3){
//                flog<<" is                       ";
//            }
//            else{
//                flog<<" is                  ";
//            }
//            flog<<":\t"<<var_name_values[var_name[i]]<<endl;
//            }
//        flog<<"The variance component tolerance is                :\t"<<p.iter_tolerance<<endl;
//        flog<<"The heritability tolerance is                      :\t"<<p.iter_tolerance_her<<endl;
//	}
}