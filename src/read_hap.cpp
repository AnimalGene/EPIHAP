#include "read_geno.h"
void Greml_ce::each_W_hap_generate(int &block,MatrixXi &geno_hap,MatrixXd &W_hap,hap_block_allele &hb_allele,int n_row,int hap_max_key_j){
    int i,k,hap_snp1,hap_snp2;
    hap_snp1=block*2;
    hap_snp2=hap_snp1+1;
#pragma omp parallel for private(i)
    for(i=0;i<n_row;i++){
        if(geno_hap(i,hap_snp1)==0 || geno_hap(i,hap_snp2)==0){
            for(k=0;k<hb_allele[block].size()-1;k++){
                W_hap(i,k)=0.0;
            }
        }
        else{
            for(auto it=hb_allele[block].begin();it!=hb_allele[block].end();++it){
                if(it->second.idx >= 0){
                    W_hap(i,(it->second.idx))=(it->second.count)*2.0;
                }
            }
            if(geno_hap(i,hap_snp1)==hap_max_key_j && geno_hap(i,hap_snp2)==hap_max_key_j){
                ;
            }
            else if(geno_hap(i,hap_snp1)!=hap_max_key_j && geno_hap(i,hap_snp2)==hap_max_key_j){
                W_hap(i,hb_allele[block][geno_hap(i,hap_snp1)].idx)-=1.0;
            }
            else if(geno_hap(i,hap_snp1)==hap_max_key_j && geno_hap(i,hap_snp2)!=hap_max_key_j){
                W_hap(i,hb_allele[block][geno_hap(i,hap_snp2)].idx)-=1.0;
            }
            else if(geno_hap(i,hap_snp1)!=hap_max_key_j && geno_hap(i,hap_snp2)!=hap_max_key_j){
                W_hap(i,hb_allele[block][geno_hap(i,hap_snp1)].idx)-=1.0;
                W_hap(i,hb_allele[block][geno_hap(i,hap_snp2)].idx)-=1.0;
            }
        }
    }
}
void Greml_ce::read_hap(const parameters &p,ofstream &flog){
    clock_t time_a, time_b;
    time_a = clock();
    int i,j,ii,hapflag,n_2_hapcol,n_hapcol,b,start_b,count1,b1;
    int geno_tmptmp_0,geno_tmptmp_1;
    double hetero_total  = 0.0,hetero_total2 = 0.0,alpha = 1.0,beta  = 0.0;
    string hap_input_file=p.haplotype,haplotype_matrix_file=p.grm_file_prefix+".g.AH";
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
    flog<<"\nThe number of haplotype blocks is                  :\t"<<n_hapcol<<endl;
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
            geno_tmptmp_0             = atoi(hapwords[i*2+1].c_str());
            geno_tmptmp_1             = atoi(hapwords[i*2+2].c_str());
            geno_hap(hapflag,i*2)     = (geno_tmptmp_0 == p.missing_hap_val)?0:geno_tmptmp_0;
            geno_hap(hapflag,(i*2+1)) = (geno_tmptmp_1 == p.missing_hap_val)?0:geno_tmptmp_1;
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
    vector<MatrixXd>Ghapp;
    MatrixXd Ghap_part;
    MatrixXd W_hap;
    MatrixXd Gmat_AH;
    int hap_count,hap_max_key_j;
    for(b1=0;b1<n_hapcol;b1++){
        hap_count = hb_allele[b1].size()-1;
        W_hap.resize(ind_n,hap_count);
        hap_max_key_j=hap_max_key(b1);
        each_W_hap_generate(b1,geno_hap,W_hap,hb_allele,ind_n,hap_max_key_j);
        Ghap_part.noalias() = MatrixXd::Zero(ind_n,ind_n);
        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,ind_n,ind_n,hap_count,alpha,W_hap.data(),ind_n,W_hap.data(), ind_n,beta, Ghap_part.data(), ind_n);
        if (b1==0){
            Gmat_AH=Ghap_part;
        }
        else{
            Gmat_AH.noalias()+=Ghap_part;
        }
    }
    Ghap_part.resize(0,0);
    double kah;
    kah=Gmat_AH.diagonal().mean();
    Gmat_AH/=kah;
    write_binary(haplotype_matrix_file.c_str(),Gmat_AH);
    haplotype_matrix_file.clear();
    cerr<<"\n*******************SAH************\n"<<endl;
    cerr<<Gmat_AH.block(0,0,5,5)<<endl;
    flog<<setprecision(6)<<scientific<<"The mean of the diagonal elements for KAH is       :\t"<<kah<<endl;
    Gmat_AH.resize(0,0);
    time_b = clock();
    flog<<"The time (seconds) cost for AH GRMs inference is   :\t"<<(int)((double)(time_b-time_a)/CLOCKS_PER_SEC)<<endl;
}