#include "read_geno.h"
void Greml_ce::jiangyong_egrm_a1(MatrixXd &W,int pow_set, int block_size,MatrixXd &WW){
    MatrixXd W_new;
    double alpha = 1.0,beta  = 0.0;
    W_new=W.array().pow(pow_set).matrix();
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,ind_n,ind_n,block_size,alpha,W_new.data(),ind_n,W_new.data(),ind_n,beta,WW.data(),ind_n);
    W_new.resize(0,0);
}
void Greml_ce::jiangyong_egrm_ad(MatrixXd &WA,MatrixXd &WD,int pow_seta, int pow_setd,int block_size,MatrixXd &WW){
    MatrixXd WA_new;
    MatrixXd WD_new;
    MatrixXd WAD_new;
    int colw,roww;
    colw=WA.cols();
    roww=WA.rows();
    double alpha = 1.0,beta  = 0.0;
    if(pow_seta==0){
        WA_new=MatrixXd::Constant(roww,colw,1.0);
    }
    else{
        WA_new=WA.array().pow(pow_seta).matrix();
    }
    WD_new=WD.array().pow(pow_setd).matrix();
    WAD_new=WA_new.cwiseProduct(WD_new);
    WA_new.resize(0,0);
    WD_new.resize(0,0);
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,ind_n,ind_n,block_size,alpha,WAD_new.data(),ind_n,WAD_new.data(),ind_n,beta,WW.data(),ind_n);
    WAD_new.resize(0,0);
}
MatrixXd Greml_ce::jiangyong_recurs(vector<MatrixXd>&WWT,int degree){
    MatrixXd WWepi_mtd3;
    double dgr=(double) degree;
    double val;
    int i;
    if(degree==0){
        WWepi_mtd3=MatrixXd::Constant(ind_n,ind_n,1.0);
        return WWepi_mtd3;
    }
    else{
        WWepi_mtd3= MatrixXd::Zero(ind_n,ind_n);
        for(i=1;i<=degree;i++){
            val=(1/dgr)*pow((-1),(i-1));
            WWepi_mtd3=WWepi_mtd3+val*jiangyong_recurs(WWT,(degree-i)).cwiseProduct(WWT[i]);
        }
        return WWepi_mtd3;
    }
}
MatrixXd Greml_ce::jiangyong_recurs_2(vector<vector<MatrixXd>>&WWT,int degreea, int degreed){
    MatrixXd WWepi_mtd3;
    double dgra=(double) degreea;
    double dgrd=(double) degreed;
    double val;
    int i,j,a,d;
    if((degreea==0)&&(degreed==0)){
        WWepi_mtd3=MatrixXd::Constant(ind_n,ind_n,1.0);
        return WWepi_mtd3;
    }
    else{
        WWepi_mtd3= MatrixXd::Zero(ind_n,ind_n);
        if(degreed==0){
            a=1;
        }
        else{
            a=0;
        }
        for(i=a;i<=degreea;i++){
            if(i==0){
                if(degreed!=0) d=1;
            }
            else{
                d=0;
            }
            for(j=d;j<=degreed;j++){
                val=(1/(dgra+dgrd))*pow((-1),(i+j-1))*binomialCoefficients((i+j),i);
                WWepi_mtd3=WWepi_mtd3+val*jiangyong_recurs_2(WWT,(degreea-i),(degreed-j)).cwiseProduct(WWT[i][j]);
            }
        }
        return WWepi_mtd3;
    }
}
void Greml_ce::jiangyong_wwt_mk(vector<MatrixXd> &WW_A_all,vector<MatrixXd> &WW_D_all,vector<vector<MatrixXd>> &WW_AD_all,MatrixXi &geno,MatrixXd &geno_stat,int is_hwd,int order,int allmarkn,int block_size){
    MatrixXd WA;
    MatrixXd WD;
    MatrixXd WW_A;
    MatrixXd WW_D;
    MatrixXd WW_AD;
    int j,k,b,start_b;
    int remainders = allmarkn % block_size;
    if(allmarkn>block_size){
        for (b=block_size;b<=allmarkn;b+=block_size){
            WA = MatrixXd::Zero(ind_n,block_size);
            WD = MatrixXd::Zero(ind_n,block_size);
            geno_transform(geno.block(0,(b-block_size),ind_n,block_size),geno_stat.block(0,(b-block_size),12,block_size),var_name[0],is_hwd,WA);
            geno_transform(geno.block(0,(b-block_size),ind_n,block_size),geno_stat.block(0,(b-block_size),12,block_size),var_name[1],is_hwd,WD);
            for(j=0;j<=order;j++){
                if(j>0){
                    WW_A.noalias() = MatrixXd::Zero(ind_n,ind_n);
                    WW_D.noalias() = MatrixXd::Zero(ind_n,ind_n);
                    jiangyong_egrm_a1(WA,j,block_size,WW_A);
                    jiangyong_egrm_a1(WD,j,block_size,WW_D);
                    if((b-block_size)==0){
                        WW_A_all.push_back(WW_A);
                        WW_D_all.push_back(WW_D);
                        WW_AD_all[j][0]=WW_A;
                    }
                    else{
                        WW_A_all[j].noalias()+=WW_A;
                        WW_D_all[j].noalias()+=WW_D;
                        WW_AD_all[j][0].noalias()+=WW_A;
                    }
                }
                for(k=1;k<=order;k++){
                    WW_AD.noalias() = MatrixXd::Zero(ind_n,ind_n);
                    jiangyong_egrm_ad(WA,WD,j,k,block_size,WW_AD);
                    if((b-block_size)==0){
                        WW_AD_all[j][k]=WW_AD;
                    }
                    else{
                        WW_AD_all[j][k].noalias()+=WW_AD;
                    }
                }
            }
        }
        if(remainders>0){
            start_b=allmarkn-remainders;
            WA = MatrixXd::Zero(ind_n,remainders);
            WD = MatrixXd::Zero(ind_n,remainders);
            geno_transform(geno.block(0,start_b,ind_n,remainders),geno_stat.block(0,start_b,12,remainders),var_name[0],is_hwd,WA);
            geno_transform(geno.block(0,start_b,ind_n,remainders),geno_stat.block(0,start_b,12,remainders),var_name[1],is_hwd,WD);
            for(j=0;j<=order;j++){
                if(j>0){
                    WW_A.noalias() = MatrixXd::Zero(ind_n,ind_n);
                    WW_D.noalias() = MatrixXd::Zero(ind_n,ind_n);
                    jiangyong_egrm_a1(WA,j,remainders,WW_A);
                    jiangyong_egrm_a1(WD,j,remainders,WW_D);
                    WW_A_all[j].noalias()+=WW_A;
                    WW_D_all[j].noalias()+=WW_D;
                    WW_AD_all[j][0].noalias()+=WW_A;
                }
                for(k=1;k<=order;k++){
                    WW_AD.noalias() = MatrixXd::Zero(ind_n,ind_n);
                    jiangyong_egrm_ad(WA,WD,j,k,remainders,WW_AD);
                    WW_AD_all[j][k].noalias()+=WW_AD;
                }
            }
        }
        else{
            ;
        }
    }
    else{
        WA = MatrixXd::Zero(ind_n,allmarkn);
        WD = MatrixXd::Zero(ind_n,allmarkn);
        geno_transform(geno,geno_stat,var_name[0],is_hwd,WA);
        geno_transform(geno,geno_stat,var_name[1],is_hwd,WD);
        for(j=0;j<=order;j++){
            if(j>0){
                WW_A.noalias() = MatrixXd::Zero(ind_n,ind_n);
                WW_D.noalias() = MatrixXd::Zero(ind_n,ind_n);
                jiangyong_egrm_a1(WA,j,allmarkn,WW_A);
                jiangyong_egrm_a1(WD,j,allmarkn,WW_D);
                WW_A_all.push_back(WW_A);
                WW_D_all.push_back(WW_D);
                WW_AD_all[j][0]=WW_A;
            }
            for(k=1;k<=order;k++){
                WW_AD.noalias() = MatrixXd::Zero(ind_n,ind_n);
                jiangyong_egrm_ad(WA,WD,j,k,allmarkn,WW_AD);
                WW_AD_all[j][k]=WW_AD;
            }

        }
    }
    WA.resize(0,0);
    WD.resize(0,0);
    WW_A.resize(0,0);
    WW_D.resize(0,0);
    WW_AD.resize(0,0);
}
void Greml_ce::henderson_chrsplit_WWTAA_WWTDD(MatrixXd &WWTepi,MatrixXi &geno,MatrixXd &geno_stat,string varname,int is_hwd){
    /////////////
    MatrixXd W,WWT,genochr_stat;
    MatrixXi genochr;
    double k;
    double alpha = 1.0,beta  = 0.0;
    int chr,startp,genol;
    chr=1;
    MatrixXd WWTepitmp;
    for (auto chri:mychrindexmap){
        startp=chri.second[0];
        genol=chri.second[1]-startp+1;
        genochr=geno.block(0,startp,ind_n,genol);
        genochr_stat=geno_stat.block(0,startp,12,genol);
        W = MatrixXd::Zero(ind_n,genol);
        geno_transform(genochr,genochr_stat,varname,is_hwd,W);
        WWT.noalias() = MatrixXd::Zero(ind_n,ind_n);
        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,ind_n,ind_n,genol,alpha,W.data(),ind_n,W.data(),ind_n,beta,WWT.data(),ind_n);
        if(chr==1){
            WWTepi=WWT.cwiseProduct(WWT);
        }
        else{
            WWTepitmp=WWT.cwiseProduct(WWT);
            WWTepi.noalias()+=WWTepitmp;
        }
        chr+=1;

    }
    genochr.resize(0,0);
    genochr_stat.resize(0,0);
    W.resize(0,0);
    WWT.resize(0,0);
    WWTepitmp.resize(0,0);
}
void Greml_ce::jiangyong_chrsplit_wwt(MatrixXd &WAAWAAT_intra,MatrixXd &WADWADT_intra,MatrixXd &WDDWDDT_intra,MatrixXi &geno,MatrixXd &geno_stat,int startp,int genol,int block_size,int is_hwd){
    /////////////
    MatrixXi genochr;
    MatrixXd genochr_stat;
    int i,j,ordera,orderd,orderc;
    ordera=1;
    orderd=1;
    orderc=2;
    vector<MatrixXd>WW_A_all;
    vector<MatrixXd>WW_D_all;
    vector<vector<MatrixXd>>WW_AD_all((orderc+1));
    for(i=0;i<=orderc;i++) WW_AD_all[i]=vector<MatrixXd>((orderc+1));
    MatrixXd W_init;
    W_init=  MatrixXd::Constant(ind_n,ind_n,1.0);
    WW_A_all.push_back(W_init);
    WW_D_all.push_back(W_init);
    WW_AD_all[0][0]=W_init;
    W_init.resize(0,0);
    genochr=geno.block(0,startp,ind_n,genol);
    genochr_stat=geno_stat.block(0,startp,12,genol);
    jiangyong_wwt_mk(WW_A_all,WW_D_all,WW_AD_all,genochr,genochr_stat,is_hwd,orderc,genol,block_size);
    genochr.resize(0,0);
    genochr_stat.resize(0,0);
//A
    WAAWAAT_intra=jiangyong_recurs(WW_A_all,2);
    WADWADT_intra=jiangyong_recurs_2(WW_AD_all,ordera,orderd);
    WDDWDDT_intra=jiangyong_recurs(WW_D_all,2);
    for(i=0;i<WW_A_all.size();i++)  WW_A_all[i].resize(0,0);
    for(i=0;i<WW_D_all.size();i++)  WW_D_all[i].resize(0,0);
    for(i=0;i<WW_AD_all.size();i++){
        for(j=0;j<WW_AD_all[i].size();j++) WW_AD_all[i][j].resize(0,0);
    }
}
void Greml_ce::henderson_chrsplit_WWTAD(MatrixXd &WWTepi,MatrixXi &geno,MatrixXd &geno_stat,int is_hwd){
    /////////////
    MatrixXd WA,WD,WWTA,WWTD,genochr_stat;
    MatrixXi genochr;
    string varnameA="A",varnameD="D";
    double alpha = 1.0,beta  = 0.0;
    int chr,startp,genol;
    chr=1;
    MatrixXd WWTepitmp;
    for (auto chri:mychrindexmap){
        startp=chri.second[0];
        genol=chri.second[1]-startp+1;
        genochr=geno.block(0,startp,ind_n,genol);
        genochr_stat=geno_stat.block(0,startp,12,genol);
        WA = MatrixXd::Zero(ind_n,genol);
        WD = MatrixXd::Zero(ind_n,genol);
        geno_transform(genochr,genochr_stat,varnameA,is_hwd,WA);
        geno_transform(genochr,genochr_stat,varnameD,is_hwd,WD);
        WWTA.noalias() = MatrixXd::Zero(ind_n,ind_n);
        WWTD.noalias() = MatrixXd::Zero(ind_n,ind_n);
        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,ind_n,ind_n,genol,alpha,WA.data(),ind_n,WA.data(),ind_n,beta,WWTA.data(),ind_n);
        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,ind_n,ind_n,genol,alpha,WD.data(),ind_n,WD.data(),ind_n,beta,WWTD.data(),ind_n);
        if(chr==1){
            WWTepi=WWTA.cwiseProduct(WWTD);
        }
        else{
            WWTepitmp=WWTA.cwiseProduct(WWTD);
            WWTepi.noalias()+=WWTepitmp;
        }
        chr+=1;
    }
    genochr.resize(0,0);
    genochr_stat.resize(0,0);
    WA.resize(0,0);
    WWTA.resize(0,0);
    WD.resize(0,0);
    WWTD.resize(0,0);
    WWTepitmp.resize(0,0);
}
void Greml_ce::matrix_rank_print(MatrixXd &Gmat,string mname){
     ////
    int mrank,n,info;
    double   rcond  = 1e-10;
    n = Gmat.cols();
    MatrixXd Gmat_inv = MatrixXd::Identity(n,n);
    double   *ss     = (double*)malloc(n*n*sizeof(double));
    info = LAPACKE_dgelss(LAPACK_COL_MAJOR,n,n,n,Gmat.data(),n,Gmat_inv.data(),n,ss,rcond,&mrank);
    cerr<<"*******************Gmat_"<<mname<<"_inv************"<<endl<<endl;
    cerr<<Gmat_inv.block(0,0,5,5)<<endl;
    cerr<<endl<<"Rank of GRM "<<mname<<" is : "<<mrank<<" ."<<endl<<endl;
    Gmat_inv.resize(0,0);
    ////
}