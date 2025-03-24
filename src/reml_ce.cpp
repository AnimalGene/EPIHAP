#include "read_geno.h"
void Greml_ce::update_dif_max(VectorXd &var,VectorXd &var_tmp,VectorXd &dif_tmp,double &dif_max){
    for(int i=0;i<var.size();i++){
	    dif_tmp(i) = abs(var[i]-var_tmp[i]);
	    if(dif_max<dif_tmp(i)){
		    dif_max = dif_tmp(i);
		}
    }
}
void Greml_ce::get_h2(MatrixXd &IAI,VectorXd &h2,VectorXd &h2se){
    int i,j;
    double    vp = var.sum();
	h2(var.size()-1)=0.0;
    for(i=0;i<var.size()-1;i++){
        if(var_name[i]!="PE"){
            h2(i) = var(i)/vp;
            h2(var.size()-1) += h2(i);
            VectorXd dl_h2(var.size());
            for(j=0;j<var.size();j++){
		        dl_h2(j)=(i==j)?((1-h2(i))/vp):(-1.0*h2(i)/vp);
		    }
            h2se(i) = sqrt(dl_h2.dot(IAI*dl_h2));
        }else{
            h2(i) =0.0;
            h2se(i)=0.0;
        }
    }
	i = var.size()-1;
	VectorXd dl_h2(var.size());
	for(j=0;j<var.size();j++){
	    dl_h2(j)=(i==j)?(-1.0*h2(i)/vp):((1-h2(i))/vp);
    }
	h2se(i) = sqrt(dl_h2.dot(IAI*dl_h2));
}
void Greml_ce::output_greml_iter_title(ofstream &var_output){
    int i;
	var_output<<setw(10)<<left<<"Iteration";
    for(i=0;i<var.size();i++){
        string name = "V" + var_name[i];
	    var_output <<left<< setw(22) << name;
//	    var_output<<left<<(name.size()<5 ?  setw(30) :  setw(30))<<name;
	    name = "Tolerance_V" + var_name[i];
	    var_output<<left << setw(22) << name;
    }
    var_output << endl;
}
void Greml_ce::output_greml_iter_var(ofstream &var_output,int iter,VectorXd &dif_tmp,ofstream &flog){
    int i;
    var_output << setw(10) << iter+1;
    for(i=0;i<var.size();i++){
        flog<<"V"<<setw(10)<<left<<var_name[i]<<"   = "<<setprecision(6)<<scientific<<setw(15)\
		<<var[i]<<"Tolerance_V"<<setw(10)<<left<<var_name[i]<<"   = "<<setprecision(6)<<scientific<<setw(15)<<dif_tmp[i]<<endl;
        var_output<<left<<setw(22)<<setprecision(6)<<scientific<<var[i];
		var_output<<left<<setw(22)<<setprecision(6)<<scientific<<dif_tmp[i];
    }
    var_output<<endl;
}
void Greml_ce::output_greml_inv_ai(ofstream &var_output,MatrixXd &IAI,ofstream &flog){
    int i,j;
    var_output<<setw(10)<<"SE";
    for(i=0;i<var.size();i++){
	    var_output<<left<<setw(22)<<setprecision(6)<<scientific<<sqrt(IAI(i,i));
	    var_output<<left<<setw(22)<<setprecision(6)<<"";
    }
    flog<<"\nInverse of AI matrix:\n\n";
    for(i=0;i<var.size();i++){
	    for(j=0;j<var.size();j++){
            flog<<setw(15)<<right<<setprecision(6)<<scientific<<IAI(i,j);
	    }
	    flog<<endl;
    }
    flog<<endl;
    var_output<<endl<<endl;
}
void Greml_ce::output_greml_her(VectorXd &h2,VectorXd &h2dif_tmp,ofstream &flog){
    int i;
    flog<<endl;
    for(i=0;i<var.size();i++){
        if(var_name[i]!="PE"){
	        if(var_name[i]=="E"){
			    flog<<"H2"<<setw(11)<<left<<""<<" = "<<setprecision(6)<<scientific<<setw(15)<<h2[i]<<"Tolerance_H2"\
			    <<setw(11)<<left<<""<<" = "<<setprecision(6)<<scientific<<setw(15)<<h2dif_tmp[i]<<endl;
            }
            else{
	            flog<<"h2_"<<setw(10)<<left<<var_name[i]<<" = "<<setprecision(6)<<scientific<<setw(15)<<h2[i]\
			    <<"Tolerance_h2_"<<setw(10)<<left<<var_name[i]<<" = "<<setprecision(6)<<scientific<<setw(15)<<h2dif_tmp[i]<<endl;
            }
        }
    }
}
void Greml_ce::output_greml_her_a(ofstream &var_output,VectorXd &h2,VectorXd &h2se,ofstream &flog){
    int i;
    flog<<endl;
    for(i=0;i<var.size();i++){
        if(var_name[i]!="PE"){
	        if(var_name[i]=="A"){
		    flog<<"Additive heritability, SE                          : ";
      var_output<<"Additive heritability, SE                          : ";
	        }
	        else if(var_name[i]=="D"){
		        flog<<"Dominance heritability, SE                         : ";
		  var_output<<"Dominance heritability, SE                         : ";
	        }
		    else if(var_name[i]=="AA"){
		        flog<<"Additive X Additive heritability, SE               : ";
		  var_output<<"Additive X Additive heritability, SE               : ";
		    }
		    else if(var_name[i]=="AD"){
	            flog<<"Additive X Dominance heritability, SE              : ";
	      var_output<<"Additive X Dominance heritability, SE              : ";
		    }
		    else if(var_name[i]=="DD"){
		        flog<<"Dominance X Dominance heritability, SE             : ";
		  var_output<<"Dominance X Dominance heritability, SE             : ";
		    }
		    else if(var_name[i]=="AA-intra"){
		        flog<<"Additive X Additive intra-chr heritability, SE     : ";
		  var_output<<"Additive X Additive intra-chr heritability, SE     : ";
		    }
		    else if(var_name[i]=="AD-intra"){
	            flog<<"Additive X Dominance intra-chr heritability, SE    : ";
	      var_output<<"Additive X Dominance intra-chr heritability, SE    : ";
		    }
		    else if(var_name[i]=="DD-intra"){
		        flog<<"Dominance X Dominance intra-chr heritability, SE   : ";
		  var_output<<"Dominance X Dominance intra-chr heritability, SE   : ";
		    }
		    else if(var_name[i]=="AA-inter"){
		        flog<<"Additive X Additive inter-chr heritability, SE     : ";
		  var_output<<"Additive X Additive inter-chr heritability, SE     : ";
		    }
		    else if(var_name[i]=="AD-inter"){
	            flog<<"Additive X Dominance inter-chr heritability, SE    : ";
	      var_output<<"Additive X Dominance inter-chr heritability, SE    : ";
		    }
		    else if(var_name[i]=="DD-inter"){
		        flog<<"Dominance X Dominance inter-chr heritability, SE   : ";
		  var_output<<"Dominance X Dominance inter-chr heritability, SE   : ";
		    }
		    else if(var_name[i]=="AAA"){
		        flog<<"Additive X Additive X Additive heritability, SE    : ";
		  var_output<<"Additive X Additive X Additive heritability, SE    : ";
		    }
		    else if(var_name[i]=="AAD"){
	            flog<<"Additive X Additive X Dominance heritability, SE   : ";
	      var_output<<"Additive X Additive X Dominance heritability, SE   : ";
		    }
		    else if(var_name[i]=="ADD"){
		        flog<<"Additive X Dominance X Dominance heritability, SE  : ";
		  var_output<<"Additive X Dominance X Dominance heritability, SE  : ";
		    }
		    else if(var_name[i]=="DDD"){
		        flog<<"Dominance X Dominance X Dominance heritability, SE : ";
		  var_output<<"Dominance X Dominance X Dominance heritability, SE : ";
		    }
		    else if(var_name[i]=="AH"){
	            flog<<"Additive Haplotype heritability, SE                : ";
	      var_output<<"Additive Haplotype heritability, SE                : ";
	        }
		    else if(var_name[i]=="E"){
			    flog<<"Heritability in the broad sense, SE                : ";
	      var_output<<"Heritability in the broad sense, SE                : ";
		    }
	        flog<<h2(i)<< ", "<<h2se(i)<<endl;
	        var_output<<h2(i)<< ", "<<h2se(i)<<endl;
	    }
	}
}
void Greml_ce::variance_component_calculation(const parameters &p,SparseMatrix<double> &Xt,SparseMatrix<double> &Zt,int ind_n,int ve_i,ofstream &flog){
    int i,j,m,n;
    clock_t time_a, time_b,time_a_tmp,time_b_tmp;
    string var_out_file=p.outprefix+"_greml.txt";
    ofstream var_output(var_out_file.c_str());
    output_greml_iter_title(var_output);
//    cerr<<"\nThe initial value for each type of variance components:"<<endl;
//    cerr<<"V"<<var_name[i]<<"\t"<<(var_name[i].size()<5 ? "\t" : "")<<var[i]<<endl;
    const double alpha_tol = 0.00001;
    const double var_tol   = 0.00001;
    double alpha = 1.0;
	double beta = 0.0;
    VectorXd  var_tmp;
    VectorXd  var_ai;
    VectorXd  h2_tmp(var.size());
    int       *ipiv;
    int       lda;
    int       info;
    MatrixXd  V;
    MatrixXd  IVX;
    MatrixXd  XIVX;
    VectorXd  PY;
    VectorXd  ZPY;
    VectorXd  Z1PY;
    MatrixXd  PZGZ;
    MatrixXd Identy;
    MatrixXd Z1t;
    if(p.modelpe_f){
        int z1cols;
        z1cols=Z1.cols();
        Identy=MatrixXd::Identity(z1cols,z1cols);
        Z1t=Z1.transpose();
    }
    MatrixXd  IAI(var.size(),var.size());
    time_a = clock();
    if(p.iter_n>0) fprintf(stderr,"\n======== GREML begins ======\n\n");
    MatrixXd  G_tmp;
    //==============Iteriation start=======================================================
    for(int iter=0;iter<p.iter_n;iter++){
        flog<<"\nIteration: "<<iter+1<<endl;
        fprintf(stderr,"Iteration %d finished!\n",iter+1);
        clock_t time_a_iter;
        time_a_iter = clock();
        var_tmp = var;
        var_ai  = var;
//=================P matrix begin=======================================================
        if(p.modelpe_f){
            MatrixXd Gpe;
            MatrixXd GpeZ1_tmp;
            G_tmp = MatrixXd::Zero(ind_n,ind_n);
            for(i=0;i<ve_i;i++){
                if(var_name[i]=="PE"){
                    Gpe=var_tmp[i]*Identy;
                }else{
                    G_tmp.noalias()+=var_tmp[i]*G[i];
                }
            }
            MatrixXd GZ_tmp = G_tmp*Zt;
            G_tmp.resize(0,0);
            V.noalias() = Z * GZ_tmp;
            GZ_tmp.resize(0,0);
            GpeZ1_tmp=Gpe*Z1t;
            Gpe.resize(0,0);
            V.noalias()+=Z1*GpeZ1_tmp;
            GpeZ1_tmp.resize(0,0);

        }else{
            G_tmp = MatrixXd::Zero(ind_n,ind_n);
            for(i=0;i<ve_i;i++) G_tmp.noalias()+=var_tmp[i]*G[i];
            MatrixXd GZ_tmp = G_tmp*Zt;
            G_tmp.resize(0,0);
            V.noalias() = Z*GZ_tmp;
            GZ_tmp.resize(0,0);
        }
#pragma omp parallel for
        for(j=0;j<Y.size();j++) V(j,j)+=var(ve_i);
        double *V_data  = V.data();
        lda  = n = V.cols();
        ipiv = (int*)malloc(n*sizeof(int));
        info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, V_data, lda, ipiv);
        info = LAPACKE_dgetri(LAPACK_COL_MAJOR, n, V_data, lda, ipiv);
        free(ipiv);
        IVX.noalias() =  V * X;          //chunkao wang, 2013-04-01
        XIVX.noalias() =  Xt * IVX;
        lda = n = XIVX.cols();
        MatrixXd XIVX_B = MatrixXd::Identity(n,n);
        double   *s     = (double*)malloc(n*n*sizeof(double));
        double   rcond  = 1e-10;
        int      rank;
        info = LAPACKE_dgelss(LAPACK_COL_MAJOR,n,n,n,XIVX.data(),n,XIVX_B.data(),n,s,rcond,&rank);
        XIVX = XIVX_B;
        XIVX_B.resize(0,0);
        free(s);
        MatrixXd IVXIXIVX_tmp = IVX * XIVX;
        V.noalias() -= IVXIXIVX_tmp * IVX.transpose(); // V used as P
        IVX.resize(0,0);
        IVXIXIVX_tmp.resize(0,0);
//=================IAI matrix begin====================================================
        PY.resize(0);
        ZPY.resize(0);
        Z1PY.resize(0);
        PY  = V*Y;
        ZPY = Zt*PY;
        if(p.modelpe_f) Z1PY = Z1t*PY;
        double   numerator,trace;
        VectorXd dl(var.size());
        for(i=0;i<var.size();i++){
	        PZGZ.resize(0,0);
	        if(i<ve_i){
	            if(var_name[i]!="PE"){
	                VectorXd GZPY_tmp = G[i] * ZPY;
	                numerator = ZPY.dot(GZPY_tmp);
	                GZPY_tmp.resize(0);
	                MatrixXd GZ_tmp =  G[i] * Zt;
	                MatrixXd ZGZ_tmp = Z * GZ_tmp;
	                GZ_tmp.resize(0,0);
	                int size_P = V.cols();
	                MatrixXd PZGZ_tmp = MatrixXd::Zero(size_P,size_P);
	                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,size_P, size_P, size_P, alpha,V.data(), size_P, ZGZ_tmp.data(),size_P,beta, PZGZ_tmp.data(), size_P);
	                trace = PZGZ_tmp.trace();
	                ZGZ_tmp.resize(0,0);
	            }else{
	                VectorXd IZ1PY_tmp = Identy * Z1PY;
	                numerator = Z1PY.dot(IZ1PY_tmp);
	                IZ1PY_tmp.resize(0);
	                MatrixXd IZ1_tmp=Identy*Z1t;
	                MatrixXd Z1IZ1_tmp = Z1*IZ1_tmp;
	                IZ1_tmp.resize(0,0);
	                int size_P = V.cols();
	                MatrixXd PZ1IZ1_tmp = MatrixXd::Zero(size_P,size_P);
	                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,size_P, size_P, size_P, alpha,V.data(), size_P, Z1IZ1_tmp.data(),size_P,beta, PZ1IZ1_tmp.data(), size_P);
	                trace = PZ1IZ1_tmp.trace();
	                Z1IZ1_tmp.resize(0,0);
	            }
	        }
	        else if(i==ve_i){
	            numerator = PY.dot(PY);
	            trace = V.trace();   // V used as P
	        }
	        var[i]=var_tmp[i]*numerator/trace;
	        dl[i]=0.5*(numerator-trace);
        }
        double   dif_max = 0.0;
        VectorXd dif_tmp(var.size());
        update_dif_max(var,var_tmp,dif_tmp,dif_max);
        if((iter>=(p.iter_ai_start-1))||iter==(p.iter_n-1)||dif_max<p.iter_tolerance){
	        vector<VectorXd> GPY;
	        for(i=0;i<ve_i;i++){
	            if(var_name[i]=="PE"){
                    VectorXd Z1PY_tmp = Z1t*PY;
	                VectorXd IZ1PY_tmp = Identy*Z1PY_tmp;
	                Z1PY_tmp.resize(0);
	                VectorXd IPY_tmp = Z1*IZ1PY_tmp;
	                IZ1PY_tmp.resize(0);
	                GPY.push_back(IPY_tmp);
	            }else{
	                VectorXd ZPY_tmp = Zt*PY;
	                VectorXd GZPY_tmp = G[i]*ZPY_tmp;
	                ZPY_tmp.resize(0);
	                VectorXd GPY_tmp = Z*GZPY_tmp;
	                GZPY_tmp.resize(0);
	                GPY.push_back(GPY_tmp);
	            }
	        }
	        GPY.push_back(PY);
	        MatrixXd AI(var.size(),var.size());
	        for(i=0;i<var.size();i++){
	            VectorXd YPGP_i = (GPY[i]).transpose()*V;
	            AI(i,i) = YPGP_i.dot(GPY[i]);
	            for(j=i+1;j<var.size();j++){
	                AI(i,j) = YPGP_i.dot(GPY[j]);
	                AI(j,i) = AI(i,j);
	            }
	            YPGP_i.resize(0);
	        }
	        AI = 0.5 * AI;
	        IAI = AI.inverse();
//	        cout<<"AI\n"<<AI<<endl;
//	        cout<<"IAI\n"<<IAI<<endl;
//	        cout<<"DL\n"<<dl<<endl;
        }
        if(iter>=(p.iter_ai_start-1)){
	        VectorXd pp = IAI * dl;
	        double   alpha = 1.0;
	        var_ai = var_tmp + alpha * pp;
//	        var = var_ai;
	        if(var_ai.minCoeff()>=var_tol){
	            flog<<"AI-REML\n";
	            var = var_ai;
	        }
	        else{
	            flog<<"EM-REML\n";
	        }
	        update_dif_max(var,var_tmp,dif_tmp,dif_max);
        }
        else{
            flog<<"EM-REML\n";
        }
        output_greml_iter_var(var_output,iter,dif_tmp,flog);
        if(iter==p.iter_n-1 || dif_max <p.iter_tolerance){
//            output_greml_inv_ai(var_output,IAI,flog);
        ;
        }
//    cout<<G_tmp.block(0,0,6,6)<<"\n\n";
        if (p.iter_tolerance_her>=0){
            VectorXd  h2se(var.size());
            if (iter==0){
		        get_h2(IAI,h2_tmp,h2se);
		    }
		    else{
		        double   dif_max_her = 0.0;
			    VectorXd dif_tmp_her(var.size());
			    VectorXd  h2(var.size());
			    get_h2(IAI,h2,h2se);
			    update_dif_max(h2,h2_tmp,dif_tmp_her,dif_max_her);
                output_greml_her(h2,dif_tmp_her,flog);
			    if(dif_max_her<p.iter_tolerance_her){
//			        output_greml_inv_ai(var_output,IAI,flog);
                    clock_t time_b_iter;
                    time_b_iter = clock();
                    flog<<"The time (seconds) cost for this iteration is      :\t"<<(int)((double)(time_b_iter-time_a_iter)/CLOCKS_PER_SEC)<<endl;
			        break;
			    }
		        h2_tmp=h2;
		    }
        }
        clock_t time_b_iter;
        time_b_iter = clock();
        flog<<"The time (seconds) cost for this iteration is      :\t"<<(int)((double)(time_b_iter-time_a_iter)/CLOCKS_PER_SEC)<<endl;
        if(dif_max < p.iter_tolerance) break;
    }
    time_b = clock();
    VectorXd  h2(var.size());
    VectorXd  h2se(var.size());
	get_h2(IAI,h2,h2se);
	output_greml_inv_ai(var_output,IAI,flog);
	flog<<"The time (seconds) cost for all iterations is      :\t"<<(int)((double)(time_b-time_a)/CLOCKS_PER_SEC)<<endl;
    output_greml_her_a(var_output,h2,h2se,flog);
    var_output.close();
    fprintf(stderr,"\n======== GREML finished ======\n\n");
}
void Greml_ce::gblup_ce(const parameters &p,SparseMatrix<double> &Xt,SparseMatrix<double> &Zt,int ind_n,int ve_i,ofstream &flog){
    fprintf(stderr,"\n========= GBLUP begins =========\n\n");
    int i,j,m,n;
    clock_t time_a,time_b,time_a_gblup, time_b_gblup;
    time_a_gblup= clock();
    time_a= clock();
    string gblup_file=p.outprefix+"_gblup.csv",fixed_effect_file=p.outprefix+"_fixed_effect.txt";
    VectorXd  var_tmp;
    var_tmp = var;
    int       *ipiv;
    int       lda;
    int       info;
    MatrixXd  V;
    MatrixXd  IVX;
    MatrixXd  XIVX;
    VectorXd  PY;
    VectorXd  ZPY;
    MatrixXd  PZGZ;
    VectorXd  Z1PY;
    MatrixXd  PZ1GZ1;
    MatrixXd G_tmp;
    if(p.modelpe_f){
        MatrixXd Identy;
        MatrixXd Gpe;
        MatrixXd GpeZ1_tmp;
        int z1cols;
        z1cols=Z1.cols();
        Identy=MatrixXd::Identity(z1cols,z1cols);
        G_tmp = MatrixXd::Zero(ind_n,ind_n);
        for(i=0;i<ve_i;i++){
            if(var_name[i]=="PE"){
                Gpe=var_tmp[i]*Identy;
                Identy.resize(0,0);
            }else{
                G_tmp.noalias()+=var_tmp[i]*G[i];
            }
        }
        MatrixXd GZ_tmp = G_tmp*Zt;
        G_tmp.resize(0,0);
        V = Z * GZ_tmp;
        GZ_tmp.resize(0,0);
        GpeZ1_tmp=Gpe*Z1.transpose();
        Gpe.resize(0,0);
        V+=Z1*GpeZ1_tmp;
        GpeZ1_tmp.resize(0,0);

    }else{
        G_tmp = MatrixXd::Zero(ind_n,ind_n);
        for(i=0;i<ve_i;i++) G_tmp.noalias()+=var_tmp[i]*G[i];
        MatrixXd GZ_tmp = G_tmp*Zt;
        G_tmp.resize(0,0);
        V = Z * GZ_tmp;
        GZ_tmp.resize(0,0);
    }
#pragma omp parallel for
    for(j=0;j<Y.size();j++){
      V(j,j) += var[ve_i];
    }
    double *V_data = V.data();
    lda = n = V.cols();
    ipiv = (int*)malloc(n*sizeof(int));
    info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, V_data, lda, ipiv);
    // Returns the determinant of A. 08/08/2014
    double det = 1.0;
    int    neg = 0;
	if(info>0){
	    det = 0.0;
	}
	else{
	    for(int i=0;i<n;i++){
	        det *= V_data[i*n+i];
	        if(ipiv[i]!=(i+1)){
			    neg = !neg;
		    }
		}
	    det = neg?-det:det;
	}
		// end of determinant
    info = LAPACKE_dgetri(LAPACK_COL_MAJOR, n, V_data, lda, ipiv);
    free(ipiv);
    IVX  = V * X;
    XIVX = Xt * IVX;
    lda = n = XIVX.cols();
    MatrixXd XIVX_B = MatrixXd:: Identity(n,n);
    double *s = (double*)malloc(n*n*sizeof(double));
    double rcond = 1e-10;
    int rank;
    info = LAPACKE_dgelss(LAPACK_COL_MAJOR,n,n,n,XIVX.data(),n,XIVX_B.data(),n,s,rcond,&rank);
    XIVX = XIVX_B;
    XIVX_B.resize(0,0);
    free(s);
    V.noalias() -= IVX * XIVX * IVX.transpose();  // V used ad P
    PY  = V * Y;   // V used as P
    ZPY = Zt * PY;
    vector<VectorXd> ebv;
    VectorXd ebv_g;
    if(p.modelpe_f){
        int pecount,ind1;
        VectorXd ebv_tmpe;
        VectorXd ebv_pe;
        ebv_g = VectorXd::Zero(ind_n);
        ebv_pe= VectorXd::Zero(ind_n);
        Z1PY=Z1.transpose()*PY;
        for(i=0;i<ve_i;i++){
            if(var_name[i]=="PE"){
                ebv_tmpe = var_tmp[i]*Z1PY;
            }else{
                VectorXd GZPY_tmp = G[i]*ZPY;
                VectorXd ebv_tmp = var_tmp[i]*GZPY_tmp;
                GZPY_tmp.resize(0);
                ebv_g.noalias()+=ebv_tmp;
                ebv.push_back(ebv_tmp);
            }
        }
        pecount=0;
        IDX *ind_tmp1, *current1;
        HASH_ITER(hh,ind_geno,current1,ind_tmp1){
            ind1 = current1->index;
            IDX *ind_idx_tmp1;
            HASH_FIND_STR(ind_idx, current1->name,ind_idx_tmp1);
            if (ind_idx_tmp1!=NULL){
                ebv_pe(ind1)=ebv_tmpe(pecount);
                pecount+=1;
            }
        }
        ebv_g.noalias()+=ebv_pe;
        ebv.push_back(ebv_pe);
    }else{
        ebv_g = VectorXd::Zero(ind_n);
        for(i=0;i<ve_i;i++){
            VectorXd GZPY_tmp = G[i]*ZPY;
            VectorXd ebv_tmp = var_tmp[i]*GZPY_tmp;
            GZPY_tmp.resize(0);
            ebv_g.noalias()+=ebv_tmp;
            ebv.push_back(ebv_tmp);
        }
    }
    ebv.push_back(ebv_g);
    ebv_g.resize(0);
    if(p.reml_ce_rel_f){
        fprintf(stderr,"\nCalculating reliability!\n");
        vector<VectorXd> rel;
		for(i=0;i<ve_i;i++){
		    if(var_name[i]!="PE"){
		        VectorXd rel_tmp(ind_n);
#pragma omp parallel for
			    for(j=0;j<ind_n;j++){
				    VectorXd ZG_col = Z * G[i].col(j);
				    VectorXd PZG_tmp = V * ZG_col;
				    rel_tmp(j) = ZG_col.dot(PZG_tmp);
				    rel_tmp(j) = var_tmp(i)*rel_tmp(j)/(G[i])(j,j);
			    }
			rel.push_back(rel_tmp);
		    }
		}
		if(var.size()>2){
			VectorXd rel_tmp(ind_n);
#pragma omp parallel for private(j)
			for(int j=0;j<ind_n;j++){
				VectorXd G_col = VectorXd::Zero(ind_n);
				double   G_diag = 0.0;
				for(int i=0;i<ve_i;i++){
				    if(var_name[i]!="PE"){
					    G_col.noalias() += var_tmp(i)*G[i].col(j);
					    G_diag += var_tmp(i) * (G[i])(j,j);
					}
				}
				VectorXd ZG_col  = Z * G_col;
				G_col.resize(0);
				VectorXd PZG_tmp = V * ZG_col;
				rel_tmp(j) = ZG_col.dot(PZG_tmp);
				rel_tmp(j) /= G_diag;
			}
			rel.push_back(rel_tmp);
		}
		else if (var.size()==2){
			rel.push_back(rel[0]);
		}

		ofstream ebv_output(gblup_file.c_str());
        fprintf(stderr,"Printing GBLUP results from GREML_CE.\n");
		ebv_output<<"ID,";
		for(i=0;i<ve_i;i++){
		    if(var_name[i]!="PE"){
		        ebv_output<<"GBLUP_"+var_name[i]<<","<< "Reliability_" + var_name[i]<<",";
		    }else{
		        ebv_output<<"GBLUP_"+var_name[i]<<",";
		    }
		}
        ebv_output<<"GBLUP_G,Reliability_G,Train./Valid.\n";
        IDX *ind_tmp, *current;
        HASH_ITER(hh,ind_geno,current,ind_tmp){
            ebv_output<<current->name<<",";
            int ind = current->index;
            for(i=0;i<var.size();i++) ebv_output<<setprecision(6)<<(ebv[i])(ind)<<","<<(rel[i])(ind)<<",";
            IDX *ind_idx_tmp;
            HASH_FIND_STR(ind_idx, current->name,ind_idx_tmp);
            ind_idx_tmp==NULL?ebv_output<<"V"<<endl:ebv_output<<"T"<<endl;
        }
        ebv_output.close();
        for(i=0;i<ebv.size();i++) ebv[i].resize(0);
    }
    else{
        ofstream ebv_output(gblup_file.c_str());
        fprintf(stderr,"Output of GBLUP results from GREML_CE.\n");
		ebv_output<<"ID,";
        for(i=0;i<ve_i;i++) ebv_output<<"GBLUP_"+var_name[i]<<",";
        ebv_output<<"GBLUP_G,Train./Valid.\n";
        IDX *ind_tmp, *current;
        HASH_ITER(hh,ind_geno,current,ind_tmp){
            ebv_output<<current->name<<",";
            int ind = current->index;
            for(i=0;i<var.size();i++) ebv_output<<setprecision(6)<<(ebv[i])(ind)<<",";
            IDX *ind_idx_tmp;
            HASH_FIND_STR(ind_idx, current->name,ind_idx_tmp);
            ind_idx_tmp==NULL?ebv_output<<"V"<<endl:ebv_output<<"T"<<endl;
        }
        ebv_output.close();
        for(i=0;i<ebv.size();i++) ebv[i].resize(0);
    }
    VectorXd y_tmp = IVX.transpose() * Y;
	VectorXd b     = XIVX * y_tmp;
	y_tmp.resize(0);
	ofstream fixed_effect_output(fixed_effect_file);
	fixed_effect_output << setw(20) << right << "Fixed_effect";
	fixed_effect_output << setw(STR_LEN) << right << "Level_name";
	fixed_effect_output << setw(15) << right << "Level";
	fixed_effect_output << setw(15) << right << "Value" << endl;
	int count = 0;
	int effect_count = 0;
	fixed_effect_output << setw(20) << right << effect_count;
	fixed_effect_output << setw(STR_LEN) << right<<"mu";
	fixed_effect_output << setw(15)<<right<<1;
	fixed_effect_output << setw(15) << right<< setprecision(6) << scientific<<b(0)<<endl;
	count++;
	effect_count++;
	for(i=0;i<cross.size();i++){
	    IDX *cross_tmp,*current;
	    HASH_ITER(hh,cross_idx[i],current,cross_tmp){
		fixed_effect_output << setw(20) << right << effect_count;
		fixed_effect_output << setw(STR_LEN) << right << current->name;
		fixed_effect_output << setw(15) << right<< fixed << setprecision(0) << current->index;
		fixed_effect_output << setw(15) << right << setprecision(6) << scientific << b(count) << endl;
		count++;
		}
	    effect_count++;
	}
	for(i=0;i<covar.size();i++){
	    fixed_effect_output << setw(20) << right << effect_count;
		fixed_effect_output << setw(STR_LEN) << right << "Covariable";
		fixed_effect_output << setw(15) << right << fixed << setprecision(0) << 1;
		fixed_effect_output << setw(15) << right<< setprecision(6) << scientific<< b(count) << endl;
		count++;
		effect_count++;
	}
	VectorXd Xb     = X * b;
	vector<string> current_name;
	fixed_effect_output<<"Xb"<<endl;
	IDX *ind_tmp, *current;
	HASH_ITER(hh,ind_geno,current,ind_tmp){
        int ind = current->index;
        IDX *ind_idx_tmp;
        HASH_FIND_STR(ind_idx, current->name,ind_idx_tmp);
        if (ind_idx_tmp!=NULL) current_name.push_back(current->name);
    }
    for(i=0;i<Xb.size();i++) fixed_effect_output<<current_name[i]<<"\t"<<setprecision(6)<<scientific<<Xb(i)<<endl;
	fixed_effect_output.close();
	current_name.clear();
	Xb.resize(0);
	V.resize(0,0);
    IVX.resize(0,0);
    XIVX.resize(0,0);
    PY.resize(0);
    PZGZ.resize(0,0);
    b.resize(0);
    for(i=0;i<G.size();i++) G[i].resize(0,0);
    time_b= clock();
    flog<<"\nThe time (seconds) cost for GBLUP estimation is    :\t"<<(int)((double)(time_b-time_a)/CLOCKS_PER_SEC)<<endl;
	fprintf(stderr,"\n========= GBLUP finished =========\n\n");
	if(p.output_mrk_effect){
	    fprintf(stderr,"Calculating marker effects and heritabilities.\n");
	    get_SNP_eff(p,ZPY);
	}
}
void Greml_ce::reml(const parameters &p,ofstream &flog){
    SparseMatrix<double> Xt = X.transpose();
    SparseMatrix<double> Zt = Z.transpose();
    ind_n=HASH_COUNT(ind_geno);
    int ve_i=var.size()-1;
    double tolerance=p.iter_tolerance+1.0;
    flog<<"The number of individuals in genotypes file is     :\t"<<ind_n<<endl;
	flog<<"The number of individuals in training dataset is   :\t"<<HASH_COUNT(ind_idx)<<endl;
	flog<<"The number of individuals in validation dataset is :\t"<<ind_n-HASH_COUNT(ind_idx)<<endl;
    if(p.cin_var_f){
        gblup_ce(p,Xt,Zt,ind_n,ve_i,flog);
    }
    else{
        variance_component_calculation(p,Xt,Zt,ind_n,ve_i,flog);
        gblup_ce(p,Xt,Zt,ind_n,ve_i,flog);
    }
}