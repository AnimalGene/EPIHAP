#include "tools.h"
vector<string> str_split(const string &delimiter,const string &str){
    vector<string> words;
    int str_len = str.length();
    int del_len = delimiter.length();
    if (str_len==0 || del_len==0) return words;
    char *line = new char [str_len+1];
    strcpy(line, str.c_str());
    char *sep  = new char [del_len+1];
    strcpy(sep, delimiter.c_str());
    string tmp_str;
    char   *token = strtok(line,sep);
    while(token!=NULL){
      tmp_str = token;
      words.push_back(tmp_str);
      tmp_str.clear();
      token = strtok(NULL,sep);
    }
    delete [] line;
    delete [] sep;
    delete token;
    return words;
  }
void word2geno(vector<string> &words,int row,int col_start,MatrixXi &geno,MatrixXd &geno_stat,int round_geno){
int genoValue;
#pragma omp parallel for
    for(int j=1;j<words.size();j++){
      int col = j - 1+col_start;
      genoValue=atoi(words[j].c_str());
      geno(row,col)=genoValue;
      double geno_tmp = (double)genoValue;
      if(round_geno==1){
        geno_tmp = floor(geno_tmp+0.5);
      }
      if(geno_tmp==0.0){
        geno_stat(0,col) += 1.0;
        geno_stat(1,col) += geno_tmp;
        geno_stat(9,col) += 1.0;
      }
      else if(geno_tmp==1.0){
        geno_stat(0,col) += 1.0;
        geno_stat(1,col) += geno_tmp;
        geno_stat(8,col) += 1.0;
      }
      else if(geno_tmp==2.0){
        geno_stat(0,col) += 1.0;
        geno_stat(1,col) += geno_tmp;
        geno_stat(10,col) += 1.0;
      }
      else{
        geno_stat(11,col) += 1.0;
      }
    }
}
void get_ind_list_geno(string word,int row,IDX **ind_geno){
    char name[STR_LEN];
    strcpy(name,word.c_str());
    IDX *ind_geno_tmp;
    HASH_FIND_STR(*ind_geno,name,ind_geno_tmp);
    if(ind_geno_tmp==NULL){
      ind_geno_tmp=(IDX*)malloc(sizeof(IDX));
      strcpy(ind_geno_tmp->name,name);
      ind_geno_tmp->index=row;
      HASH_ADD_STR(*ind_geno,name,ind_geno_tmp);
    }
    else{
      fprintf(stderr,"ERROR: Individual %s has been in genotype file!\n",name);
      throw;
    }
  }
void genostat_update(MatrixXd &geno_stat){
    int mrk_num = geno_stat.cols();
#pragma omp parallel for
    for(int j=0;j<mrk_num;j++){
        double freq=geno_stat(1,j)/(2*geno_stat(0,j));
        geno_stat(1,j) = freq;
        geno_stat(2,j) = 2*freq;
        geno_stat(3,j) = 2*freq-1;
        geno_stat(4,j) = 2*freq-2;
        geno_stat(5,j) = -2*freq*freq;
        geno_stat(6,j) = 2*freq*(1-freq);
        geno_stat(7,j) = -2*(1-freq)*(1-freq);
        geno_stat(9,j) /= geno_stat(0,j);
        geno_stat(8,j) /= geno_stat(0,j);
        geno_stat(10,j) /= geno_stat(0,j);
    }
 }
void genostat_update_epi(MatrixXd &geno_stat,MatrixXd &geno_stat_epi, int snp1_index,int mrk_c){
    int j2;
#pragma omp parallel for
    for(j2=snp1_index+1;j2<mrk_c;j2++){
//        A by A,A by D, D by A, D by D;
        VectorXd SNP1  = VectorXd::Zero(36);
        VectorXd SNP2  = VectorXd::Zero(36);
        VectorXd SNP12 = VectorXd::Zero(36);

        SNP1<<geno_stat(2,snp1_index),geno_stat(2,snp1_index),geno_stat(2,snp1_index),geno_stat(3,snp1_index),geno_stat(3,snp1_index),geno_stat(3,snp1_index),geno_stat(4,snp1_index),geno_stat(4,snp1_index),geno_stat(4,snp1_index),\
        geno_stat(2,snp1_index),geno_stat(2,snp1_index),geno_stat(2,snp1_index),geno_stat(3,snp1_index),geno_stat(3,snp1_index),geno_stat(3,snp1_index),geno_stat(4,snp1_index),geno_stat(4,snp1_index),geno_stat(4,snp1_index),\
        geno_stat(5,snp1_index),geno_stat(5,snp1_index),geno_stat(5,snp1_index),geno_stat(6,snp1_index),geno_stat(6,snp1_index),geno_stat(6,snp1_index),geno_stat(7,snp1_index),geno_stat(7,snp1_index),geno_stat(7,snp1_index),\
        geno_stat(5,snp1_index),geno_stat(5,snp1_index),geno_stat(5,snp1_index),geno_stat(6,snp1_index),geno_stat(6,snp1_index),geno_stat(6,snp1_index),geno_stat(7,snp1_index),geno_stat(7,snp1_index),geno_stat(7,snp1_index);

        SNP2<<geno_stat(2,j2),geno_stat(3,j2),geno_stat(4,j2),geno_stat(2,j2),geno_stat(3,j2),geno_stat(4,j2),geno_stat(2,j2),geno_stat(3,j2),geno_stat(4,j2),\
        geno_stat(5,j2),geno_stat(6,j2),geno_stat(7,j2),geno_stat(5,j2),geno_stat(6,j2),geno_stat(7,j2),geno_stat(5,j2),geno_stat(6,j2),geno_stat(7,j2),\
        geno_stat(2,j2),geno_stat(3,j2),geno_stat(4,j2),geno_stat(2,j2),geno_stat(3,j2),geno_stat(4,j2),geno_stat(2,j2),geno_stat(3,j2),geno_stat(4,j2),\
        geno_stat(5,j2),geno_stat(6,j2),geno_stat(7,j2),geno_stat(5,j2),geno_stat(6,j2),geno_stat(7,j2),geno_stat(5,j2),geno_stat(6,j2),geno_stat(7,j2);

        SNP12=SNP1.cwiseProduct(SNP2);
        geno_stat_epi.row(j2-snp1_index-1)=SNP12.transpose();
    }
 }
void get_het_sum(MatrixXd &geno_stat,int is_hwd,double *hetero_total, double *hetero_total2){
    double het0=0.0;
    double het1=0.0;
    int mrk_num = geno_stat.cols();
    if(is_hwd==0){
#pragma omp parallel for reduction(+:het0,het1)
      for(int j=0;j<mrk_num;j++){
	    het0 += geno_stat(6,j);
	    het1 += geno_stat(6,j)*geno_stat(6,j);
      }
    }
    else if(is_hwd==3){
#pragma omp parallel for reduction(+:het0,het1)
      for(int j=0;j<mrk_num;j++){
	    het0 += geno_stat(6,j);
	    het1 += geno_stat(6,j)*(1.0-geno_stat(6,j));
      }
    }
    *hetero_total  += het0;
    *hetero_total2 += het1;
}
void geno_transform(const MatrixXi& geno_org,const MatrixXd& geno_stat,string geno_type,int is_hwd,MatrixXd& geno_new){
    int i,j;
    int row_n = geno_org.rows();
    int col_n = geno_org.cols();
    if(is_hwd==0){
      // hwe genotype transform
      if(geno_type.compare("A")==0){
#pragma omp parallel for schedule(dynamic,1) collapse(2)
        for(j=0;j<col_n;j++){
          for(i=0;i<row_n;i++){
            if(geno_org(i,j)==0)      geno_new(i,j) = geno_stat(2,j);
            else if(geno_org(i,j)==1) geno_new(i,j) = geno_stat(3,j);
            else if(geno_org(i,j)==2) geno_new(i,j) = geno_stat(4,j);
            else                        geno_new(i,j) = 0.0;
          }
        }
      }
      if(geno_type.compare("D")==0){
#pragma omp parallel for schedule(dynamic,1) collapse(2)
        for(j=0;j<col_n;j++){
          for(i=0;i<row_n;i++){
            if(geno_org(i,j)==0)      geno_new(i,j) = geno_stat(5,j);
            else if(geno_org(i,j)==1) geno_new(i,j) = geno_stat(6,j);
            else if(geno_org(i,j)==2) geno_new(i,j) = geno_stat(7,j);
            else                        geno_new(i,j) = 0.0;
          }
        }
      }
    }
    else if(is_hwd==1){
      if(geno_type.compare("A")==0){
        for(j=0;j<col_n;j++){
          double q     = geno_stat(1,j);
          double p     = 1-geno_stat(1,j);
          double denom = p*q-0.25*geno_stat(8,j);
          double w0    = p*q*q/denom;
          double w1    = 0.5*p*q*(q-p)/denom;
          double w2    = -1.0*p*q*p/denom;
#pragma omp parallel for private(i)
          for(i=0;i<row_n;i++){
            if(geno_org(i,j)==0.0)      geno_new(i,j) = w0;
            else if(geno_org(i,j)==1.0) geno_new(i,j) = w1;
            else if(geno_org(i,j)==2.0) geno_new(i,j) = w2;
            else                        geno_new(i,j) = 0.0;
          }
        }
      }
      else if(geno_type.compare("D")==0){
        for(j=0;j<col_n;j++){
          double q     = geno_stat(1,j);
          double p     = 1-geno_stat(1,j);
          double qg    = geno_stat(10,j);
          double rg    = geno_stat(8,j);
          double pg    = geno_stat(9,j);
          double denom = p*q-0.25*geno_stat(8,j);
          double w0    = -0.5*qg*rg/denom;
          double w1    = qg*pg/denom;
          double w2    = -0.5*pg*rg/denom;
#pragma omp parallel for private(i)
          for(i=0;i<row_n;i++){
            if(geno_org(i,j)==0)      geno_new(i,j) = w0;
            else if(geno_org(i,j)==1) geno_new(i,j) = w1;
            else if(geno_org(i,j)==2) geno_new(i,j) = w2;
            else                        geno_new(i,j) = 0.0;
          }
        }
      }
    }
    else if(is_hwd==2){
      if(geno_type.compare("A")==0){
        for(j=0;j<col_n;j++){
          double q     = geno_stat(1,j);
          double p     = 1-geno_stat(1,j);
          //double denom = p*q-0.25*geno_stat(8,j);
          double w0    = p*q*q;         // /denom;
          double w1    = 0.5*p*q*(q-p); // /denom;
          double w2    = -1.0*p*q*p;    // /denom;
#pragma omp parallel for private(i)
          for(i=0;i<row_n;i++){
            if(geno_org(i,j)==0)      geno_new(i,j) = w0;
            else if(geno_org(i,j)==1) geno_new(i,j) = w1;
            else if(geno_org(i,j)==2) geno_new(i,j) = w2;
            else                        geno_new(i,j) = 0.0;
          }
          double qg    = geno_stat(10,j);
          double rg    = geno_stat(8,j);
          double pg    = geno_stat(9,j);
          double ss    = pg*q*q+0.25*rg*(q-p)*(q-p)+qg*p*p;
          double sd    = p*q*sqrt(ss);
          geno_new.col(j) /= sd;
        }
      }
      else if(geno_type.compare("D")==0){
        for(j=0;j<col_n;j++){
          double q     = geno_stat(1,j);
          double p     = 1-geno_stat(1,j);
          double qg    = geno_stat(10,j);
          double rg    = geno_stat(8,j);
          double pg    = geno_stat(9,j);
          double w0    = -0.5*qg*rg; // /denom;
          double w1    = qg*pg;      // /denom;
          double w2    = -0.5*pg*rg; // /denom;
#pragma omp parallel for private(i)
          for(i=0;i<row_n;i++){
            if(geno_org(i,j)==0)      geno_new(i,j) = w0;
            else if(geno_org(i,j)==1) geno_new(i,j) = w1;
            else if(geno_org(i,j)==2) geno_new(i,j) = w2;
            else                        geno_new(i,j) = 0.0;
          }
          double ss = 0.0;
          if(pg==0.0 || rg==0.0 || qg==0.0){
            ss  = (4*p*p*q*q)*(pg*qg+0.25*rg*(pg+qg));
          }
          else{
            ss  = pg*rg*qg*(pg*qg+0.25*rg*(pg+qg));
          }
          double sd    = sqrt(ss);
          geno_new.col(j) /= sd;
        }
      }
    }
    else if(is_hwd==3){
      // Su Guosheng's methods
      if(geno_type.compare("A")==0){
#pragma omp parallel for schedule(dynamic,1) collapse(2)
        for(j=0;j<col_n;j++){
          for(i=0;i<row_n;i++){
            if(geno_org(i,j)==0)      geno_new(i,j) = geno_stat(2,j);
            else if(geno_org(i,j)==1) geno_new(i,j) = geno_stat(3,j);
            else if(geno_org(i,j)==2) geno_new(i,j) = geno_stat(4,j);
            else                        geno_new(i,j) = 0.0;
          }
        }
      }
      if(geno_type.compare("D")==0){
#pragma omp parallel for schedule(dynamic,1) collapse(2)
        for(j=0;j<col_n;j++){
          for(i=0;i<row_n;i++){
            if(geno_org(i,j)==0)      geno_new(i,j) = -1.0*geno_stat(6,j);
            else if(geno_org(i,j)==1) geno_new(i,j) = 1-geno_stat(6,j);
            else if(geno_org(i,j)==2) geno_new(i,j) = -1.0*geno_stat(6,j);
            else                        geno_new(i,j) = 0.0;
          }
        }
      }
    }
    else if(is_hwd==4){
      // hwd genotype transform and normalize for each marker
      if(geno_type.compare("A")==0){
        for(j=0;j<col_n;j++){
          double q     = geno_stat(1,j);
          double p     = 1-geno_stat(1,j);
          double denom = p*q-0.25*geno_stat(8,j);
          double w0    = p*q*q/denom;
          double w1    = 0.5*p*q*(q-p)/denom;
          double w2    = -1.0*p*q*p/denom;
#pragma omp parallel for private(i)
          for(i=0;i<row_n;i++){
            if(geno_org(i,j)==0)      geno_new(i,j) = w0;
            else if(geno_org(i,j)==1) geno_new(i,j) = w1;
            else if(geno_org(i,j)==2) geno_new(i,j) = w2;
            else                        geno_new(i,j) = 0.0;
          }
        }
      }
      else if(geno_type.compare("D")==0){
        for(j=0;j<col_n;j++){
          double q     = geno_stat(1,j);
          double p     = 1-geno_stat(1,j);
          double qg    = geno_stat(10,j);
          double rg    = geno_stat(8,j);
          double pg    = geno_stat(9,j);
          double denom = p*q-0.25*geno_stat(8,j);
          double w0    = -0.5*qg*rg /denom;
          double w1    = qg*pg /denom;
          double w2    = -0.5*pg*rg /denom;
#pragma omp parallel for private(i)
          for(i=0;i<row_n;i++){
            if(geno_org(i,j)==0)      geno_new(i,j) = w0;
            else if(geno_org(i,j)==1) geno_new(i,j) = w1;
            else if(geno_org(i,j)==2) geno_new(i,j) = w2;
            else                        geno_new(i,j) = 0.0;
          }
        }
      }
    }
    else if(is_hwd==5){
      // Peter Visshier's methods'
      if(geno_type.compare("A")==0){
#pragma omp parallel for schedule(dynamic,1) collapse(2)
        for(j=0;j<col_n;j++){
          for(i=0;i<row_n;i++){
            if(geno_org(i,j)==0)      geno_new(i,j) = geno_stat(2,j);
            else if(geno_org(i,j)==1) geno_new(i,j) = geno_stat(3,j);
            else if(geno_org(i,j)==2) geno_new(i,j) = geno_stat(4,j);
            else                        geno_new(i,j) = 0.0;
          }
        }
#pragma omp parallel for
        for(j=0;j<col_n;j++){
          if(geno_stat(6,j)<=0.0){
            geno_new.col(j) *=0.0;
          }
          else{
            geno_new.col(j) /= sqrt(geno_stat(6,j));
          }
	}
      }
      if(geno_type.compare("D")==0){
#pragma omp parallel for schedule(dynamic,1) collapse(2)
        for(j=0;j<col_n;j++){
          for(i=0;i<row_n;i++){
            if(geno_org(i,j)==0)      geno_new(i,j) = geno_stat(5,j);
            else if(geno_org(i,j)==1) geno_new(i,j) = geno_stat(6,j);
            else if(geno_org(i,j)==2) geno_new(i,j) = geno_stat(7,j);
            else                        geno_new(i,j) = 0.0;
          }
        }
#pragma omp parallel for
        for(j=0;j<col_n;j++){
          if(geno_stat(6,j)<=0.0){
            geno_new.col(j) *=0.0;
          }
          else{
            geno_new.col(j) /= sqrt((geno_stat(6,j)*geno_stat(6,j)));
          }
        }
      }
    }
    else if(is_hwd==6){
      // Peter Visshier's methods with observed genotype variance
      if(geno_type.compare("A")==0){
        VectorXi size = VectorXi::Constant(col_n,row_n);
#pragma omp parallel for schedule(dynamic,1) collapse(2)
        for(j=0;j<col_n;j++){
          for(i=0;i<row_n;i++){
            if(geno_org(i,j)==0)      geno_new(i,j) = geno_stat(2,j);
            else if(geno_org(i,j)==1) geno_new(i,j) = geno_stat(3,j);
            else if(geno_org(i,j)==2) geno_new(i,j) = geno_stat(4,j);
            else{
              geno_new(i,j) = 0.0;
              size(j) -= 1;
            }
          }
        }
#pragma omp parallel for
        for(j=0;j<col_n;j++){
          double norm = 0.0;
          for(i=0;i<row_n;i++){
            norm += geno_new(i,j)*geno_new(i,j);
          }
          if(norm==0.0){
            fprintf(stderr,"ERROR: tool_geno::geno_transform::is_hwd=6::A\n");
            throw;
          }
	  geno_new.col(j) /= sqrt(norm/size(j));
	}
      }
      if(geno_type.compare("D")==0){
        VectorXi size = VectorXi::Constant(col_n,row_n);
#pragma omp parallel for schedule(dynamic,1) collapse(2)
        for(j=0;j<col_n;j++){
          for(i=0;i<row_n;i++){
            if(geno_org(i,j)==0)      geno_new(i,j) = geno_stat(5,j);
            else if(geno_org(i,j)==1) geno_new(i,j) = geno_stat(6,j);
            else if(geno_org(i,j)==2) geno_new(i,j) = geno_stat(7,j);
            else{
              geno_new(i,j) = 0.0;
              size(j) -= 1;
            }
          }
        }
#pragma omp parallel for
        for(j=0;j<col_n;j++){
          double norm = geno_new.col(j).norm();
          if(norm==0.0){
            fprintf(stderr,"ERRPR: tool_geno::geno_transform::is_hwd=6::D\n");
            throw;
          }
          geno_new.col(j) /= (norm/sqrt(size(j)*1.0));
        }
      }
    }
  }
void save_ind_geno(string filename,IDX **ind_geno){
    ofstream indgeno(filename,ios::out|ios::binary|ios::trunc);
	int hashcnt = HASH_COUNT(*ind_geno);
	indgeno.write((char*)&hashcnt, sizeof(int));
	char name[STR_LEN];
	int indx;
	IDX *ind_tmp, *current;
	HASH_ITER(hh, *ind_geno, current, ind_tmp){
	    strcpy(name,current->name);
	    indx = current->index;
	    indgeno.write((char*)&name, sizeof(name));
		indgeno.write((char*)&indx, sizeof(int));
    }
	indgeno.close();
  }
void load_ind_geno(string filename,IDX **ind_geno){
    ifstream indgeno(filename,ios::in | ios::binary);
    int hashcnt = HASH_COUNT(*ind_geno);
	indgeno.read((char*)&hashcnt, sizeof(int));
	char name[STR_LEN];
	int indx;
	IDX *ind_tmp, *current;
	for (int i=0;i<hashcnt;i++){
	    indgeno.read((char*)&name, sizeof(name));
	    indgeno.read((char*)&indx, sizeof(int));
		get_ind_list_geno(name,indx,ind_geno);
		}
	indgeno.close();
  }
int binomialCoefficients(int n, int k) {
     if (k == 0||k == n) return 1;
     return binomialCoefficients(n-1,k-1)+binomialCoefficients(n-1,k);
}
int calc_sample_size(const string &filename){
    int num;
    num=0;
    ifstream input(filename);
    string line;
    getline(input,line);
    vector<string>words=str_split(" \r\t\n\x0b",line);
    while(getline(input,line).good()){
        num++;
    }
    input.close();
    return num;
}
void read_plain_txt(const char* filename, MatrixXd& matrix){
    vector<double> matrixEntries;
    ifstream matrixDataFile(filename);
    string matrixRowString, matrixEntry;
    int matrixRowNumber = 0;
	while (getline(matrixDataFile, matrixRowString)){
        stringstream matrixRowStringStream(matrixRowString);
        while (getline(matrixRowStringStream, matrixEntry,'\t')){
            matrixEntries.push_back(stod(matrixEntry));
        }
        matrixRowNumber++;
    }
   matrix=Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(matrixEntries.data(), matrixRowNumber, matrixEntries.size()/matrixRowNumber);
   matrixDataFile.close();
}
void saveData(string fileName, MatrixXd &matrix){
    int i,j;
    ofstream fileo(fileName);
    for(i=0;i<matrix.rows();i++){
        for(j=0;j<matrix.cols();j++){
            if(j==(matrix.cols()-1)){
                fileo<<matrix(i,j);
            }
            else{
                fileo<<matrix(i,j)<<"\t";
            }
        }
        fileo <<"\n";
    }
    fileo.close();
}