#include "read_geno.h"
Greml::~Greml() {}
void Greml::read_phen_pre(const parameters &p,ofstream &flog){
    int i,j;
    ifstream input_phen(p.phenotype.c_str());
    if(!input_phen){
      fprintf(stderr,"ERROR: Can not open the file %s\n",p.phenotype.c_str());
      throw;
    }
    int cross_n=cross.size();
    cross_idx = (IDX**)malloc(cross_n*sizeof(IDX*));
    for(i=0;i<cross_n;i++) cross_idx[i] = NULL;
    ind_idx = NULL;
    idx_ind = NULL;
    int pos_max = 0;
    if(p.trait_pos>pos_max) pos_max=p.trait_pos;
    for(i=0;i<cross.size();i++){
      if(cross[i]>pos_max) pos_max = cross[i];
    }
    for(i=0;i<covar.size();i++){
      if(covar[i]>pos_max) pos_max=covar[i];
    }
    string line;
    getline(input_phen,line);
    line.clear();
    int rec_n = 0;
    while(getline(input_phen,line)){
      vector<string> words = str_split(" \t",line);
      if(words.size()==0) continue;
      if(words.size()<pos_max){
	    fprintf(stderr,"ERROR: Number of fields  in line %d in phenotype file is too few!\n",rec_n +1);
	    throw;
      }
      double y_tmp = atof(words[p.trait_pos-1].c_str());
      if(y_tmp == p.missing_phen_val) continue;
      for(i=0;i<cross.size();i++){
	    IDX *cross_tmp;
	    char name[STR_LEN];
	    strcpy(name,words[cross[i]-1].c_str());
	    if(atof(name)!=p.missing_phen_val){
	    HASH_FIND_STR(cross_idx[i],name,cross_tmp);
	        if(cross_tmp==NULL){
	            cross_tmp = (IDX*)malloc(sizeof(IDX));
	            strcpy(cross_tmp->name,name);
	            cross_tmp->index = HASH_COUNT(cross_idx[i]);
	            HASH_ADD_STR(cross_idx[i],name,cross_tmp);
	        }
	    }
    }
      char name[STR_LEN];
      strcpy(name,words[0].c_str());
      IDX *ind_tmp;
      HASH_FIND_STR(ind_idx,name,ind_tmp);
      if(ind_tmp==NULL){
        ind_tmp = (IDX*)malloc(sizeof(IDX));
	    strcpy(ind_tmp->name,name);
	    ind_tmp->index = HASH_COUNT(ind_idx);
	    HASH_ADD_STR(ind_idx,name,ind_tmp);
      }
      IDX *idx_tmp;
      idx_tmp = (IDX*)malloc(sizeof(IDX));
      strcpy(idx_tmp->name,name);
      idx_tmp->index = rec_n;
      HASH_ADD_INT(idx_ind,index,idx_tmp);
      words.clear();
      line.clear();
      rec_n++;
    }
    input_phen.close();
    phen_n = rec_n;
    flog<<endl;
    for(i=0;i<cross.size();i++) flog<<"The number of levels for descrete variable "<<i+1<<" is    :\t"<<HASH_COUNT(cross_idx[i])<<endl;
}
void Greml::read_phen(const parameters &p){
    int i,j;
    EFF eff_tmp;
    int count = 0;
    eff_tmp.start = count;
    eff_tmp.len = 1;
    eff_fixed.push_back(eff_tmp);
    count += eff_tmp.len;
    string ok="ok";
    for(i=0;i<cross.size();i++){
      eff_tmp.start = count;
      //eff_tmp.len = HASH_COUNT(cross_idx[i])-1;
      eff_tmp.len = HASH_COUNT(cross_idx[i]);
      eff_fixed.push_back(eff_tmp);
      count += eff_tmp.len;
    }
    for(i=0;i<covar.size();i++){
      eff_tmp.start = count;
      eff_tmp.len = 1;
      eff_fixed.push_back(eff_tmp);
      count += eff_tmp.len;
    }
    Y.resize(phen_n);
    Y_index.resize(phen_n);
    X.resize(phen_n,count);
    count = 0;
    eff_tmp.start = count;
    eff_tmp.len = HASH_COUNT(ind_geno);
    eff_rand.push_back(eff_tmp);
    count += eff_tmp.len;
    Z.resize(phen_n,count);
    ifstream input_phen(p.phenotype.c_str());
    if(!input_phen){
      fprintf(stderr,"ERROR: Can not open the file %s!\n",p.phenotype.c_str());
      throw;
    }
    string line;
    getline(input_phen,line);
    line.clear();
    int rec_n = 0;
    while(getline(input_phen,line)){
        vector<string> words = str_split(" \t", line);
        if(words.size()==0) continue;
        double y_tmp  = atof(words[p.trait_pos-1].c_str());
        if(y_tmp == p.missing_phen_val) continue;
        Y(rec_n) = y_tmp; count = 0; X.coeffRef(rec_n,count) += 1.0; count++;
        for(i=0;i<cross.size();i++){
            IDX *cross_tmp;char name[STR_LEN];
            strcpy(name,words[cross[i]-1].c_str());
            if(atof(name)!=p.missing_phen_val){
                HASH_FIND_STR(cross_idx[i],name,cross_tmp);int pos = cross_tmp->index;
                //if(pos < HASH_COUNT(cross_idx[i])-1)
                if(pos < HASH_COUNT(cross_idx[i])){
                    int col = eff_fixed[count].start + pos;
                    X.coeffRef(rec_n,col) += 1.0;
                }
                count++;
            }
        }
        for(i=0;i<covar.size();i++){
            int col = eff_fixed[count].start;
            double covar_tmp = atof(words[covar[i]-1].c_str());
            if(covar_tmp == p.missing_phen_val) covar_tmp = 0.0;
            X.coeffRef(rec_n,col) += covar_tmp; count++;
        }
        IDX *ind_tmp; char name[STR_LEN];strcpy(name,words[0].c_str());
        HASH_FIND_STR(ind_geno,name,ind_tmp);
        if(ind_tmp==NULL){
            cerr<<"ERROR: Individual "<<name<<" is not genotyped!"<<endl;
            throw;
        }
        int pos = ind_tmp->index; Z.coeffRef(rec_n,pos) += 1.0;
        Y_index(rec_n)=pos; my_yindex[pos]=ok; rec_n++;
    }
    input_phen.close();
    int zrows,zcols,z1cols,sumzcol,z1colcount;
    zrows=Z.rows();
    zcols=Z.cols();
    cerr<<"***The matrix size of Z ******"<<endl;
    cerr<<"Z rows: "<<zrows<<"   Z cols: "<<zcols<<endl;
    cerr<<"***The matrix size of X ******"<<endl;
    cerr<<"X rows: "<<X.rows()<<"   X cols: "<<X.cols()<<endl;
    cerr<<"***The matrix size of Y ******"<<endl;
    cerr<<"Y rows: "<<Y.rows()<<"   Y cols: "<<Y.cols()<<endl;
 }