#include "read_geno.h"
void Greml_ce::read_map(const string &mapFile){
    int index;
    index=-1;
    ifstream biminput(mapFile.c_str());
    string line;
    while(getline(biminput,line).good()){
        vector<string> words = str_split(" \t",line);
        Snp_info one_row;
        one_row.chr=words[0].c_str();
        one_row.pos=words[2].c_str();
        one_row.mrk=words[1].c_str();
        if(index>-1){
            mysnpmap[index]=one_row;
            mychrindexmap[one_row.chr].push_back(index);
            if(mychrindexmap[one_row.chr].size()>2){
                mychrindexmap[one_row.chr][1]=mychrindexmap[one_row.chr][2];
                mychrindexmap[one_row.chr].pop_back();
            }
        }
        index=index+1;
        mrk_n++;
        words.clear();
        line.clear();
    }
    biminput.close();
}
void Greml_ce::read_rgf(const string &rawgFile,MatrixXi &geno,MatrixXd &geno_stat,int round_geno){
    int flag;
    ifstream input_geno(rawgFile);
    if(!input_geno){
	    fprintf(stderr,"\nERROR: Can not open the file %s\n",rawgFile.c_str());
	    throw;
    }
    else{
	    fprintf(stderr,"\nReading genotype file %s ...\n",rawgFile.c_str());
    }
    string   line;
    getline(input_geno,line);
    vector<string> words = str_split(" \r\t\n\x0b", line);
    words.erase(words.begin());
    line.clear();
    flag=0;
    while(getline(input_geno,line).good()){
        vector<string>words=str_split(" \r\t\n\x0b",line);
        if(words.size()==0) continue;
        if(ind_geno==NULL){
            get_ind_list_geno(words[0],flag,&ind_geno);
        }
        else if(HASH_COUNT(ind_geno)!=ind_n){
            get_ind_list_geno(words[0],flag,&ind_geno);
        }
        else{
            char name[STR_LEN];
            strcpy(name,words[0].c_str());
            IDX *ind_geno_tmp;
            HASH_FIND_STR(ind_geno,name,ind_geno_tmp);
            if(ind_geno_tmp==NULL){
                fprintf(stderr,"\nERROR: Individual is not in genotype file\n");
                throw;
            }
            else{
                if(ind_geno_tmp->index!=flag){
                    fprintf(stderr,"\nERROR: The postion of Individual\n");
                    throw;
                }
            }
        }
        word2geno(words,flag,0,geno,geno_stat,round_geno);
        line.clear();
        flag++;
    }
    input_geno.close();
};
void Greml_ce::read_all_genofiles(const parameters &p,MatrixXi &geno,MatrixXd &geno_stat,ofstream &flog){
        if (p.rawgeno_f==true){
            fprintf(stderr,"\nThe format of genotypic data used is genotypic plain text file!\n");
            string rawgenoFile,rawmapFile;
            mrk_n=-1;
            ind_n=0;
            rawgenoFile=p.genotype;
            if(p.rawmap_f){
                rawmapFile=p.genomap;
                read_map(rawmapFile);
            }
        else{
            fprintf(stderr,"\nERROR: Missing SNP map file!\n");
        }
//        flog<<"Genotype plain text file                           :\t"<<rawgenoFile<<endl;
//        flog<<"SNP map file                                       :\t"<<rawmapFile<<endl;
        ind_n     = calc_sample_size(rawgenoFile);
        geno      = MatrixXi::Zero(ind_n,mrk_n);
        geno_stat = MatrixXd::Zero(12,mrk_n);
        read_rgf(rawgenoFile,geno,geno_stat,round_geno);
    }
    else{
        fprintf(stderr,"\nERROR: Missing genotype file!");
    }
}