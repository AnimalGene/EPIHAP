#ifndef __READ_GENO_H
#define __READ_GENO_H
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <omp.h>
#include <bitset>
#include <cmath>
#include <unistd.h>
#include <pthread.h>
#include "matrix.h"
#include "parameters.h"
#include "tools.h"
using namespace Eigen;
using namespace std;
typedef struct _EFF{
    int            start;
    int            pos;
    int            len;
    double         val;
  }EFF;
typedef struct _Pairwise_result_one{
    string effect_type;
    int snp1;
    int snp2;
    double real_value;
    double effect;
}Pairwise_result_one;
typedef struct _Snp_info{
    string chr;
    string pos;
    string mrk;
}Snp_info;
typedef struct _IDXCOUNT{int idx;double count;} idxcount;
typedef unordered_map<int,idxcount> allele_list;
typedef vector<allele_list> hap_block_allele;
typedef array<int,2> geno_code;
typedef unordered_map<geno_code,idxcount> het_geno_list;
typedef vector<het_geno_list> hap_block_het_geno;
class Greml {
    public:
        int mrk_n;
        int mrk_epi_n;
        int ind_n;
        IDX **cross_idx;
        IDX *ind_idx;
        IDX *idx_ind;
        IDX *ind_geno;
        int round_geno;
        vector<int> cross;
        vector<int> covar;
        int phen_n;
        VectorXd var;
        VectorXd G_diag_mean;
        VectorXd G_diag_mean_sqrt;
        vector<string> var_name;
        vector<int> var_flag;
        vector<EFF> eff_fixed;
        vector<EFF> eff_rand;
        vector<string> hap_name;
        hap_block_allele    hb_allele;
		hap_block_het_geno  hb_het_geno;
		VectorXi            hap_max_key;
        SparseMatrix<double> X;
        SparseMatrix<double> Z;
        SparseMatrix<double> Z1;
        VectorXd Y;
        VectorXi Y_index;
    public:
        unordered_map<int,string>my_yindex;
        void read_phen_pre(const parameters &params,ofstream &flog);
        void read_phen(const parameters &params);
        virtual ~Greml();
};
class Greml_ce : public Greml{
  public:
    double hetero_total_tot;
    double hetero_total2_tot;
    vector<VectorXd>G_diag_tot;
    vector<MatrixXd>G;
    vector<MatrixXi>my_geno_stat_epi_idx;
    unordered_map<string,int>my_ma_type;
    unordered_map<int,Snp_info>mysnpmap;
    map<string,vector<int>>mychrindexmap;
    unordered_map<string,double>var_name_values;
    unordered_map<string,int>var_name_flags;
    unordered_map<string,int>::iterator var_name_flags_iter;
  public:
    void init(const parameters &params);
    void read_map(const string &mapFile);
    void read_rgf(const string &rawgFile,MatrixXi &geno,MatrixXd &geno_stat,int round_geno);
    void read_all_genofiles(const parameters &params,MatrixXi &geno,MatrixXd &geno_stat,ofstream &flog);
    void read_genotype_method_1(const parameters &params,ofstream &flog);
    void each_W_hap_generate(int &block,MatrixXi &geno_hap,MatrixXd &W_hap,hap_block_allele &hb_allele,int n_row,int hap_max_key_j);
    void read_hap(const parameters &params,ofstream &flog);
    void save_wt(const parameters &params,const MatrixXi &geno_org,const MatrixXd &geno_stat,int block_size);
    void jiangyong_egrm_a1(MatrixXd &W,int pow_set, int block_size,MatrixXd &WW);
    void jiangyong_egrm_ad(MatrixXd &WA,MatrixXd &WD,int pow_seta, int pow_setd,int block_size,MatrixXd &WW);
    MatrixXd jiangyong_recurs(vector<MatrixXd>&WWT,int degree);
    MatrixXd jiangyong_recurs_2(vector<vector<MatrixXd>>&WWT,int degreea, int degreed);
    void jiangyong_wwt_mk(vector<MatrixXd> &WW_A_all,vector<MatrixXd> &WW_D_all,vector<vector<MatrixXd>> &WW_AD_all,MatrixXi &geno,MatrixXd &geno_stat,int is_hwd,int order,int allmarkn,int block_size);
    void jiangyong_chrsplit_wwt(MatrixXd &WAAWAAT_intra,MatrixXd &WADWADT_intra,MatrixXd &WDDWDDT_intra,MatrixXi &geno,MatrixXd &geno_stat,int startp,int genol,int block_size,int is_hwd);
    void henderson_chrsplit_WWTAA_WWTDD(MatrixXd &WWTepi,MatrixXi &geno,MatrixXd &geno_stat,string varname,int is_hwd);
    void henderson_chrsplit_WWTAD(MatrixXd &WWTepi,MatrixXi &geno,MatrixXd &geno_stat,int is_hwd);
    void matrix_rank_print(MatrixXd &Gmat,string mname);
    void read_genotype_method_2(const parameters &params,ofstream &flog);
    void load_genotype(const parameters &params,ofstream &flog);
    void variance_component_calculation(const parameters &params,SparseMatrix<double> &Xt,SparseMatrix<double> &Zt,int ind_n,int ve_i,ofstream &flog);
    void update_dif_max(VectorXd &var,VectorXd &var_tmp,VectorXd &dif_tmp,double &dif_max);
    void get_h2(MatrixXd &IAI,VectorXd &h2,VectorXd &h2se);
    void output_greml_iter_title(ofstream &var_output);
    void output_greml_iter_var(ofstream &var_output,int iter,VectorXd &dif_tmp,ofstream &flog);
    void output_greml_inv_ai(ofstream &var_output,MatrixXd &IAI,ofstream &flog);
    void output_greml_her(VectorXd &h2,VectorXd &h2dif_tmp,ofstream &flog);
    void output_greml_her_a(ofstream &var_output,VectorXd &h2,VectorXd &h2se,ofstream &flog);
    void gblup_ce(const parameters &params,SparseMatrix<double> &Xt,SparseMatrix<double> &Z2t,int ind_n,int ve_i,ofstream &flog);
    MatrixXd get_epi_eff_each_type(string &W1Tfile,string &W2Tfile,string var_epi_type,VectorXd &r);
    void get_single_mrk_eff(const parameters &params,vector<VectorXd> &single_mrk_eff_all);
    void get_mrk_eff_hap(const parameters &params,VectorXd &r);
    void get_epi_mrk_her_m_std(const parameters &params,MatrixXd &epi_mrk_eff,VectorXd &her_all, int indx);
    void get_epi_mrk_eff_AA_or_DD(const parameters &params,string eff_type,MatrixXd &epi_mrk_eff,VectorXd &her_all,Pairwise_result_one *epistasis_effects);
    void get_epi_mrk_eff_AD(const parameters &params,string eff_type,MatrixXd &epi_mrk_eff1,MatrixXd &epi_mrk_eff2,VectorXd &her_all1,VectorXd &her_all2,Pairwise_result_one *epistasis_effects);
    vector<double> mean_std_single_snp(VectorXd &vec);
    vector<double> mean_std_epi_snp(MatrixXd &vec);
    void insert_epistasis_effect(const parameters &params,Pairwise_result_one *epistasis_effects,int &n_inserted,Pairwise_result_one epistasis_effect);
    void calc_each_epi_eff(const parameters &params,MatrixXd &mrk_epi_eff,string eff_type,Pairwise_result_one *epistasis_effects);
    void get_SNP_eff(const parameters &params,VectorXd &r);
    void reml(const parameters &params,ofstream &flog);
    ~Greml_ce(){};
  };
#endif