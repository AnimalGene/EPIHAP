#ifndef __TOOLS_H
#define __TOOLS_H
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <unordered_map>
#include <map>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "uthash.h"
#include "omp.h"
#if defined MKL
#include <mkl_cblas.h>
#include <mkl_lapacke.h>
#else
#include <cblas.h>
#include <lapacke.h>
#endif
#define  STR_LEN  50
using namespace std;
using namespace Eigen;
typedef struct _IDX{
    char           name[STR_LEN];
    double         index;
    UT_hash_handle hh;
  }IDX;
vector<string>str_split(const string &delimiter,const string &explodeme);
void output_geno_summary_header(string file);
typedef Eigen::Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXd_r;
void get_ind_list_geno(string word,int row,IDX **ind_geno);
void genostat_update(MatrixXd &geno_stat);
void genostat_update_epi(MatrixXd &geno_stat,MatrixXd &geno_stat_epi,int snp1_index,int mrk_c);
void get_het_sum(MatrixXd &geno_stat,int is_hwd,double *hetero_total,double *hetero_total2);
void geno_transform(const MatrixXi& geno_org,const MatrixXd& geno_stat,string geno_type,int is_hwd,MatrixXd& geno_new);
void save_ind_geno(string filename,IDX **ind_geno);
void load_ind_geno(string filename,IDX **ind_geno);
void read_plain_txt(const char* filename, MatrixXd& matrix);
void word2geno(vector<string> &words,int row,int col_start,MatrixXi &geno,MatrixXd &geno_stat,int round_geno);
int binomialCoefficients(int n, int k);
int calc_sample_size(const string &filename);
void saveData(string fileName, MatrixXd &matrix);
template<class Matrix>
void write_binary(const char* filename, const Matrix& matrix){
    ofstream out(filename,ios::out | ios::binary | ios::trunc);
    typename Matrix::Index rows=matrix.rows(), cols=matrix.cols();
	out.write((char*) (&rows), sizeof(typename Matrix::Index));
	out.write((char*) (&cols), sizeof(typename Matrix::Index));
	out.write((char*) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar));
	out.close();
}
template<class Matrix>
void read_binary(const char* filename, Matrix& matrix){
    ifstream in(filename,ios::in | ios::binary);
	typename Matrix::Index rows=0, cols=0;
	in.read((char*) (&rows),sizeof(typename Matrix::Index));
	in.read((char*) (&cols),sizeof(typename Matrix::Index));
	matrix.resize(rows, cols);
	in.read( (char *) matrix.data() , rows*cols*sizeof(typename Matrix::Scalar));
	in.close();
}
template<>
struct std::hash<array<int, 2>>{
    typedef array<int, 2> argument_type;
	typedef size_t result_type;
	result_type operator()(const argument_type& a)const{
	hash<int> hasher;
	result_type h = 0;
	for (result_type i = 0; i < 2; ++i){
	    h = h * 31 + hasher(a[i]);
	}
	return h;
    }
};
#endif