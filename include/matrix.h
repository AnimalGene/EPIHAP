#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
typedef struct {
        unsigned int rows;
        unsigned int cols;
        double **mat_data;
} matrix_struct;
matrix_struct *get_matrix_struct(char matrix[]);
void free_matrix(matrix_struct *matrix_to_free);