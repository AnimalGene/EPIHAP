
#include "matrix.h"

matrix_struct *get_matrix_struct(char matrix[]) {
    matrix_struct *m =(matrix_struct*)malloc(sizeof(matrix_struct));
    m->rows = 0;
    m->cols = 0;
    FILE* myfile = fopen(matrix, "r");
    if(myfile == NULL) {
        printf("Error: The file you entered could not be found.\n");
        exit(EXIT_FAILURE);
    }
    // get the rows and columns
    int ch = 0;
    do {
        ch = fgetc(myfile);
        // count the columns at the first line (looking for "\t")
        if(m->rows == 0 && ch == '\t')
            m->cols++;
        // count the rows with "\n"
        if(ch == '\n')
            m->rows++;
            
    } while (ch != EOF);
    // write rows and cols to struct
    m->cols++;
    // allocate memory for matrix data
    m->mat_data = (double**)calloc(m->rows, sizeof(double*));
    int i;
    for(i=0; i < m->rows; ++i) m->mat_data[i]=(double*)calloc(m->cols, sizeof(double));
    rewind(myfile);
    int x,y;
    // fill matrix with data
    for(x = 0; x < m->rows; x++) {
        for(y = 0; y < m->cols; y++) {
            if (!fscanf(myfile, "%lf", &m->mat_data[x][y]))
            break;
        }
    }
    fclose(myfile);
    return m;
}
void free_matrix(matrix_struct *matrix_to_free) {
    for(int i = 0; i < matrix_to_free->rows; i++) {
        free(matrix_to_free->mat_data[i]);
    }
    free(matrix_to_free->mat_data);
    free(matrix_to_free);
}

