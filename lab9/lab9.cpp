#include <iostream>
#include <cmath>
#include <fstream>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string>
// #include "/usr/include/gsl/gsl_math.h"
// #include "/usr/include/gsl/gsl_linalg.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

#define n_x 40
#define n_y 40
#define n ((n_x + 1)*(n_y + 1))
#define delta 1
#define delta_t 1
#define T_A 40
#define T_B 0
#define T_C 30
#define T_D 0
#define k_B 0.1
#define k_D 0.6
#define IT_MAX 2000


gsl_vector* multiplyMatrixVector(gsl_matrix* m, gsl_vector* v){
    gsl_vector *temp = gsl_vector_calloc(n + 1);
    for(int i = 0; i <= n; i++){
        double val = 0;
        for(int j = 0; j <= n; j++){
            val += gsl_matrix_get(m, i, j) * gsl_vector_get(v, j);
        }
        gsl_vector_set(temp, i, val);
    }

    return temp;
}

gsl_vector* addVectors(gsl_vector* v1, gsl_vector* v2){
    gsl_vector* temp = gsl_vector_calloc(n + 1);
    for(int i = 0; i <= n; i++){
        gsl_vector_set(temp, i, gsl_vector_get(v1, i) + gsl_vector_get(v2, i));
    }
    return temp;
}

int main(){
    mkdir("data", 0777);
    mkdir("plots", 0777);
    
    int s;
    gsl_matrix *A = gsl_matrix_calloc(n + 1, n + 1);
    gsl_matrix *B = gsl_matrix_calloc(n + 1, n + 1);
    gsl_vector *c = gsl_vector_calloc(n + 1);
    gsl_vector *d = gsl_vector_calloc(n + 1);
    gsl_vector *T = gsl_vector_calloc(n + 1);
    gsl_permutation *p = gsl_permutation_calloc(n + 1);

    for(int i = 1; i < n_x; i++){
        for(int j = 1; j < n_y; j++){
            int l = i + j * (n_x + 1);
            gsl_matrix_set(A, l, l - n_x - 1,   delta_t / (2 * delta * delta));
            gsl_matrix_set(A, l, l - 1,         delta_t / (2 * delta * delta));
            gsl_matrix_set(A, l, l + 1,         delta_t / (2 * delta * delta));
            gsl_matrix_set(A, l, l + n_x + 1,   delta_t / (2 * delta * delta));
            gsl_matrix_set(A, l, l,             - (2 * delta_t) / (delta * delta) - 1);

            gsl_matrix_set(B, l, l - n_x - 1,   - delta_t / (2 * delta * delta));
            gsl_matrix_set(B, l, l - 1,         - delta_t / (2 * delta * delta));
            gsl_matrix_set(B, l, l + 1,         - delta_t / (2 * delta * delta));
            gsl_matrix_set(B, l, l + n_x + 1,   - delta_t / (2 * delta * delta));
            gsl_matrix_set(B, l, l,             (2 * delta_t) / (delta * delta) - 1);
        }
            int j = n_y;
            int l = i + j * (n_x + 1);
            gsl_matrix_set(A, l, l - n_x - 1,   - 1 / (k_B * delta));
            gsl_matrix_set(A, l, l,             1 + 1 / (k_B * delta));
            gsl_vector_set(c, l, T_B);
            for(j = 0; j <= n; j++){
                gsl_matrix_set(B, l, j, 0);
            }

            j = 0;
            l = i + j * (n_x + 1);
            gsl_matrix_set(A, l, l + n_x + 1,   1 + 1 / (k_D * delta));
            gsl_matrix_set(A, l, l,             1 + 1 / (k_D * delta));
            gsl_vector_set(c, l, T_D);
            for(j = 0; j <= n; j++){
                gsl_matrix_set(B, l, j, 0);
            }
    }

    for(int j = 0; j <= n_y; j++){
        int i = 0;
        int l = i + j * (n_x + 1);
        gsl_matrix_set(A, l, l, 1);
        gsl_matrix_set(B, l, l, 1);
        gsl_vector_set(c, l, 0);

        i = n_x;
        l = i + j * (n_x + 1);
        gsl_matrix_set(A, l, l, 1);
        gsl_matrix_set(B, l, l, 1);
        gsl_vector_set(c, l, 0);
    }

    for(int j = 0; j <= n_y; j++){
        int i = 0;
        int l = i + j * (n_x + 1);
        gsl_vector_set(T, l, T_A);
    }

    for(int j = 0; j <= n_y; j++){
        int i = n_x;
        int l = i + j * (n_x + 1);
        gsl_vector_set(T, l, T_C);
    }

    for(int i = 1; i < n_x; i++){
        for(int j = 1; j < n_y; j++){
            int l = i + j * (n_x + 1);
            gsl_vector_set(T, l, 0);
        }
    }

    gsl_linalg_LU_decomp(A, p, &s);
    p = gsl_permutation_calloc(n + 1);

    for(int it = 0; it <= IT_MAX; it++){
        d = addVectors(multiplyMatrixVector(B, T), c);

        // gsl_linalg_LU_solve(A, p, T, d);

        if(it % 200){
            std::string name = "data/map" + std::to_string(it/200) + ".txt";
                
            std::ofstream file_map(name);
            for(int l = 0; l <= n; l++){
                int j = floor(l / (n_x + 1));
                int i = l - j * (n_x + 1);
                
                file_map << i * delta << " " << j * delta << " " << gsl_vector_get(T, l) << "\n";
            }
        }
    }

    return 0;
}