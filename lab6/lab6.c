# include <stdio.h>
# include <stddef.h>
# include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include "mgmres.h"

#define itr_max 500
#define mr 500
#define tol_abs (1e-8)
#define tol_rel (1e-8)


int j_val(int l, int n_x){
    return floor(l / (n_x + 1));
}

int i_val(int l, int j, int n_x){
    return l - j*(n_x + 1);
}

double el(int l, int n_x, double e_1, double e_2){
    int j = j_val(l, n_x);
    int i = i_val(l, j, n_x);

    if(i <= n_x / 2)
        return e_1;
    return e_2;
}

double rho1(double x, double y, int x_max, int y_max){
    double rho = x_max / 10;
    return exp(- pow((x - 0.25 * x_max), 2)/pow(rho, 2) - pow((y - 0.5 * y_max), 2)/pow(rho, 2));
}
double rho2(double x, double y, int x_max, int y_max){
    double rho = x_max / 10;
    return -exp(- pow((x - 0.75 * x_max), 2)/pow(rho, 2) - pow((y - 0.5 * y_max), 2)/pow(rho, 2));
}

void solve(int ex, double delta, double e_1, double e_2, int V_1, int V_2, int n_x, const char filename[]){
    int N = ((n_x + 1) * (n_x + 1));
    int ia[N + 1];
    memset(ia, -1, N + 1);
    int ja[5 * N];
    double a[5 * N];
    double b[N];
    double V[N];
    memset(V, 0, N);
    int x_max = n_x * delta;
    int y_max = n_x * delta;

    int k = -1;
    for(int l = 0; l < N; l++){
        int edge = 0; 
        int vb = 0;
        int j = j_val(l, n_x);
        int i = i_val(l, j, n_x);

        if(i == 0 || i == n_x){
            edge = 1;
            vb = V_1;
        }
        if(j == n_x || j == 0){
            edge = 1;
            vb = V_2;
        }
        
        if(ex == 6)
            b[l] = -(rho1(i * delta, j * delta, x_max, y_max) + rho2(i * delta, j * delta, x_max, y_max));
        else{
            b[l] = 0;
        }
        if(edge == 1)
            b[l] = vb;
        
        ia[l] = -1;
        double e_l = el(l, n_x, e_1, e_2);
        double el_1 = el(l + 1, n_x, e_1, e_2);
        double el_n = el(l + n_x + 1, n_x, e_1, e_2);
        if(l - n_x - 1 >= 0 && edge == 0){
            k++;
            if(ia[l] < 0)
                ia[l] = k;
            a[k] = (e_l * 1.0 / delta / delta);
            ja[k] = l - n_x - 1;
        }
        if(l  - 1 >= 0 && edge == 0){
            k++;
            if(ia[l] < 0)
                ia[l] = k;
            a[k] = (e_l * 1.0 / delta / delta);
            ja[k] = l - 1;
        }
        k++;
        if(ia[l] < 0)
            ia[l] = k;
        if(edge == 0){
            a[k] = (-(2 * e_l + el_1 + el_n) * 1.0 / delta / delta);
        }
            
        else{
            a[k] = 1;
        }
        ja[k] = l;
        
        if(l < N && edge == 0){
            k++;
            a[k] = (el_1 / delta / delta);
            ja[k] = l + 1;
        }

        if(l < N - n_x - 1 && edge == 0){
            k++;
            a[k] = (el_n / delta / delta);
            ja[k] = l + n_x + 1;
        }
    }
    int nz_num = k + 1;
    ia[N] = nz_num;

    if(ex == 4){
        FILE *file_a = fopen("data/b.txt", "w");
        FILE *file_b = fopen("data/a.txt", "w");

        fprintf(file_b, "# l \t i_l \t j_l \t b[l] \n");
        fprintf(file_a, "# k \t a[k] \n");

        for(int l = 0; l < N; l++){
            int j = j_val(l, n_x);
            int i = i_val(l, j, n_x);
            fprintf(file_b, "# %d \t %d \t %d \t %lf \n", l, i, j, b[l]);
        }
        for(k = 0; k < N * 5; k++){
            fprintf(file_a, "# %d \t %lf \n", k, a[k]);
        }
    }

    if(ex != 4){
        pmgmres_ilu_cr ( N, nz_num, ia, ja, a, V, b, itr_max, mr, tol_abs, tol_rel);
        FILE *file = fopen(filename, "w");
        for(int l = 0; l < N; l++){
            int j = j_val(l, n_x);
            int i = i_val(l, j, n_x);
            fprintf(file, "%lf \t %lf \t %lf \n", i*delta, j*delta, V[l]);
        }
        
    }
}

int main(){
    
    mkdir("data", 0777);
    mkdir("plots", 0777);

    solve(4, 0.1, 1.0, 1.0, 10, -10, 4, "");
    
    solve(5, 0.1, 1.0, 1.0, 10, -10, 50, "data/5_a.txt");
    solve(5, 0.1, 1.0, 1.0, 10, -10, 100, "data/5_b.txt");
    solve(5, 0.1, 1.0, 1.0, 10, -10, 200, "data/5_c.txt");

    solve(6, 0.1, 1.0, 1.0, 0, 0, 100, "data/6_a.txt");
    solve(6, 0.1, 1.0, 2.0, 0, 0, 100, "data/6_b.txt");
    solve(6, 0.1, 1.0, 10.0, 0, 0, 100, "data/6_c.txt");
    return 0;
}