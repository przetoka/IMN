#include <iostream>
#include <cmath>
#include <fstream>

#define e       1.0
#define n_x     150
#define n_y     100
#define delta   0.1
#define x_max   delta * n_x
#define y_max   delta * n_y
#define sigma_x 0.1 * x_max
#define sigma_y 0.1 * y_max
#define TOL     1e-8
#define V1      10.0
#define V2      0.0
double rho[n_x + 1][n_y + 1];


double rho_1(double i, double j){
    return exp(
        -((i * delta - 0.35 * x_max) * (i * delta - 0.35 * x_max)) / (sigma_x * sigma_x) 
        -((j * delta - 0.5 * y_max) *(j * delta - 0.5 * y_max)) / (sigma_y * sigma_y));
}

double rho_2(double i, double j){
    return -exp(
        -((i * delta - 0.65 * x_max) * (i * delta - 0.65 * x_max)) / (sigma_x * sigma_x) 
        -((j * delta - 0.5 * y_max) * (j * delta - 0.5 * y_max)) / (sigma_y * sigma_y));
}

double s(double V[][n_y + 1]){
    double S = 0.0;
    for(int i = 0; i < n_x; i ++){
        for(int j = 0; j < n_y; j++){
            S += delta * delta * (0.5 * ((V[i + 1][j] - V[i][j]) / delta) * ((V[i + 1][j] - V[i][j]) / delta) 
                + 0.5 * ((V[i][j + 1] - V[i][j]) / delta) * ((V[i][j + 1] - V[i][j]) / delta)
                - rho[i][j] * V[i][j]);
        }
    }
    return S;
}

void global_relaxation(double omega_g, std::string filename_map, std::string filename_err, std::string filename_s){
    double V_n[n_x + 1][n_y + 1] = {0.0};
    double err[n_x + 1][n_y + 1] = {0.0};
    double V[n_x + 1][n_y + 1]   = {0.0};
    std::ofstream file_map("data/" + filename_map);
    std::ofstream file_err("data/" + filename_err);
    std::ofstream file_s("data/" + filename_s);

    for(int i = 0; i <= n_x; i++){
        V_n[i][0] = V1;
        V_n[i][n_y] = V2;
        V[i][0] = V1;
        V[i][n_y] = V2;
    }

    double s_prev;
    double s_next = 0.0;
    int iter = 0;

    do{
        for(int i = 1; i < n_x; i++){
            for(int j = 1; j < n_y; j++)
                V_n[i][j] = 0.25 * (V[i + 1][j] + V[i - 1][j] + V[i][j + 1] + V[i][j - 1] + delta * delta * rho[i][j] / e);
        }
        for(int j = 1; j < n_y; j++){ 
            V_n[0][j] = V_n[1][j];
            V_n[n_x][j] = V_n[n_x - 1][j];
        }
        for(int i = 0; i <= n_x; i++){
            for(int j = 0; j <= n_y; j++)
                V[i][j] = (1 - omega_g) * V[i][j] + omega_g * V_n[i][j];
        }
        
        s_prev = s_next;
        s_next = s(V);
        iter ++;
        file_s << iter << " " << s_next << "\n";
    }while(fabs((s_next - s_prev)/s_prev) > TOL);

    for(int i = 1; i < n_x; i++){
        for(int j = 1; j < n_y; j++){
            err[i][j] = (V[i + 1][j] - 2.0 * V[i][j] + V[i - 1][j])/(delta * delta) + 
                (V[i][j + 1] - 2.0 * V[i][j] + V[i][j - 1])/(delta * delta) + 
                rho[i][j] / e;
            file_map << i * delta << " " << j * delta << " " << V[i][j] << "\n";
            file_err << i * delta << " " << j * delta << " " << err[i][j]<< "\n";
        }
    }

}

void local_relaxation(double omega_l, std::string filename){
    double V[n_x + 1][n_y + 1] = {0.0};
    std::ofstream file("data/" + filename);

    for(int i = 0; i <= n_x; i++){
        V[i][0] = V1;
        V[i][n_y] = V2;
    }

    double s_prev;
    double s_next = 0.0;
    int iter = 0;

    do{
        for(int i = 1; i < n_x; i++){
            for(int j = 1; j < n_y; j++)
                V[i][j] = (1.0 - omega_l) * V[i][j] + omega_l * 0.25 * (V[i + 1][j] + V[i - 1][j] + V[i][j + 1] + V[i][j - 1] + delta * delta * rho[i][j] / e);
        }
        for(int j = 1; j < n_y; j++){ 
            V[0][j] = V[1][j];
            V[n_x][j] = V[n_x - 1][j];
        }

        s_prev = s_next;
        s_next = s(V);
        iter ++;
        file << iter << " " << s_next << "\n";
    }while(fabs((s_next - s_prev)/s_prev) > TOL);
}

int main(){
    double omega_L[] = {1.0, 1.4, 1.8, 1.9};
    double omega_G[] = {0.6, 1.0};
    

    for (int i = 0; i < n_x + 1; i++){
        for(int j = 0; j < n_y + 1; j++){
            rho[i][j] = rho_1(i, j) + rho_2(i, j);
        }
    }

    std::string sum_l[] = {"local_sum1.txt", "local_sum2.txt", "local_sum3.txt", "local_sum4.txt"};
    std::string map_g[] = {"global_map1.txt", "global_map2.txt"};
    std::string err_g[] = {"global_err1.txt", "global_err2.txt"};
    std::string sum_g[] = {"global_sum1.txt", "global_sum2.txt"};

    for(int i = 0; i < 4; i++){
        local_relaxation(omega_L[i], sum_l[i]);
    }
    for(int i = 0; i < 2; i++){
        global_relaxation(omega_G[i],  map_g[i], err_g[i], sum_g[i]);
    }


    return 0;
}