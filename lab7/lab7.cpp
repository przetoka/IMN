#include <iostream>
#include <cmath>
#include <fstream>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>

#define delta 0.01
#define rho 1
#define mi 1
#define n_x 200
#define n_y 90
#define IT_MAX 20000
#define i_1 50
#define j_1 55
#define y_j1 (delta * j_1)
#define y_ny (delta * n_y)
#define j_2 (j_1 + 2)

void WB_psi(double psi[][n_y + 1], double Q_in, double Q_out){
    for(int j = j_1; j <= n_y; j++)
        psi[0][j] = Q_in / (2 * mi) * (pow((j * delta), 3) / 3.0 - pow((j * delta), 2) * (y_j1 + y_ny) / 2.0 + (j * delta * y_j1 * y_ny));

    for(int j = 0; j <= n_y; j++)
        psi[n_x][j] = Q_out / (2 * mi) * (pow((j * delta), 3) / 3.0 - pow((j * delta), 2) * y_ny / 2.0 ) 
                    + (Q_in * pow((y_j1), 2) *(3 * y_ny - y_j1)) / (12.0 * mi);

    for(int i = 1; i < n_x; i++)
        psi[i][n_y] = psi[0][n_y];

    for(int i = i_1; i < n_x + 1; i++)
        psi[i][0] = psi[0][j_1];

    for(int j = 1; j <= j_1; j++)
        psi[i_1][j] = psi[0][j_1];

    for(int i = 1; i <= i_1; i++)
        psi[i][j_1] = psi[0][j_1];
}

void WB_zeta(double zeta[][n_y + 1], double psi[][n_y + 1], double Q_in, double Q_out){
    for(int j = j_1; j <= n_y; j++)
        zeta[0][j] = Q_in * (2.0 * j * delta - y_j1 - y_ny) / (2.0 * mi);

    for(int j = 0; j <= n_y; j++)
        zeta[n_x][j] = Q_out * (2.0 * j * delta - y_ny) / (2.0 * mi);

    for(int i = 1; i < n_x; i++)
        zeta[i][n_y] = 2.0 * (psi[i][n_y - 1] - psi[i][n_y]) / (delta * delta);

    for(int i = i_1 + 1; i < n_x; i++)
        zeta[i][0] = 2.0 * (psi[i][1] - psi[i][0]) / (delta * delta);

    for(int j = 1; j < j_1; j++)
        zeta[i_1][j] = 2.0 * (psi[i_1 + 1][j] - psi[i_1][j]) / (delta * delta);

    for(int i = 1; i <= i_1; i++)
        zeta[i][j_1] = 2.0 * (psi[i][j_1 + 1] - psi[i][j_1]) / (delta * delta);

    zeta[i_1][j_1] = 0.5 * (zeta[i_1 - 1][j_1] + zeta[i_1][j_1-1]);
}
double get_Q_out(double Q_in){
    return Q_in * (pow(y_ny, 3) - pow(y_j1, 3) - (3.0 * y_j1 * pow(y_ny, 2)) + (3.0 * pow(y_j1, 2) * y_ny)) / pow(y_ny, 3);
}

bool fill(int i, int j){
    if(i <= i_1 && j <= j_1)
        return true;
    if (i == 0 & j_1 < j && j < n_y) 
        return true;
    if (j == n_y)
        return true;
    if(i == n_x)
        return true;
    if(i_1 < i && i < n_x && j == 0)
        return true;
    if(i == i_1 && 0 < j && j < j_1)
        return true;
    if(0 < i && i < i_1 && j ==j_1)
        return true;
    
    return false;
}

void relaxation_NS(double Q_in, std::string filename_psi_zeta, std::string filename_u_v = " "){
    double psi[n_x + 1][n_y + 1]; // ψ
    double zeta[n_x + 1][n_y + 1]; // ζ
    double u[n_x + 1][n_y + 1]; // u
    double v[n_x + 1][n_y + 1]; // v
    int omega;
    double Q_out = get_Q_out(Q_in);

    WB_psi(psi, Q_in, Q_out);
    for(int it = 0; it < IT_MAX; it++){
        if(it < 2000) omega = 0;
        else omega = 1;

        for(int i = 1; i < n_x; i ++){
            for(int j = 1; j < n_y; j++){
                if(!fill(i, j)){
                    psi[i][j] = 0.25 * (psi[i + 1][j] + psi[i - 1][j] + psi[i][j + 1] + psi[i][j - 1] - delta * delta * zeta[i][j]);

                    zeta[i][j] = 0.25 * (zeta[i + 1][j] + zeta[i - 1][j] + zeta[i][j + 1] + zeta[i][j - 1]) 
                    - omega * rho / (16.0 * mi) * ((psi[i][j + 1] - psi[i][j - 1]) * (zeta[i + 1][j] - zeta[i - 1][j])- (psi[i + 1][j] - psi[i - 1][j]) * (zeta[i][j + 1] - zeta[i][j - 1]));

                    u[i][j] =  (psi[i][j + 1] - psi[i][j - 1]) / (2.0 * delta);
                    v[i][j] = -(psi[i + 1][j] - psi[i - 1][j]) / (2.0 * delta);
                }
            }
        }

        WB_zeta(zeta, psi, Q_in, Q_out);
        double gamma = 0;
        for(int i = 1; i < n_x; i++){
            gamma += psi[i + 1][j_2] + psi[i - 1][j_2] + psi[i][j_2 + 1] + psi[i][j_2 - 1] - 4 * psi[i][j_2] - delta * delta * zeta[i][j_2];
        }
        std::cout << gamma << " ";
    }

    std::ofstream file_psi_zeta("data/" + filename_psi_zeta);
    std::ofstream file_u_v("data/" + filename_u_v);

    for(int i = 0; i < n_x; i++){
        for(int j = 0; j < n_y; j++){
            file_psi_zeta << i * delta << "\t" << j * delta << "\t" << psi[i][j] << "\t" << zeta[i][j] <<"\n";
            if(filename_u_v != " ")
                file_u_v << i * delta << "\t" << j * delta << "\t" << u[i][j] << "\t" << v[i][j] << "\n";
        }
    }    
}

int main(){
    mkdir("data", 0777);
    mkdir("plots", 0777);
    
    
    relaxation_NS(-1000, "psi_zeta_-1000.txt", "u_v_-1000.txt");
    relaxation_NS(-4000, "psi_zeta_-4000.txt", "u_v_-4000.txt");
    relaxation_NS(4000, "psi_zeta_4000.txt");

    return 0;
}