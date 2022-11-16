#include <iostream>
#include <cmath>
#include <fstream>

#define e       1.0
#define n_x     128
#define n_y     128
#define delta   0.2
#define x_max   delta * n_x
#define y_max   delta * n_y
#define sigma_x 0.1 * x_max
#define sigma_y 0.1 * y_max
#define TOL     1e-8
#define V1      10.0
#define V2      0.0
double rho[n_x + 1][n_y + 1];
int k[] = {16, 8, 4, 2, 1};

double V_B1(double y){ return sin(M_PI * y / y_max); }
double V_B2(double x){ return - sin(2* M_PI * x / x_max); }
double V_B3(double y){ return sin(M_PI * y / y_max); }
double V_B4(double x){ return sin(2* M_PI * x / x_max); }

double s(double V[][n_y + 1], int k){
    double S = 0.0;
    for(int i = 0; i < n_x - k ; i += k){
        for(int j = 0; j < n_y - k; j += k){
            S += (delta * k) * (delta * k) / 2.0 * ( 
                  pow(((V[i + k][j] - V[i][j]) / (2.0 * k * delta) + (V[i + k][j + k] - V[i][j + k]) / (2.0 * k * delta)), 2)
                + pow(((V[i][j + k] - V[i][j]) / (2.0 * k * delta) + (V[i + k][j + k] - V[i + k][j]) / (2.0 * k * delta)), 2));
        }
    }
    return S;
}

void multimesh_relaxation(int k, std::string filename_sum, std::string filename_rel){

}


int main(){

    return 0;
}