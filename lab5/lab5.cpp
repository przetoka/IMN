#include <iostream>
#include <cmath>
#include <fstream>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>

#define n_x     128
#define n_y     128
#define delta   0.2
#define x_max   (delta * n_x)
#define y_max   (delta * n_y)
#define TOL     1e-8
int it = 0;
int k[] = {16, 8, 4, 2, 1};

double V_B1(double y){ return sin(M_PI * y / y_max); }
double V_B2(double x){ return - sin(2* M_PI * x / x_max); }
double V_B3(double y){ return sin(M_PI * y / y_max); }
double V_B4(double x){ return sin(2* M_PI * x / x_max); }

double s(double V[][n_y + 1], int k){
    double S = 0.0;
    for(int i = 0; i <= n_x - k ; i += k){
        for(int j = 0; j <= n_y - k; j += k){
            S +=    pow((delta * k), 2) * 0.5 * ( 
                    pow(((V[i + k][j] - V[i][j] + V[i + k][j + k] - V[i][j + k]) / (2.0 * k * delta)), 2)
                +   pow(((V[i][j + k] - V[i][j] + V[i + k][j + k] - V[i + k][j]) / (2.0 * k * delta)), 2));
        }
    }
    return S;
}

int main(){
    if (mkdir("data", 0777) == -1)   std::cerr << "Error :  " << strerror(errno) << std::endl;
    else    std::cout << "Directory created\n";

    if (mkdir("plots", 0777) == -1)  std::cerr << "Error :  " << strerror(errno) << std::endl;
    else    std::cout << "Directory created\n";
    
    double V[n_x + 1][n_y + 1]   = {0.0};
    for(int i = 0; i <= n_x; i++){
        V[i][n_y] = V_B2(i * delta);
        V[i][0] = V_B4(i * delta);
    }

    for(int j = 0; j <= n_y; j++){
        V[0][j] = V_B1(j * delta);
        V[n_x][j] = V_B3(j * delta);
    }    


    int k = 16;
    while(k >= 1){
        std::ofstream file_V("data/map" + std::to_string(k) + ".txt");
        std::ofstream file_s("data/sum" + std::to_string(k) + ".txt");
        
        double s_prev;
        double s_next = s(V, k);

        do{
            for(int i = k; i <= n_x - k ; i += k){
                for(int j = k; j <= n_y - k ; j += k){
                    V[i][j] = 0.25 * (V[i + k][j] + V[i - k][j] + V[i][j + k] + V[i][j - k]);
                }
            }
            
            s_prev = s_next;
            s_next = s(V, k);
            
            it ++;
            file_s << it << " " << s_next << "\n";
        }while(fabs((s_next - s_prev)/s_prev) > TOL);

        for(int i = 0; i <= n_x; i += k){
            for(int j = 0; j <= n_y; j += k){
                file_V << i * delta << " " << j * delta << " " << V[i][j] << "\n";
            }
        }
        
        for(int i = 0; i < n_x; i += k){
            for(int j = 0; j < n_y; j += k){
                V[i + k/2][j + k/2] = 0.25 * (V[i][j] + V[i + k][j] + V[i][j + k] + V[i + k][j + k]);
                V[i + k][j + k/2]   = 0.5 * (V[i + k][j] + V[i + k][j + k]);
                V[i + k/2][j + k]   = 0.5 * (V[i][j + k] + V[i + k][j + k]);
                V[i + k/2][j]       = 0.5 * (V[i][j] + V[i + k][j]);
                V[i][j + k/2]       = 0.5 * (V[i][j] + V[i][j + k]);
            }
        }
        
        for(int i = 0; i <= n_x; i++){
            V[i][n_y] = V_B2(i * delta);
            V[i][0] = V_B4(i * delta);
        }

        for(int j = 0; j <= n_y; j++){
            V[0][j] = V_B1(j * delta);
            V[n_x][j] = V_B3(j * delta);
        }

        k = k / 2;
    }
    std::cout << it << "\n";
    return 0;
}