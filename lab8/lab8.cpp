#include <iostream>
#include <cmath>
#include <fstream>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string>

#define delta 0.01
#define n_x 400
#define n_y 90
#define i_1 200
#define i_2 210
#define j_1 50
#define rho 10
#define x_A 0.45
#define y_A 0.45
#define sigma 0.1
#define it_max 10000
#define x_max (n_x * delta)
#define y_max (n_y * delta)

double psi[n_x + 1][n_y + 1] = {0.0};
double v_x[n_x + 1][n_y + 1] = {0.0};
double v_y[n_x + 1][n_y + 1] = {0.0};
double u_0[n_x + 1][n_y + 1] = {0.0};
double u_1[n_x + 1][n_y + 1] = {0.0};

double delta_t, v_max;

void get_data(){
    double val;
    int i, j;
    std::ifstream file;
    file.open("psi.dat");
    if (!file)
        std::cout << "File Not Found." << std::endl;
    else {
        while(file >> val){
            i = val;
            file >> val;
            j = val;
            file >> val;
            psi[i][j] = fabs(val);
        }
    }
}

double v(){
    for(int i = 1; i < n_x; i++){
        for(int j = 1; j < n_y; j++){
            v_x[i][j] = (psi[i][j + 1] - psi[i][j - 1]) / (delta * 2);
            v_y[i][j] = - (psi[i + 1][j] - psi[i - 1][j]) / (delta * 2);
        }
        v_x[i][0] = v_y[i][n_y] = 0;
    }

    for(int i = i_1; i <= i_2; i++){
        for(int j = 0; j <= j_1; j++){
            v_x[i][j] = v_y[i][j] = 0;
        }
    }

    for(int j = 0; j <= n_y; j++){
        v_x[0][j] = v_x[1][j];
        v_x[n_x][j] = v_x[n_x - 1][j];
    }
    std::ofstream file_v_x_y("data/v_x_y.txt");
    for(int i = 0; i <= n_x; i++){
        for(int j = 0; j <= n_y; j++){
            file_v_x_y << i*delta << "\t" << j*delta << "\t" << v_x[i][j] << "\t" <<  v_y[i][j] << "\n";
        }
    }

    v_max = sqrt(v_x[0][0] * v_x[0][0] + v_y[0][0] * v_y[0][0]);
    int i_max = 0, j_max = 0;
    for(int i = 0; i < n_x; i++){
        for(int j = 0; j < n_y; j++){
            double temp = sqrt(v_x[i][j] * v_x[i][j] + v_y[i][j] * v_y[i][j]);
            if(temp > v_max){
                v_max = temp;
                i_max = i;
                j_max = j;
            }
        }
    }
    std::cout << "Punkt i, j, dla którego moduł prędkości |v_i_j| jest największy: v[" << i_max << ", " << j_max << "] = " << v_max << "\n";
    return v_max;
}


void advection_diffusion(double D, std::string filename_c, std::string filename_x_sr, std::string f_name){
    std::ofstream file_c("data/" + filename_c);
    std::ofstream file_x_sr("data/" + filename_x_sr);


    for(int i = 0; i <= n_x; i++){
        for(int j = 0; j <= n_y; j++)
            u_0[i][j] = exp(-((i * delta - x_A) * (i * delta - x_A) + (j * delta - y_A) * (j * delta- y_A)) / (2 * sigma * sigma)) / (2 * M_PI * sigma * sigma);
    }

    for(int it = 0; it < it_max; it++){

        for(int i = 0; i <= n_x; i++){
            for(int j = 0; j <= n_y; j++)
                u_1[i][j] = u_0[i][j];
        }

        for(int k = 1; k <= 20; k++){
            for(int i = 0; i <= n_x; i++){
                for(int j = 0; j <= n_y; j++)
                    if(i_1 <= i && i <= i_2 && j <=j_1){
                        continue;
                    }

                    else if (i == 0)
                        u_1[i][j] = 1 / (1 + ((2 * D * delta_t) / (delta * delta))) * (u_0[i][j] 
                            - (delta_t / 2 * v_x[i][j] * ((u_0[i + 1][j] - u_0[n_x][j]) / (2 * delta) + (u_1[i + 1][j] - u_1[n_x][j]) / (2 * delta)))
                            - (delta_t / 2 * v_y[i][j] * ((u_0[i][j + 1] - u_0[i][j - 1]) / (2 * delta) + (u_1[i][j + 1] - u_1[i][j - 1]) / (2 * delta)))
                            + (delta_t / 2 * D * ((u_0[i + 1][j] + u_0[n_x][j] + u_0[i][j + 1] + u_0[i][j - 1] - 4 * u_0[i][j]) / (delta * delta)
                                              + ((u_1[i + 1][j] + u_1[n_x][j] + u_1[i][j + 1] + u_1[i][j - 1])/ (delta * delta)))));

                    else if (i == n_x)
                        u_1[i][j] = 1/(1 + ((2 * D *delta_t) / (delta * delta))) * (u_0[i][j] 
                            - (delta_t / 2 * v_x[i][j] * ((u_0[0][j] - u_0[i - 1][j]) / (2 * delta) + (u_1[0][j] - u_1[i - 1][j]) / (2 * delta)))
                            - (delta_t / 2 * v_y[i][j] * ((u_0[i][j + 1] - u_0[i][j - 1]) / (2 * delta) + (u_1[i][j + 1] - u_1[i][j - 1]) / (2 * delta)))
                            + (delta_t / 2 * D * ((u_0[0][j] + u_0[i - 1][j] + u_0[i][j + 1] + u_0[i][j - 1] - 4 * u_0[i][j]) / (delta * delta)
                                              + ((u_1[0][j] + u_1[i - 1][j] + u_1[i][j + 1] + u_1[i][j - 1])/ (delta * delta)))));

                    else
                        u_1[i][j] = 1/(1 + ((2 * D *delta_t) / (delta * delta))) * (u_0[i][j] 
                            - (delta_t / 2 * v_x[i][j] * ((u_0[i + 1][j] - u_0[i - 1][j]) / (2 * delta) + (u_1[i + 1][j] - u_1[i - 1][j]) / (2 * delta))) 
                            - (delta_t / 2 * v_y[i][j] * ((u_0[i][j + 1] - u_0[i][j - 1]) / (2 * delta) + (u_1[i][j + 1] - u_1[i][j - 1]) / (2 * delta)))
                            + (delta_t / 2 * D * ((u_0[i + 1][j] + u_0[i - 1][j] + u_0[i][j + 1] + u_0[i][j - 1] - 4 * u_0[i][j]) / (delta * delta)
                                              + ((u_1[i + 1][j] + u_1[i - 1][j] + u_1[i][j + 1] + u_1[i][j - 1])/ (delta * delta)))));
            }
        }
        double c = 0, x_sr = 0;

        for(int i = 0; i <= n_x; i++){
            for(int j = 0; j <= n_y; j++){
                c += (u_0[i][j] * delta * delta);
                x_sr += ((i * delta) * u_0[i][j] * delta * delta);
            }
        }

        file_c << (it * delta_t) << "\t" << c << "\n";
        file_x_sr << (it * delta_t) << "\t" << x_sr << "\n";
        for(int i = 0; i <= n_x; i++){
            for(int j = 0; j <= n_y; j++)
                u_0[i][j] = u_1[i][j];
        }
        
        if(it % 200){
            std::string name = "data/maps/map" + std::to_string(it/200) + f_name;
                
            std::ofstream file_map(name);
            for(int i = 0; i <= n_x; i++){
                for(int j = 0; j <= n_y; j++)
                    file_map << i * delta << "\t" << j * delta << "\t" << u_1[i][j] << "\n";
            }
        }
    }
}

int main(){
    mkdir("data", 0777);
    mkdir("data/maps", 0777);
    mkdir("plots", 0777);
    mkdir("plots/maps", 0777);
    
    get_data();
    v_max = v();
    delta_t = (delta / (4 * v_max));
    std::cout << delta_t;
    

    advection_diffusion(0.0, "c_1.txt", "x_sr_1.txt", "_00.txt");
    advection_diffusion(0.1, "c_2.txt", "x_sr_2.txt", "_01.txt");
    
    return 0;
}