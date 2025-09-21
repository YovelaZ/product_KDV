#include "SLAM.h"
#include <cmath>
#include <vector>



double compute_SLWmatrix_2(int m, double p_x, double p_y, double q_y, int k_type_y, double b_y){
    double value;
    if(p_x == 0 && m==0){
        value = compute_density(q_y,p_y,k_type_y,b_y);
    } else {
        value = pow(p_x, m)*compute_density(q_y,p_y,k_type_y,b_y);
    }
    return value;
}
    
double compute_SLWmatrix_3(int m, double p_x, double p_y, double p_z, double q_y, double q_z, int k_type_y, int k_type_z, double b_y, double b_z){
    double value;
    if(p_x == 0 && m==0){
        value = compute_density(q_y,p_y,k_type_y,b_y)*compute_density(q_z,p_z,k_type_z,b_z);
    } else{
        value = pow(p_x, m)*compute_density(q_y,p_y,k_type_y,b_y)*compute_density(q_z,p_z,k_type_z,b_z);
    }
    return value;
}



void SLAM(statistics& stat){
    if(stat.dim==2){
        // boundary_line_array
        if(stat.k_type_x==stat.k_type_y && stat.k_type_y==3){
            vector<vector<double>> SLWmatrix2_L_0(stat.n_x, vector<double>(stat.n_y, 0.0));
            vector<vector<double>> SLWmatrix2_U_0(stat.n_x, vector<double>(stat.n_y, 0.0));
            
            for(int m = 0; m<stat.n ; m++){
                int l = ceil((stat.featureVector[m][0]-stat.x_L+stat.b_x)/stat.incr_x);
                l = min(l,stat.n_x);
                l=max(l,0);
                if(l!=0){
                    l=l-1;
                }
                for(int j = 0; j<stat.n_y; j++){
                    double L_y = stat.y_L + j*stat.incr_y - stat.b_y;
                    double U_y = stat.y_L + j*stat.incr_y + stat.b_y;
                    if(stat.featureVector[m][1]>=L_y && stat.featureVector[m][1]<U_y){
                        SLWmatrix2_L_0[l][j]= SLWmatrix2_L_0[l][j]+1/stat.b_x/stat.b_y;
                    }
                }
            }
            for(int m = 0; m<stat.n ; m++){
                int u = ceil((stat.featureVector[m][0]-stat.x_L-stat.b_x)/stat.incr_x);
                u = min(u,stat.n_x);
                u = max(u,0);
                if(u!=0){
                    u=u-1;
                }
                for(int j = 0; j<stat.n_y; j++){
                    double L_y = stat.y_L + j*stat.incr_y - stat.b_y;
                    double U_y = stat.y_L + j*stat.incr_y + stat.b_y;
                    if(stat.featureVector[m][1]>=L_y && stat.featureVector[m][1]<U_y){
                        SLWmatrix2_U_0[u][j]= SLWmatrix2_U_0[u][j]+1/stat.b_x/stat.b_y;
                    }
                }
            }
            for (int i = 1; i < stat.n_x; ++i) {
                for (int j = 0; j < stat.n_y; ++j) {
                    SLWmatrix2_L_0[i][j] += SLWmatrix2_L_0[i - 1][j];
                    SLWmatrix2_U_0[i][j] += SLWmatrix2_U_0[i - 1][j];
                }
            }
            // Construct others dimension array
            for (int x_index = 0; x_index < stat.n_x; x_index++){
                for(int y_index = 0; y_index < stat.n_y; y_index++){
                    stat.out_tensor_2[x_index][y_index] = SLWmatrix2_U_0[x_index][y_index]-SLWmatrix2_L_0[x_index][y_index];
                }
            }
        }
        if(stat.k_type_x==stat.k_type_y && stat.k_type_y==1){
            vector<vector<double>> SLWmatrix2_L_0(stat.n_x, vector<double>(stat.n_y, 0.0));
            vector<vector<double>> SLWmatrix2_L_1(stat.n_x, vector<double>(stat.n_y, 0.0));
            vector<vector<double>> SLWmatrix2_L_2(stat.n_x, vector<double>(stat.n_y, 0.0));
            vector<vector<double>> SLWmatrix2_U_0(stat.n_x, vector<double>(stat.n_y, 0.0));
            vector<vector<double>> SLWmatrix2_U_1(stat.n_x, vector<double>(stat.n_y, 0.0));
            vector<vector<double>> SLWmatrix2_U_2(stat.n_x, vector<double>(stat.n_y, 0.0));
            for(int m = 0; m<stat.n ; m++){
                int l = ceil((stat.featureVector[m][0]-stat.x_L+stat.b_x)/stat.incr_x);
                l = min(l,stat.n_x);
                l=max(l,0);
                if(l!=0){
                    l=l-1;
                }
                for(int j = 0; j<stat.n_y; j++){
                    double L_y = stat.y_L + j*stat.incr_y - stat.b_y;
                    double U_y = stat.y_L + j*stat.incr_y + stat.b_y;
                    double q_y = stat.y_L + j*stat.incr_y;
                    if(stat.featureVector[m][1]>=L_y && stat.featureVector[m][1]<U_y){
                        SLWmatrix2_L_0[l][j]= SLWmatrix2_L_0[l][j]+compute_SLWmatrix_2(0,stat.featureVector[m][0],stat.featureVector[m][1],q_y,stat.k_type_y,stat.b_y);
                        SLWmatrix2_L_1[l][j]= SLWmatrix2_L_1[l][j]+compute_SLWmatrix_2(1,stat.featureVector[m][0],stat.featureVector[m][1],q_y,stat.k_type_y,stat.b_y);
                        SLWmatrix2_L_2[l][j]= SLWmatrix2_L_2[l][j]+compute_SLWmatrix_2(2,stat.featureVector[m][0],stat.featureVector[m][1],q_y,stat.k_type_y,stat.b_y);
                    }
                }
            }
            for(int m = 0; m<stat.n ; m++){
                int u = ceil((stat.featureVector[m][0]-stat.x_L-stat.b_x)/stat.incr_x);
                u = min(u,stat.n_x);
                u = max(u,0);
                if(u!=0){
                    u=u-1;
                }
                for(int j = 0; j<stat.n_y; j++){
                    double L_y = stat.y_L + j*stat.incr_y - stat.b_y;
                    double U_y = stat.y_L + j*stat.incr_y + stat.b_y;
                    double q_y = stat.y_L + j*stat.incr_y;
                    if(stat.featureVector[m][1]>=L_y && stat.featureVector[m][1]<U_y){
                        SLWmatrix2_U_0[u][j]= SLWmatrix2_U_0[u][j]+compute_SLWmatrix_2(0,stat.featureVector[m][0],stat.featureVector[m][1],q_y,stat.k_type_y,stat.b_y);
                        SLWmatrix2_U_1[u][j]= SLWmatrix2_U_1[u][j]+compute_SLWmatrix_2(1,stat.featureVector[m][0],stat.featureVector[m][1],q_y,stat.k_type_y,stat.b_y);
                        SLWmatrix2_U_2[u][j]= SLWmatrix2_U_2[u][j]+compute_SLWmatrix_2(2,stat.featureVector[m][0],stat.featureVector[m][1],q_y,stat.k_type_y,stat.b_y);
                    }
                }
            }
            for (int i = 1; i < stat.n_x; ++i) {
                for (int j = 0; j < stat.n_y; ++j) {
                    SLWmatrix2_L_0[i][j] += SLWmatrix2_L_0[i - 1][j];
                    SLWmatrix2_L_1[i][j] += SLWmatrix2_L_1[i - 1][j];
                    SLWmatrix2_L_2[i][j] += SLWmatrix2_L_2[i - 1][j];
                    SLWmatrix2_U_0[i][j] += SLWmatrix2_U_0[i - 1][j];
                    SLWmatrix2_U_1[i][j] += SLWmatrix2_U_1[i - 1][j];
                    SLWmatrix2_U_2[i][j] += SLWmatrix2_U_2[i - 1][j];
                }
            }
            cout << "SLWmatrix2_L_0"<<endl;
            // Construct others dimension array
            for (int x_index = 0; x_index < stat.n_x; x_index++){
                stat.q[0] = stat.x_L + x_index * stat.incr_x;
                for(int y_index = 0; y_index < stat.n_y; y_index++){
                    stat.out_tensor_2[x_index][y_index] = (1-stat.q[0]*stat.q[0]/stat.b_x/stat.b_x)*(SLWmatrix2_U_0[x_index][y_index]-SLWmatrix2_L_0[x_index][y_index])+2*stat.q[0]/stat.b_x/stat.b_x*(SLWmatrix2_U_1[x_index][y_index]-SLWmatrix2_L_1[x_index][y_index])-1/stat.b_x/stat.b_x*(SLWmatrix2_U_2[x_index][y_index]-SLWmatrix2_L_2[x_index][y_index]);
                }
            }
            
        }


    
    }
    if (stat.featureVector != nullptr) {
        for (int i = 0; i < stat.n; i++) {
            delete[] stat.featureVector[i]; 
            stat.featureVector[i] = nullptr; 
        }
        delete[] stat.featureVector;

        stat.featureVector = nullptr;
    }
    
}