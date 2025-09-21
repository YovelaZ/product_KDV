#include "SCAN.h"
#include <cstdlib>

double compute_density(double q, double p, int k_type, double b){
    double value;
    if(k_type == 1){
        value = 1-(1/b)*(1/b)*abs(q-p)*abs(q-p);
    }
    if(k_type == 2){
        value = 1-(1/b)*abs(q-p);
    }
    if(k_type == 3){
        value = (1/b);
    }
    return value;
}

void SCAN(statistics& stat){
    if (stat.dim==2){
        for (int x_index = 0; x_index < stat.n_x; x_index++){
			stat.q[0] = stat.x_L + x_index * stat.incr_x;
            for(int y_index = 0; y_index < stat.n_y; y_index++){
                stat.q[1] = stat.y_L + y_index * stat.incr_y;
                double density = 0;
                for (int i = 0; i < stat.n; i++){
                    if(abs(stat.q[0]-stat.featureVector[i][0])<stat.b_x && abs(stat.q[1]-stat.featureVector[i][1])<stat.b_y)
                        density = density + compute_density(stat.q[0],stat.featureVector[i][0],stat.k_type_x,stat.b_x)*compute_density(stat.q[1],stat.featureVector[i][1],stat.k_type_y,stat.b_y);
                }
                stat.out_tensor_2[x_index][y_index] = density;
            }     
		}
    } else if (stat.dim==3){
        for (int x_index = 0; x_index < stat.n_x; x_index++){
			stat.q[0] = stat.x_L + x_index * stat.incr_x;
            for(int y_index = 0; y_index < stat.n_y; y_index++){
                stat.q[1] = stat.y_L + y_index * stat.incr_y;
                for(int z_index = 0; z_index < stat.n_z; z_index++){
                    stat.q[2] = stat.z_L + z_index * stat.incr_z;
                    double density = 0;
                    for (int i = 0; i < stat.n; i++){
                        if(abs(stat.q[0]-stat.featureVector[i][0])<stat.b_x && abs(stat.q[1]-stat.featureVector[i][1])<stat.b_y && abs(stat.q[2]-stat.featureVector[i][2])<stat.b_z)
                            density = density + compute_density(stat.q[0],stat.featureVector[i][0],stat.k_type_x,stat.b_x)*compute_density(stat.q[1],stat.featureVector[i][1],stat.k_type_y,stat.b_y)*compute_density(stat.q[2],stat.featureVector[i][2],stat.k_type_z,stat.b_z);
                    }
                    stat.out_tensor_3[x_index][y_index][z_index] = density;
                }
                
            }     
		}
    }
}