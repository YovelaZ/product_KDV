#include "alg_visual.h"
#include <unistd.h>   // for sysconf
#include <chrono>
#include <iostream>


void visual_Algorithm(statistics& stat)
{
    double run_time;

    //Visualization
    auto start_s = chrono::high_resolution_clock::now();

    if(stat.method == 0){ //scan
       SCAN(stat);
    }
    if(stat.method == 1){ //SLW
       SLAM(stat);
    }
    if(stat.method == 2){ //PREFIX
       MASS_CR(stat);

    }
    if(stat.method == 3){ //PREFIX_nosort
        MASS_OPT(stat);
    }
    if(stat.method == 4){ //kd_tree
        if (stat.dim==2){
            KDTree tree(stat);
            for (int x_index = 0; x_index < stat.n_x; x_index++){
                stat.q[0] = stat.x_L + x_index * stat.incr_x;
                double qx_LB = stat.q[0]-stat.b_x;
                double qx_UB = stat.q[0]+stat.b_x;
                for(int y_index = 0; y_index < stat.n_y; y_index++){
                    stat.q[1] = stat.y_L + y_index * stat.incr_y;
                    double qy_LB = stat.q[1]-stat.b_y;
                    double qy_UB = stat.q[1]+stat.b_y;
                    stat.density = 0;
                    auto inRange = tree.rangeQuery(qx_LB, qx_UB, qy_LB, qy_UB);

//                    for (int idx : inRange) {
//                        density = density + compute_density(stat.q[0],stat.featureVector[idx][0],stat.k_type_x,stat.b_x)*compute_density(stat.q[1],stat.featureVector[idx][1],stat.k_type_y,stat.b_y);
//                    }
                    stat.out_tensor_2[x_index][y_index] = stat.density;
                }
            }
        }
    }
    if(stat.method == 5){ //RANGE_tree
        if (stat.dim==2){
            cout << "Building RangeTree..." << endl;
            RangeTree  rt(stat);
            cout << "over Build RangeTree..." << endl;
            for (int x_index = 0; x_index < stat.n_x; x_index++){
                stat.q[0] = stat.x_L + x_index * stat.incr_x;
                double qx_LB = stat.q[0]-stat.b_x;
                if(qx_LB<stat.x_L){
                    qx_LB=stat.x_L;
                }
                double qx_UB = stat.q[0]+stat.b_x;
                if(qx_UB>stat.x_U){
                    qx_UB=stat.x_U;
                }
                for(int y_index = 0; y_index < stat.n_y; y_index++){
                    stat.q[1] = stat.y_L + y_index * stat.incr_y;
                    double qy_LB = stat.q[1]-stat.b_y;
                    if(qy_LB<stat.y_L){
                        qy_LB=stat.y_L;
                    }
                    double qy_UB = stat.q[1]+stat.b_y;
                    if(qy_UB>stat.y_U){
                        qy_UB=stat.y_U;
                    }
                    stat.density = 0;
                    rt.query(qx_LB, qx_UB, qy_LB, qy_UB);
//                    cout << "query RangeTree..." << endl;
                    stat.out_tensor_2[x_index][y_index] = stat.density;
//                    cout << x_index << endl;cout << y_index << endl;cout << stat.out_tensor_2[x_index][y_index] << endl;
                }
            }
        }
    }


    auto end_s = chrono::high_resolution_clock::now();
    run_time = (chrono::duration_cast<chrono::nanoseconds>(end_s - start_s).count()) / 1000000000.0;
    cout << "method " << stat.method << ":" << run_time << endl;
    stat.runtime = run_time;
}
