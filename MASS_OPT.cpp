#include "MASS_OPT.h"



pair<int, double> find_k_and_c(double b, double d) {
    if (d <= 0) {
        cout << "Pixel size must be greater than zero." << endl;
        return make_pair(-1, -1);
    }

    double lower_bound = floor(b / (d / double(2)));
    double upper_bound = ceil((b + d / double(2)) / (d / double(2)));

    // find k
    for (int k = lower_bound; k <= upper_bound; ++k) {
        double c = b - (k * d / double(2));

        if (c >= -0.0001 && c < d / double(2)) {
            return make_pair(k, c);
        }
    }
    cout << "Can't find k and c" << endl;

    return make_pair(-1, -1);
}

double get_q_lower_bound_L(double L, double b, double d, int k, double c){
    double q_lower_bound_L;
    if(abs(c)<=0.0001){
        q_lower_bound_L = L - b + floor(k / double(2)) * d;
    } else if(k%2==0){
        if (k==0){
            q_lower_bound_L = L-b;
        }else{
            q_lower_bound_L = L-b+ceil(k / double(2))*d;
        }
    } else{
        q_lower_bound_L = L+b-ceil(k / double(2))*d;
    }
    return q_lower_bound_L;
}
int get_index(double p, double L, double b, double d, int k, double c){
    int index;
    if(abs(c)<=0.0001){
        double base_value = L - b + floor(k / double(2)) * d;
        pair<int, double> result_k_c = find_k_and_c(p-base_value, 2*d);
        int k_2 = result_k_c.first;
        double c_2 = result_k_c.second;
        if(c_2 <= 0.00001){
            index = k_2;
        } else {
            index = k_2+1;
        }
    } else if(k%2==0){
        pair<int, double> result_k_c;
        int k_2 ;
        double c_2 ;
        if (k==0){
            pair<int, double> result_k_c = find_k_and_c(p-L+b, 2*d);
            k_2 = result_k_c.first;
            c_2 = result_k_c.second;
        }else{
            pair<int, double> result_k_c = find_k_and_c(p-L+b-ceil(k / double(2))*d, 2*d);
            k_2 = result_k_c.first;
            c_2 = result_k_c.second;
        }

        if(c_2 <= 0.00001){
            index = 2*k_2;
        } else if(c_2<=2*b-k*d+0.00001){
            index = 2*k_2+1;
        } else {
            index = 2*k_2+2;
        }
    } else{
        pair<int, double> result_k_c = find_k_and_c(p-L-b+ceil(k / double(2))*d, 2*d);
        int k_2 = result_k_c.first;
        double c_2 = result_k_c.second;
        if(c_2 <= 0.00001){
            index = 2*k_2;
        }else if(c_2<=(k+1)*d-2*b+0.00001){
            index = 2*k_2+1;
        } else {
            index = 2*k_2+2;
        }
    }

    return index;
}

int get_index_T(double p, double L, double b, double d, int k, double c){
    int index;
    if(abs(c)<=0.0001){
        if(k%2==0){
            double base_value = L - b + floor(k / double(2)) * d;
            pair<int, double> result_k_c = find_k_and_c(p-base_value, 2*d);
            int k_2 = result_k_c.first;
            double c_2 = result_k_c.second;
            if(c_2 <= 0.00001){
                index = k_2;
            } else {
                index = k_2+1;
            }
        }else{
            double base_value = L - b + floor(k / double(2)) * d;
            pair<int, double> result_k_c = find_k_and_c(p-base_value, d);
            int k_2 = result_k_c.first;
            double c_2 = result_k_c.second;
            if(c_2 <= 0.00001){
                index = k_2;
            } else {
                index = k_2+1;
            }
        }

    } else if(k%2==0){
        pair<int, double> result_k_c;
        int k_2 ;
        double c_2 ;
        if (k==0){
            pair<int, double> result_k_c = find_k_and_c(p-L+b, 2*d);
            k_2 = result_k_c.first;
            c_2 = result_k_c.second;
        }else{
            pair<int, double> result_k_c = find_k_and_c(p-L+b-ceil(k / double(2))*d, 2*d);
            k_2 = result_k_c.first;
            c_2 = result_k_c.second;
        }

        if(c_2 <= 0.00001){
            index = 3*k_2;
        } else if(c_2<=(2*b-k*d)/2+0.00001){
            index = 3*k_2+1;
        }else if(c_2>(2*b-k*d)/2+0.00001 && c_2<=2*b-k*d+0.00001){
            index = 3*k_2+2;
        } else {
            index = 3*k_2+3;
        }
    } else{
        pair<int, double> result_k_c = find_k_and_c(p-L-b+ceil(k / double(2))*d, 2*d);
        int k_2 = result_k_c.first;
        double c_2 = result_k_c.second;
        if(c_2 <= 0.00001){
            index = 3*k_2;
        }else if(c_2<=((k+1)*d-2*b)/2+0.00001){
            index = 3*k_2+1;
        } else if(c_2>((k+1)*d-2*b)/2+0.00001 && c_2<=(k+1)*d-2*b+0.00001){
            index = 3*k_2+2;
        }else {
            index = 3*k_2+3;
        }
    }

    return index;
}

void MASS_OPT(statistics& stat){
    if(stat.dim==2){
        if(stat.k_type_x==1||stat.k_type_x==3){
            if(stat.k_type_x==stat.k_type_y && stat.k_type_y==3){
                vector<vector<double>> init_PREFIX_2_n(2*stat.n_x+4, vector<double>(2*stat.n_y+4, 0.0));
                int k_x;
                double c_x;
                pair<int, double> result_k_c_x = find_k_and_c(stat.b_x, stat.incr_x);
                k_x = result_k_c_x.first;
                c_x = result_k_c_x.second;
                int k_y;
                double c_y;
                pair<int, double> result_k_c_y = find_k_and_c(stat.b_y, stat.incr_y);
                k_y = result_k_c_y.first;
                c_y = result_k_c_y.second;
                stat.qx_lower_bound_L=get_q_lower_bound_L(stat.x_L,stat.b_x,stat.incr_x,k_x,c_x);
                stat.qy_lower_bound_L=get_q_lower_bound_L(stat.y_L,stat.b_y,stat.incr_y,k_y,c_y);
                for(int i=0;i<stat.n;i++){

                    Point point = {stat.featureVector[i][0], stat.featureVector[i][1], 0};

                    int x_index = get_index(point.x,stat.x_L,stat.b_x,stat.incr_x,k_x,c_x);
                    int y_index = get_index(point.y,stat.y_L,stat.b_y,stat.incr_y,k_y,c_y);

                    init_PREFIX_2_n[x_index][y_index] += 1;

                }
                for (int i = 0; i < 2*stat.n_x+4; ++i) {
                    for (int j = 0; j < 2*stat.n_y+4; ++j) {
                        if(i==0 && j==0){
                            init_PREFIX_2_n[i][j] = 0;
                        }else if (i!=0 && j==0){
                            init_PREFIX_2_n[i][j] += init_PREFIX_2_n[i-1][j];

                        }else if (i==0 && j!=0){
                            init_PREFIX_2_n[i][j] += init_PREFIX_2_n[i][j-1];

                        }else{
                            init_PREFIX_2_n[i][j] += init_PREFIX_2_n[i-1][j] + init_PREFIX_2_n[i][j-1] - init_PREFIX_2_n[i-1][j-1];
                        }
                    }
                }
                for (int x_index = 0; x_index < stat.n_x; x_index++){
                    stat.q[0] = stat.x_L + x_index * stat.incr_x;
                    double q_x_L = stat.q[0]-stat.b_x;
                    if(q_x_L<=stat.qx_lower_bound_L){
                        q_x_L=stat.qx_lower_bound_L;
                    }
                    double q_x_U = stat.q[0]+stat.b_x;
                    if(q_x_U>=stat.x_U){
                        q_x_U=stat.x_U;
                    }
                    for(int y_index = 0; y_index < stat.n_y; y_index++){
                        stat.q[1] = stat.y_L + y_index * stat.incr_y;
                        double q_y_L = stat.q[1]-stat.b_y;
                        if(q_y_L<=stat.qy_lower_bound_L){
                            q_y_L=stat.qy_lower_bound_L;
                        }
                        double q_y_U = stat.q[1]+stat.b_y;
                        if(q_y_U>=stat.y_U){
                            q_y_U=stat.y_U;
                        }
                        int index_x_l = get_index(q_x_L,stat.x_L,stat.b_x,stat.incr_x,k_x,c_x);
                        int index_y_l = get_index(q_y_L,stat.y_L,stat.b_y,stat.incr_y,k_y,c_y);
                        int index_x_u = get_index(q_x_U,stat.x_L,stat.b_x,stat.incr_x,k_x,c_x);
                        int index_y_u = get_index(q_y_U,stat.y_L,stat.b_y,stat.incr_y,k_y,c_y);

                        double PREFIX_2_n = init_PREFIX_2_n[index_x_u][index_y_u]-init_PREFIX_2_n[index_x_l][index_y_u]-init_PREFIX_2_n[index_x_u][index_y_l]+init_PREFIX_2_n[index_x_l][index_y_l];
                        stat.out_tensor_2[x_index][y_index] = PREFIX_2_n/stat.b_x/stat.b_y;
                    }
                }
            }
            if(stat.k_type_x==stat.k_type_y && stat.k_type_y==1){
                vector<vector<double>> init_PREFIX_2_n(2*stat.n_x+4, vector<double>(2*stat.n_y+4, 0.0));
                vector<vector<double>> init_PREFIX_2_px(2*stat.n_x+4, vector<double>(2*stat.n_y+4, 0.0));
                vector<vector<double>> init_PREFIX_2_py(2*stat.n_x+4, vector<double>(2*stat.n_y+4, 0.0));
                vector<vector<double>> init_PREFIX_2_px2(2*stat.n_x+4, vector<double>(2*stat.n_y+4, 0.0));
                vector<vector<double>> init_PREFIX_2_py2(2*stat.n_x+4, vector<double>(2*stat.n_y+4, 0.0));
                vector<vector<double>> init_PREFIX_2_pxy(2*stat.n_x+4, vector<double>(2*stat.n_y+4, 0.0));
                vector<vector<double>> init_PREFIX_2_pxy2(2*stat.n_x+4, vector<double>(2*stat.n_y+4, 0.0));
                vector<vector<double>> init_PREFIX_2_px2y(2*stat.n_x+4, vector<double>(2*stat.n_y+4, 0.0));
                vector<vector<double>> init_PREFIX_2_px2y2(2*stat.n_x+4, vector<double>(2*stat.n_y+4, 0.0));
                int k_x;
                double c_x;
                pair<int, double> result_k_c_x = find_k_and_c(stat.b_x, stat.incr_x);
                k_x = result_k_c_x.first;
                c_x = result_k_c_x.second;
                int k_y;
                double c_y;
                pair<int, double> result_k_c_y = find_k_and_c(stat.b_y, stat.incr_y);
                k_y = result_k_c_y.first;
                c_y = result_k_c_y.second;
                stat.qx_lower_bound_L=get_q_lower_bound_L(stat.x_L,stat.b_x,stat.incr_x,k_x,c_x);
                stat.qy_lower_bound_L=get_q_lower_bound_L(stat.y_L,stat.b_y,stat.incr_y,k_y,c_y);


                for(int i=0;i<stat.n;i++){

                    Point point = {stat.featureVector[i][0], stat.featureVector[i][1], 0};

                    int x_index = get_index(point.x,stat.x_L,stat.b_x,stat.incr_x,k_x,c_x);
                    int y_index = get_index(point.y,stat.y_L,stat.b_y,stat.incr_y,k_y,c_y);

                    init_PREFIX_2_n[x_index][y_index] += 1;
                    init_PREFIX_2_px[x_index][y_index] += point.x;
                    init_PREFIX_2_py[x_index][y_index] += point.y;
                    init_PREFIX_2_px2[x_index][y_index] += point.x*point.x;
                    init_PREFIX_2_py2[x_index][y_index] += point.y*point.y;
                    init_PREFIX_2_pxy[x_index][y_index] += point.x*point.y;
                    init_PREFIX_2_pxy2[x_index][y_index] += point.x*point.y*point.y;
                    init_PREFIX_2_px2y[x_index][y_index] += point.x*point.x*point.y;
                    init_PREFIX_2_px2y2[x_index][y_index] += point.x*point.x*point.y*point.y;
                }
                for (int i = 0; i < 2*stat.n_x+4; ++i) {
                    for (int j = 0; j < 2*stat.n_y+4; ++j) {
                        if(i==0 && j==0){
                            init_PREFIX_2_n[i][j] = 0;
                        }else if (i!=0 && j==0){
                            init_PREFIX_2_n[i][j] += init_PREFIX_2_n[i-1][j];
                            init_PREFIX_2_px[i][j] += init_PREFIX_2_px[i-1][j];
                            init_PREFIX_2_py[i][j] += init_PREFIX_2_py[i-1][j];
                            init_PREFIX_2_px2[i][j] += init_PREFIX_2_px2[i-1][j];
                            init_PREFIX_2_py2[i][j] += init_PREFIX_2_py2[i-1][j];
                            init_PREFIX_2_pxy[i][j] += init_PREFIX_2_pxy[i-1][j];
                            init_PREFIX_2_pxy2[i][j] += init_PREFIX_2_pxy2[i-1][j];
                            init_PREFIX_2_px2y[i][j] += init_PREFIX_2_px2y[i-1][j];
                            init_PREFIX_2_px2y2[i][j] += init_PREFIX_2_px2y2[i-1][j];
                        }else if (i==0 && j!=0){
                            init_PREFIX_2_n[i][j] += init_PREFIX_2_n[i][j-1];
                            init_PREFIX_2_px[i][j] += init_PREFIX_2_px[i][j-1];
                            init_PREFIX_2_py[i][j] += init_PREFIX_2_py[i][j-1];
                            init_PREFIX_2_px2[i][j] += init_PREFIX_2_px2[i][j-1];
                            init_PREFIX_2_py2[i][j] += init_PREFIX_2_py2[i][j-1];
                            init_PREFIX_2_pxy[i][j] += init_PREFIX_2_pxy[i][j-1];
                            init_PREFIX_2_pxy2[i][j] += init_PREFIX_2_pxy2[i][j-1];
                            init_PREFIX_2_px2y[i][j] += init_PREFIX_2_px2y[i][j-1];
                            init_PREFIX_2_px2y2[i][j] += init_PREFIX_2_px2y2[i][j-1];
                        }else{
                            init_PREFIX_2_n[i][j] += init_PREFIX_2_n[i-1][j] + init_PREFIX_2_n[i][j-1] - init_PREFIX_2_n[i-1][j-1];
                            init_PREFIX_2_px[i][j] += init_PREFIX_2_px[i-1][j] + init_PREFIX_2_px[i][j-1] - init_PREFIX_2_px[i-1][j-1];
                            init_PREFIX_2_py[i][j] += init_PREFIX_2_py[i-1][j] + init_PREFIX_2_py[i][j-1] - init_PREFIX_2_py[i-1][j-1];
                            init_PREFIX_2_px2[i][j] += init_PREFIX_2_px2[i-1][j] + init_PREFIX_2_px2[i][j-1] - init_PREFIX_2_px2[i-1][j-1];
                            init_PREFIX_2_py2[i][j] += init_PREFIX_2_py2[i-1][j] + init_PREFIX_2_py2[i][j-1] - init_PREFIX_2_py2[i-1][j-1];
                            init_PREFIX_2_pxy[i][j] += init_PREFIX_2_pxy[i-1][j] + init_PREFIX_2_pxy[i][j-1] - init_PREFIX_2_pxy[i-1][j-1];
                            init_PREFIX_2_pxy2[i][j] += init_PREFIX_2_pxy2[i-1][j] + init_PREFIX_2_pxy2[i][j-1] - init_PREFIX_2_pxy2[i-1][j-1];
                            init_PREFIX_2_px2y[i][j] += init_PREFIX_2_px2y[i-1][j] + init_PREFIX_2_px2y[i][j-1] - init_PREFIX_2_px2y[i-1][j-1];
                            init_PREFIX_2_px2y2[i][j] += init_PREFIX_2_px2y2[i-1][j] + init_PREFIX_2_px2y2[i][j-1] - init_PREFIX_2_px2y2[i-1][j-1];
                        }
                    }
                }
                for (int x_index = 0; x_index < stat.n_x; x_index++){
                    stat.q[0] = stat.x_L + x_index * stat.incr_x;
                    double q_x_L = stat.q[0]-stat.b_x;
                    if(q_x_L<=stat.qx_lower_bound_L){
                        q_x_L=stat.qx_lower_bound_L;
                    }
                    double q_x_U = stat.q[0]+stat.b_x;
                    if(q_x_U>=stat.x_U){
                        q_x_U=stat.x_U;
                    }
                    for(int y_index = 0; y_index < stat.n_y; y_index++){
                        stat.q[1] = stat.y_L + y_index * stat.incr_y;
                        double q_y_L = stat.q[1]-stat.b_y;
                        if(q_y_L<=stat.qy_lower_bound_L){
                            q_y_L=stat.qy_lower_bound_L;
                        }
                        double q_y_U = stat.q[1]+stat.b_y;
                        if(q_y_U>=stat.y_U){
                            q_y_U=stat.y_U;
                        }
                        int index_x_l = get_index(q_x_L,stat.x_L,stat.b_x,stat.incr_x,k_x,c_x);
                        int index_y_l = get_index(q_y_L,stat.y_L,stat.b_y,stat.incr_y,k_y,c_y);
                        int index_x_u = get_index(q_x_U,stat.x_L,stat.b_x,stat.incr_x,k_x,c_x);
                        int index_y_u = get_index(q_y_U,stat.y_L,stat.b_y,stat.incr_y,k_y,c_y);

                        double PREFIX_2_n = init_PREFIX_2_n[index_x_u][index_y_u]-init_PREFIX_2_n[index_x_l][index_y_u]-init_PREFIX_2_n[index_x_u][index_y_l]+init_PREFIX_2_n[index_x_l][index_y_l];
                        double PREFIX_2_px = init_PREFIX_2_px[index_x_u][index_y_u]-init_PREFIX_2_px[index_x_l][index_y_u]-init_PREFIX_2_px[index_x_u][index_y_l]+init_PREFIX_2_px[index_x_l][index_y_l];
                        double PREFIX_2_py = init_PREFIX_2_py[index_x_u][index_y_u]-init_PREFIX_2_py[index_x_l][index_y_u]-init_PREFIX_2_py[index_x_u][index_y_l]+init_PREFIX_2_py[index_x_l][index_y_l];
                        double PREFIX_2_px2 = init_PREFIX_2_px2[index_x_u][index_y_u]-init_PREFIX_2_px2[index_x_l][index_y_u]-init_PREFIX_2_px2[index_x_u][index_y_l]+init_PREFIX_2_px2[index_x_l][index_y_l];
                        double PREFIX_2_py2 = init_PREFIX_2_py2[index_x_u][index_y_u]-init_PREFIX_2_py2[index_x_l][index_y_u]-init_PREFIX_2_py2[index_x_u][index_y_l]+init_PREFIX_2_py2[index_x_l][index_y_l];
                        double PREFIX_2_pxy = init_PREFIX_2_pxy[index_x_u][index_y_u]-init_PREFIX_2_pxy[index_x_l][index_y_u]-init_PREFIX_2_pxy[index_x_u][index_y_l]+init_PREFIX_2_pxy[index_x_l][index_y_l];
                        double PREFIX_2_pxy2 = init_PREFIX_2_pxy2[index_x_u][index_y_u]-init_PREFIX_2_pxy2[index_x_l][index_y_u]-init_PREFIX_2_pxy2[index_x_u][index_y_l]+init_PREFIX_2_pxy2[index_x_l][index_y_l];
                        double PREFIX_2_px2y = init_PREFIX_2_px2y[index_x_u][index_y_u]-init_PREFIX_2_px2y[index_x_l][index_y_u]-init_PREFIX_2_px2y[index_x_u][index_y_l]+init_PREFIX_2_px2y[index_x_l][index_y_l];
                        double PREFIX_2_px2y2 = init_PREFIX_2_px2y2[index_x_u][index_y_u]-init_PREFIX_2_px2y2[index_x_l][index_y_u]-init_PREFIX_2_px2y2[index_x_u][index_y_l]+init_PREFIX_2_px2y2[index_x_l][index_y_l];
                        double value = PREFIX_2_n*(1-stat.q[0]*stat.q[0]/stat.b_x/stat.b_x)*(1-stat.q[1]*stat.q[1]/stat.b_y/stat.b_y)+2*stat.q[0]/stat.b_x/stat.b_x*(1-stat.q[1]*stat.q[1]/stat.b_y/stat.b_y)*PREFIX_2_px+2*stat.q[1]/stat.b_y/stat.b_y*(1-stat.q[0]*stat.q[0]/stat.b_x/stat.b_x)*PREFIX_2_py+(stat.q[1]*stat.q[1]/stat.b_y/stat.b_y-1)/stat.b_x/stat.b_x*PREFIX_2_px2+(stat.q[0]*stat.q[0]/stat.b_x/stat.b_x-1)/stat.b_y/stat.b_y*PREFIX_2_py2+4*stat.q[0]*stat.q[1]/stat.b_x/stat.b_x/stat.b_y/stat.b_y*PREFIX_2_pxy-2*stat.q[0]/stat.b_x/stat.b_x/stat.b_y/stat.b_y*PREFIX_2_pxy2-2*stat.q[1]/stat.b_x/stat.b_x/stat.b_y/stat.b_y*PREFIX_2_px2y+PREFIX_2_px2y2/stat.b_x/stat.b_x/stat.b_y/stat.b_y;
                        stat.out_tensor_2[x_index][y_index] = PREFIX_2_n*(1-stat.q[0]*stat.q[0]/stat.b_x/stat.b_x)*(1-stat.q[1]*stat.q[1]/stat.b_y/stat.b_y)+2*stat.q[0]/stat.b_x/stat.b_x*(1-stat.q[1]*stat.q[1]/stat.b_y/stat.b_y)*PREFIX_2_px+2*stat.q[1]/stat.b_y/stat.b_y*(1-stat.q[0]*stat.q[0]/stat.b_x/stat.b_x)*PREFIX_2_py+(stat.q[1]*stat.q[1]/stat.b_y/stat.b_y-1)/stat.b_x/stat.b_x*PREFIX_2_px2+(stat.q[0]*stat.q[0]/stat.b_x/stat.b_x-1)/stat.b_y/stat.b_y*PREFIX_2_py2+4*stat.q[0]*stat.q[1]/stat.b_x/stat.b_x/stat.b_y/stat.b_y*PREFIX_2_pxy-2*stat.q[0]/stat.b_x/stat.b_x/stat.b_y/stat.b_y*PREFIX_2_pxy2-2*stat.q[1]/stat.b_x/stat.b_x/stat.b_y/stat.b_y*PREFIX_2_px2y+PREFIX_2_px2y2/stat.b_x/stat.b_x/stat.b_y/stat.b_y;

                    }
                }
            }
        }else if(stat.k_type_x==2 && stat.k_type_x==stat.k_type_y){
            vector<vector<double>> init_PREFIX_2_n(3*stat.n_x+6, vector<double>(3*stat.n_y+6, 0.0));
            vector<vector<double>> init_PREFIX_2_px(3*stat.n_x+6, vector<double>(3*stat.n_y+6, 0.0));
            vector<vector<double>> init_PREFIX_2_py(3*stat.n_x+6, vector<double>(3*stat.n_y+6, 0.0));
            vector<vector<double>> init_PREFIX_2_pxy(3*stat.n_x+6, vector<double>(3*stat.n_y+6, 0.0));
            int k_x;
            double c_x;
            pair<int, double> result_k_c_x = find_k_and_c(stat.b_x, stat.incr_x);
            k_x = result_k_c_x.first;
            c_x = result_k_c_x.second;
            int k_y;
            double c_y;
            pair<int, double> result_k_c_y = find_k_and_c(stat.b_y, stat.incr_y);
            k_y = result_k_c_y.first;
            c_y = result_k_c_y.second;
            stat.qx_lower_bound_L=get_q_lower_bound_L(stat.x_L,stat.b_x,stat.incr_x,k_x,c_x);
            stat.qy_lower_bound_L=get_q_lower_bound_L(stat.y_L,stat.b_y,stat.incr_y,k_y,c_y);


            for(int i=0;i<stat.n;i++){

                Point point = {stat.featureVector[i][0], stat.featureVector[i][1], 0};

                int x_index = get_index_T(point.x,stat.x_L,stat.b_x,stat.incr_x,k_x,c_x);
                int y_index = get_index_T(point.y,stat.y_L,stat.b_y,stat.incr_y,k_y,c_y);

                init_PREFIX_2_n[x_index][y_index] += 1;
                init_PREFIX_2_px[x_index][y_index] += point.x;
                init_PREFIX_2_py[x_index][y_index] += point.y;
                init_PREFIX_2_pxy[x_index][y_index] += point.x*point.y;
            }
            for (int i = 0; i < 3*stat.n_x+6; ++i) {
                for (int j = 0; j < 3*stat.n_y+6; ++j) {
                    if(i==0 && j==0){
                        init_PREFIX_2_n[i][j] = 0;
                    }else if (i!=0 && j==0){
                        init_PREFIX_2_n[i][j] += init_PREFIX_2_n[i-1][j];
                        init_PREFIX_2_px[i][j] += init_PREFIX_2_px[i-1][j];
                        init_PREFIX_2_py[i][j] += init_PREFIX_2_py[i-1][j];
                        init_PREFIX_2_pxy[i][j] += init_PREFIX_2_pxy[i-1][j];
                    }else if (i==0 && j!=0){
                        init_PREFIX_2_n[i][j] += init_PREFIX_2_n[i][j-1];
                        init_PREFIX_2_px[i][j] += init_PREFIX_2_px[i][j-1];
                        init_PREFIX_2_py[i][j] += init_PREFIX_2_py[i][j-1];
                        init_PREFIX_2_pxy[i][j] += init_PREFIX_2_pxy[i][j-1];
                    }else{
                        init_PREFIX_2_n[i][j] += init_PREFIX_2_n[i-1][j] + init_PREFIX_2_n[i][j-1] - init_PREFIX_2_n[i-1][j-1];
                        init_PREFIX_2_px[i][j] += init_PREFIX_2_px[i-1][j] + init_PREFIX_2_px[i][j-1] - init_PREFIX_2_px[i-1][j-1];
                        init_PREFIX_2_py[i][j] += init_PREFIX_2_py[i-1][j] + init_PREFIX_2_py[i][j-1] - init_PREFIX_2_py[i-1][j-1];
                        init_PREFIX_2_pxy[i][j] += init_PREFIX_2_pxy[i-1][j] + init_PREFIX_2_pxy[i][j-1] - init_PREFIX_2_pxy[i-1][j-1];
                    }
                }
            }
            for (int x_index = 0; x_index < stat.n_x; x_index++){
                stat.q[0] = stat.x_L + x_index * stat.incr_x;
                double q_x_L = stat.q[0]-stat.b_x;
                if(q_x_L<=stat.qx_lower_bound_L){
                    q_x_L=stat.qx_lower_bound_L;
                }
                double q_x_U = stat.q[0]+stat.b_x;
                if(q_x_U>=stat.x_U){
                    q_x_U=stat.x_U;
                }
                for(int y_index = 0; y_index < stat.n_y; y_index++){
                    stat.q[1] = stat.y_L + y_index * stat.incr_y;
                    double q_y_L = stat.q[1]-stat.b_y;
                    if(q_y_L<=stat.qy_lower_bound_L){
                        q_y_L=stat.qy_lower_bound_L;
                    }
                    double q_y_U = stat.q[1]+stat.b_y;
                    if(q_y_U>=stat.y_U){
                        q_y_U=stat.y_U;
                    }
                    int index_x_l = get_index_T(q_x_L,stat.x_L,stat.b_x,stat.incr_x,k_x,c_x);
                    int index_y_l = get_index_T(q_y_L,stat.y_L,stat.b_y,stat.incr_y,k_y,c_y);
                    int index_x_u = get_index_T(q_x_U,stat.x_L,stat.b_x,stat.incr_x,k_x,c_x);
                    int index_y_u = get_index_T(q_y_U,stat.y_L,stat.b_y,stat.incr_y,k_y,c_y);
                    int index_x_q = get_index_T(stat.q[0],stat.x_L,stat.b_x,stat.incr_x,k_x,c_x);
                    int index_y_q = get_index_T(stat.q[1],stat.y_L,stat.b_y,stat.incr_y,k_y,c_y);

                    double PREFIX_2_n_1 = init_PREFIX_2_n[index_x_u][index_y_u]-init_PREFIX_2_n[index_x_q][index_y_u]-init_PREFIX_2_n[index_x_u][index_y_q]+init_PREFIX_2_n[index_x_q][index_y_q];
                    double PREFIX_2_px_1 = init_PREFIX_2_px[index_x_u][index_y_u]-init_PREFIX_2_px[index_x_q][index_y_u]-init_PREFIX_2_px[index_x_u][index_y_q]+init_PREFIX_2_px[index_x_q][index_y_q];
                    double PREFIX_2_py_1 = init_PREFIX_2_py[index_x_u][index_y_u]-init_PREFIX_2_py[index_x_q][index_y_u]-init_PREFIX_2_py[index_x_u][index_y_q]+init_PREFIX_2_py[index_x_q][index_y_q];
                    double PREFIX_2_pxy_1 = init_PREFIX_2_pxy[index_x_u][index_y_u]-init_PREFIX_2_pxy[index_x_q][index_y_u]-init_PREFIX_2_pxy[index_x_u][index_y_q]+init_PREFIX_2_pxy[index_x_q][index_y_q];

                    double PREFIX_2_n_2 = init_PREFIX_2_n[index_x_q][index_y_u]-init_PREFIX_2_n[index_x_l][index_y_u]-init_PREFIX_2_n[index_x_q][index_y_q]+init_PREFIX_2_n[index_x_l][index_y_q];
                    double PREFIX_2_px_2 = init_PREFIX_2_px[index_x_q][index_y_u]-init_PREFIX_2_px[index_x_l][index_y_u]-init_PREFIX_2_px[index_x_q][index_y_q]+init_PREFIX_2_px[index_x_l][index_y_q];
                    double PREFIX_2_py_2 = init_PREFIX_2_py[index_x_q][index_y_u]-init_PREFIX_2_py[index_x_l][index_y_u]-init_PREFIX_2_py[index_x_q][index_y_q]+init_PREFIX_2_py[index_x_l][index_y_q];
                    double PREFIX_2_pxy_2 = init_PREFIX_2_pxy[index_x_q][index_y_u]-init_PREFIX_2_pxy[index_x_l][index_y_u]-init_PREFIX_2_pxy[index_x_q][index_y_q]+init_PREFIX_2_pxy[index_x_l][index_y_q];

                    double PREFIX_2_n_3 = init_PREFIX_2_n[index_x_q][index_y_q]-init_PREFIX_2_n[index_x_l][index_y_q]-init_PREFIX_2_n[index_x_q][index_y_l]+init_PREFIX_2_n[index_x_l][index_y_l];
                    double PREFIX_2_px_3 = init_PREFIX_2_px[index_x_q][index_y_q]-init_PREFIX_2_px[index_x_l][index_y_q]-init_PREFIX_2_px[index_x_q][index_y_l]+init_PREFIX_2_px[index_x_l][index_y_l];
                    double PREFIX_2_py_3 = init_PREFIX_2_py[index_x_q][index_y_q]-init_PREFIX_2_py[index_x_l][index_y_q]-init_PREFIX_2_py[index_x_q][index_y_l]+init_PREFIX_2_py[index_x_l][index_y_l];
                    double PREFIX_2_pxy_3 = init_PREFIX_2_pxy[index_x_q][index_y_q]-init_PREFIX_2_pxy[index_x_l][index_y_q]-init_PREFIX_2_pxy[index_x_q][index_y_l]+init_PREFIX_2_pxy[index_x_l][index_y_l];

                    double PREFIX_2_n_4 = init_PREFIX_2_n[index_x_u][index_y_q]-init_PREFIX_2_n[index_x_q][index_y_q]-init_PREFIX_2_n[index_x_u][index_y_l]+init_PREFIX_2_n[index_x_q][index_y_l];
                    double PREFIX_2_px_4 = init_PREFIX_2_px[index_x_u][index_y_q]-init_PREFIX_2_px[index_x_q][index_y_q]-init_PREFIX_2_px[index_x_u][index_y_l]+init_PREFIX_2_px[index_x_q][index_y_l];
                    double PREFIX_2_py_4 = init_PREFIX_2_py[index_x_u][index_y_q]-init_PREFIX_2_py[index_x_q][index_y_q]-init_PREFIX_2_py[index_x_u][index_y_l]+init_PREFIX_2_py[index_x_q][index_y_l];
                    double PREFIX_2_pxy_4 = init_PREFIX_2_pxy[index_x_u][index_y_q]-init_PREFIX_2_pxy[index_x_q][index_y_q]-init_PREFIX_2_pxy[index_x_u][index_y_l]+init_PREFIX_2_pxy[index_x_q][index_y_l];

                    stat.out_tensor_2[x_index][y_index] = PREFIX_2_n_1 * (1+stat.q[0]/stat.b_x+stat.q[1]/stat.b_y+stat.q[0]/stat.b_x*stat.q[1]/stat.b_y) -
                                                          PREFIX_2_px_1 * (1/stat.b_x+stat.q[1]/stat.b_x/stat.b_y) -
                                                          PREFIX_2_py_1 * (1/stat.b_y+stat.q[0]/stat.b_x/stat.b_y) +
                                                          PREFIX_2_pxy_1 / stat.b_x / stat.b_y +
                                                          PREFIX_2_n_2 * (1-stat.q[0]/stat.b_x+stat.q[1]/stat.b_y-stat.q[0]/stat.b_x*stat.q[1]/stat.b_y) +
                                                          PREFIX_2_px_2 * (1/stat.b_x+stat.q[1]/stat.b_x/stat.b_y) +
                                                          PREFIX_2_py_2 * (-1/stat.b_y+stat.q[0]/stat.b_x/stat.b_y) -
                                                          PREFIX_2_pxy_2 / stat.b_x / stat.b_y +
                                                          PREFIX_2_n_3 * (1-stat.q[0]/stat.b_x-stat.q[1]/stat.b_y+stat.q[0]/stat.b_x*stat.q[1]/stat.b_y) +
                                                          PREFIX_2_px_3 * (1/stat.b_x-stat.q[1]/stat.b_x/stat.b_y) +
                                                          PREFIX_2_py_3 * (1/stat.b_y-stat.q[0]/stat.b_x/stat.b_y) +
                                                          PREFIX_2_pxy_3 / stat.b_x / stat.b_y +
                                                          PREFIX_2_n_4 * (1+stat.q[0]/stat.b_x-stat.q[1]/stat.b_y-stat.q[0]/stat.b_x*stat.q[1]/stat.b_y) +
                                                          PREFIX_2_px_4 * (-1/stat.b_x+stat.q[1]/stat.b_x/stat.b_y) +
                                                          PREFIX_2_py_4 * (1/stat.b_y+stat.q[0]/stat.b_x/stat.b_y) -
                                                          PREFIX_2_pxy_4 / stat.b_x / stat.b_y ;
                }
            }
        }


    }
}
