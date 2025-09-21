#include "MASS_CR.h"
#include <iostream>
#include <vector>
#include <algorithm>


int findGridIndex(const vector<double>& boundary_array, double value) {
    auto it = lower_bound(boundary_array.begin(), boundary_array.end(), value);
    if (it == boundary_array.begin()) {
        return 0;
    } else if (it == boundary_array.end()) {
        return boundary_array.size()-1;
    } else {
        return distance(boundary_array.begin(), it);
    }
}

pair<int, int> findGridCell_2(const Point& point, const vector<double>& boundary_line_array_x, const vector<double>& boundary_line_array_y) {
    int x_index = findGridIndex(boundary_line_array_x, point.x);
    int y_index = findGridIndex(boundary_line_array_y, point.y);
    return {x_index, y_index};
}

tuple<int, int, int> findGridCell_3(const Point& point,const vector<double>& boundary_line_array_x,const vector<double>& boundary_line_array_y,const vector<double>& boundary_line_array_z) {
    int x_index = findGridIndex(boundary_line_array_x, point.x);
    int y_index = findGridIndex(boundary_line_array_y, point.y);
    int z_index = findGridIndex(boundary_line_array_z, point.z);
    return {x_index, y_index, z_index};
}

void MASS_CR(statistics& stat){
    if(stat.k_type_x==1||stat.k_type_x==3){
        vector<double> boundary_line_array_x;
        boundary_line_array_x.reserve(2*stat.n_x);
        double base_value_L_x = stat.x_L-stat.b_x;
        double base_value_U_x = stat.x_L+stat.b_x;
        for (int i = 0; i < stat.n_x; ++i) {
            boundary_line_array_x.push_back(base_value_L_x + i * stat.incr_x);
        }
        for (int i = 0; i < stat.n_x; ++i) {
            boundary_line_array_x.push_back(base_value_U_x + i * stat.incr_x);
        }
        sort(boundary_line_array_x.begin(), boundary_line_array_x.end());

        vector<double> boundary_line_array_y;
        boundary_line_array_y.reserve(2*stat.n_y);
        double base_value_L_y = stat.y_L-stat.b_y;
        double base_value_U_y = stat.y_L+stat.b_y;
        for (int i = 0; i < stat.n_y; ++i) {
            boundary_line_array_y.push_back(base_value_L_y + i * stat.incr_y);
        }
        for (int i = 0; i < stat.n_y; ++i) {
            boundary_line_array_y.push_back(base_value_U_y + i * stat.incr_y);
        }
        sort(boundary_line_array_y.begin(), boundary_line_array_y.end());

        if(stat.dim==2){
            if(stat.k_type_x==stat.k_type_y && stat.k_type_y==3){
                vector<vector<double>> init_PREFIX_2_n(2*stat.n_x, vector<double>(2*stat.n_y, 0.0));
                for(int i=0;i<stat.n;i++){
                    Point point = {stat.featureVector[i][0], stat.featureVector[i][1], 0};
                    auto grid_cell = findGridCell_2(point, boundary_line_array_x, boundary_line_array_y);
                    int x_index = get<0>(grid_cell);
                    int y_index = get<1>(grid_cell);
                    init_PREFIX_2_n[x_index][y_index] += 1;
                }
                for (int i = 0; i < 2*stat.n_x; ++i) {
                    for (int j = 0; j < 2*stat.n_y; ++j) {
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
                for (int x_index = 0; x_index < stat.n_x; x_index++) {
                    stat.q[0] = stat.x_L + x_index * stat.incr_x;
                    double q_x_L = stat.q[0] - stat.b_x;
                    double q_x_U = stat.q[0] + stat.b_x;
                    for (int y_index = 0; y_index < stat.n_y; y_index++) {
                        stat.q[1] = stat.y_L + y_index * stat.incr_y;
                        double q_y_L = stat.q[1] - stat.b_y;
                        double q_y_U = stat.q[1] + stat.b_y;
                        Point point = {q_x_U, q_y_U, 0};
                        auto grid_cell = findGridCell_2(point, boundary_line_array_x, boundary_line_array_y);
                        Point point_x = {q_x_L, q_y_U, 0};
                        auto grid_cell_x = findGridCell_2(point_x, boundary_line_array_x, boundary_line_array_y);
                        Point point_y = {q_x_U, q_y_L, 0};
                        auto grid_cell_y = findGridCell_2(point_y, boundary_line_array_x, boundary_line_array_y);
                        Point point_xy = {q_x_L, q_y_L, 0};
                        auto grid_cell_xy = findGridCell_2(point_xy, boundary_line_array_x, boundary_line_array_y);
                        double PREFIX_2_n = init_PREFIX_2_n[get<0>(grid_cell)][get<1>(grid_cell)] -
                                            init_PREFIX_2_n[get<0>(grid_cell_x)][get<1>(grid_cell_x)] -
                                            init_PREFIX_2_n[get<0>(grid_cell_y)][get<1>(grid_cell_y)] +
                                            init_PREFIX_2_n[get<0>(grid_cell_xy)][get<1>(grid_cell_xy)];
                        stat.out_tensor_2[x_index][y_index] = PREFIX_2_n / stat.b_x / stat.b_y ;
                    }
                }
            }
            if(stat.k_type_x==stat.k_type_y && stat.k_type_y==1){
                vector<vector<double>> init_PREFIX_2_n(2*stat.n_x, vector<double>(2*stat.n_y, 0.0));
                vector<vector<double>> init_PREFIX_2_px(2*stat.n_x, vector<double>(2*stat.n_y, 0.0));
                vector<vector<double>> init_PREFIX_2_py(2*stat.n_x, vector<double>(2*stat.n_y, 0.0));
                vector<vector<double>> init_PREFIX_2_px2(2*stat.n_x, vector<double>(2*stat.n_y, 0.0));
                vector<vector<double>> init_PREFIX_2_py2(2*stat.n_x, vector<double>(2*stat.n_y, 0.0));
                vector<vector<double>> init_PREFIX_2_pxy(2*stat.n_x, vector<double>(2*stat.n_y, 0.0));
                vector<vector<double>> init_PREFIX_2_pxy2(2*stat.n_x, vector<double>(2*stat.n_y, 0.0));
                vector<vector<double>> init_PREFIX_2_px2y(2*stat.n_x, vector<double>(2*stat.n_y, 0.0));
                vector<vector<double>> init_PREFIX_2_px2y2(2*stat.n_x, vector<double>(2*stat.n_y, 0.0));
                for(int i=0;i<stat.n;i++){
                    Point point = {stat.featureVector[i][0], stat.featureVector[i][1], 0};
                    auto grid_cell = findGridCell_2(point, boundary_line_array_x, boundary_line_array_y);
                    int x_index = get<0>(grid_cell);
                    int y_index = get<1>(grid_cell);
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
                for (int i = 0; i < 2*stat.n_x; ++i) {
                    for (int j = 0; j < 2*stat.n_y; ++j) {
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
                for (int x_index = 0; x_index < stat.n_x; x_index++) {
                    stat.q[0] = stat.x_L + x_index * stat.incr_x;
                    double q_x_L = stat.q[0] - stat.b_x;
                    double q_x_U = stat.q[0] + stat.b_x;
                    for (int y_index = 0; y_index < stat.n_y; y_index++) {
                        stat.q[1] = stat.y_L + y_index * stat.incr_y;
                        double q_y_L = stat.q[1] - stat.b_y;
                        double q_y_U = stat.q[1] + stat.b_y;
                        Point point = {q_x_U, q_y_U, 0};
                        auto grid_cell = findGridCell_2(point, boundary_line_array_x, boundary_line_array_y);
                        Point point_x = {q_x_L, q_y_U, 0};
                        auto grid_cell_x = findGridCell_2(point_x, boundary_line_array_x, boundary_line_array_y);
                        Point point_y = {q_x_U, q_y_L, 0};
                        auto grid_cell_y = findGridCell_2(point_y, boundary_line_array_x, boundary_line_array_y);
                        Point point_xy = {q_x_L, q_y_L, 0};
                        auto grid_cell_xy = findGridCell_2(point_xy, boundary_line_array_x, boundary_line_array_y);
                        double PREFIX_2_n = init_PREFIX_2_n[get<0>(grid_cell)][get<1>(grid_cell)] -
                                            init_PREFIX_2_n[get<0>(grid_cell_x)][get<1>(grid_cell_x)] -
                                            init_PREFIX_2_n[get<0>(grid_cell_y)][get<1>(grid_cell_y)] +
                                            init_PREFIX_2_n[get<0>(grid_cell_xy)][get<1>(grid_cell_xy)];
                        double PREFIX_2_px = init_PREFIX_2_px[get<0>(grid_cell)][get<1>(grid_cell)] -
                                             init_PREFIX_2_px[get<0>(grid_cell_x)][get<1>(grid_cell_x)] -
                                             init_PREFIX_2_px[get<0>(grid_cell_y)][get<1>(grid_cell_y)] +
                                             init_PREFIX_2_px[get<0>(grid_cell_xy)][get<1>(grid_cell_xy)];
                        double PREFIX_2_py = init_PREFIX_2_py[get<0>(grid_cell)][get<1>(grid_cell)] -
                                             init_PREFIX_2_py[get<0>(grid_cell_x)][get<1>(grid_cell_x)] -
                                             init_PREFIX_2_py[get<0>(grid_cell_y)][get<1>(grid_cell_y)] +
                                             init_PREFIX_2_py[get<0>(grid_cell_xy)][get<1>(grid_cell_xy)];
                        double PREFIX_2_px2 = init_PREFIX_2_px2[get<0>(grid_cell)][get<1>(grid_cell)] -
                                              init_PREFIX_2_px2[get<0>(grid_cell_x)][get<1>(grid_cell_x)] -
                                              init_PREFIX_2_px2[get<0>(grid_cell_y)][get<1>(grid_cell_y)] +
                                              init_PREFIX_2_px2[get<0>(grid_cell_xy)][get<1>(grid_cell_xy)];
                        double PREFIX_2_py2 = init_PREFIX_2_py2[get<0>(grid_cell)][get<1>(grid_cell)] -
                                              init_PREFIX_2_py2[get<0>(grid_cell_x)][get<1>(grid_cell_x)] -
                                              init_PREFIX_2_py2[get<0>(grid_cell_y)][get<1>(grid_cell_y)] +
                                              init_PREFIX_2_py2[get<0>(grid_cell_xy)][get<1>(grid_cell_xy)];
                        double PREFIX_2_pxy = init_PREFIX_2_pxy[get<0>(grid_cell)][get<1>(grid_cell)] -
                                              init_PREFIX_2_pxy[get<0>(grid_cell_x)][get<1>(grid_cell_x)] -
                                              init_PREFIX_2_pxy[get<0>(grid_cell_y)][get<1>(grid_cell_y)] +
                                              init_PREFIX_2_pxy[get<0>(grid_cell_xy)][get<1>(grid_cell_xy)];
                        double PREFIX_2_pxy2 = init_PREFIX_2_pxy2[get<0>(grid_cell)][get<1>(grid_cell)] -
                                               init_PREFIX_2_pxy2[get<0>(grid_cell_x)][get<1>(grid_cell_x)] -
                                               init_PREFIX_2_pxy2[get<0>(grid_cell_y)][get<1>(grid_cell_y)] +
                                               init_PREFIX_2_pxy2[get<0>(grid_cell_xy)][get<1>(grid_cell_xy)];
                        double PREFIX_2_px2y = init_PREFIX_2_px2y[get<0>(grid_cell)][get<1>(grid_cell)] -
                                               init_PREFIX_2_px2y[get<0>(grid_cell_x)][get<1>(grid_cell_x)] -
                                               init_PREFIX_2_px2y[get<0>(grid_cell_y)][get<1>(grid_cell_y)] +
                                               init_PREFIX_2_px2y[get<0>(grid_cell_xy)][get<1>(grid_cell_xy)];
                        double PREFIX_2_px2y2 = init_PREFIX_2_px2y2[get<0>(grid_cell)][get<1>(grid_cell)] -
                                                init_PREFIX_2_px2y2[get<0>(grid_cell_x)][get<1>(grid_cell_x)] -
                                                init_PREFIX_2_px2y2[get<0>(grid_cell_y)][get<1>(grid_cell_y)] +
                                                init_PREFIX_2_px2y2[get<0>(grid_cell_xy)][get<1>(grid_cell_xy)];
                        stat.out_tensor_2[x_index][y_index] =
                                PREFIX_2_n * (1 - stat.q[0] * stat.q[0] / stat.b_x / stat.b_x) *
                                (1 - stat.q[1] * stat.q[1] / stat.b_y / stat.b_y) +
                                2 * stat.q[0] / stat.b_x / stat.b_x * (1 - stat.q[1] * stat.q[1] / stat.b_y / stat.b_y) *
                                PREFIX_2_px +
                                2 * stat.q[1] / stat.b_y / stat.b_y * (1 - stat.q[0] * stat.q[0] / stat.b_x / stat.b_x) *
                                PREFIX_2_py +
                                (stat.q[1] * stat.q[1] / stat.b_y / stat.b_y - 1) / stat.b_x / stat.b_x * PREFIX_2_px2 +
                                (stat.q[0] * stat.q[0] / stat.b_x / stat.b_x - 1) / stat.b_y / stat.b_y * PREFIX_2_py2 +
                                4 * stat.q[0] * stat.q[1] / stat.b_x / stat.b_x / stat.b_y / stat.b_y * PREFIX_2_pxy -
                                2 * stat.q[0] / stat.b_x / stat.b_x / stat.b_y / stat.b_y * PREFIX_2_pxy2 -
                                2 * stat.q[1] / stat.b_x / stat.b_x / stat.b_y / stat.b_y * PREFIX_2_px2y +
                                PREFIX_2_px2y2 / stat.b_x / stat.b_x / stat.b_y / stat.b_y;

                    }
                }
            }
        }
    } else if(stat.k_type_x==2){
        vector<double> boundary_line_array_x;
        boundary_line_array_x.reserve(3*stat.n_x);
        double base_value_L_x = stat.x_L-stat.b_x;
        double base_value_U_x = stat.x_L+stat.b_x;
        for (int i = 0; i < stat.n_x; ++i) {
            boundary_line_array_x.push_back(base_value_L_x + i * stat.incr_x);
        }
        for (int i = 0; i < stat.n_x; ++i) {
            boundary_line_array_x.push_back(base_value_U_x + i * stat.incr_x);
        }
        for (int i = 0; i < stat.n_x; ++i) {
            boundary_line_array_x.push_back(stat.x_L + i * stat.incr_x);
        }
        sort(boundary_line_array_x.begin(), boundary_line_array_x.end());

        vector<double> boundary_line_array_y;
        boundary_line_array_y.reserve(3*stat.n_y);
        double base_value_L_y = stat.y_L-stat.b_y;
        double base_value_U_y = stat.y_L+stat.b_y;
        for (int i = 0; i < stat.n_y; ++i) {
            boundary_line_array_y.push_back(base_value_L_y + i * stat.incr_y);
        }
        for (int i = 0; i < stat.n_y; ++i) {
            boundary_line_array_y.push_back(base_value_U_y + i * stat.incr_y);
        }
        for (int i = 0; i < stat.n_y; ++i) {
            boundary_line_array_y.push_back(stat.y_L + i * stat.incr_y);
        }
        sort(boundary_line_array_y.begin(), boundary_line_array_y.end());

        vector<double> boundary_line_array_z;
        

        if(stat.dim==2) {
            if (stat.k_type_x == stat.k_type_y && stat.k_type_y == 2) {
                vector<vector<double>> init_PREFIX_2_n(3 * stat.n_x, vector<double>(3 * stat.n_y, 0.0));
                vector<vector<double>> init_PREFIX_2_px(3 * stat.n_x, vector<double>(3 * stat.n_y, 0.0));
                vector<vector<double>> init_PREFIX_2_py(3 * stat.n_x, vector<double>(3 * stat.n_y, 0.0));
                vector<vector<double>> init_PREFIX_2_pxy(3 * stat.n_x, vector<double>(3 * stat.n_y, 0.0));
                for (int i = 0; i < stat.n; i++) {
                    Point point = {stat.featureVector[i][0], stat.featureVector[i][1], 0};
                    auto grid_cell = findGridCell_2(point, boundary_line_array_x, boundary_line_array_y);
                    int x_index = get<0>(grid_cell);
                    int y_index = get<1>(grid_cell);
                    init_PREFIX_2_n[x_index][y_index] += 1;
                    init_PREFIX_2_px[x_index][y_index] += point.x;
                    init_PREFIX_2_py[x_index][y_index] += point.y;
                    init_PREFIX_2_pxy[x_index][y_index] += point.x * point.y;

                }
                for (int i = 0; i < 3 * stat.n_x; ++i) {
                    for (int j = 0; j < 3 * stat.n_y; ++j) {
                        if (i == 0 && j == 0) {
                            init_PREFIX_2_n[i][j] = 0;
                        } else if (i != 0 && j == 0) {
                            init_PREFIX_2_n[i][j] += init_PREFIX_2_n[i - 1][j];
                            init_PREFIX_2_px[i][j] += init_PREFIX_2_px[i - 1][j];
                            init_PREFIX_2_py[i][j] += init_PREFIX_2_py[i - 1][j];
                            init_PREFIX_2_pxy[i][j] += init_PREFIX_2_pxy[i - 1][j];
                        } else if (i == 0 && j != 0) {
                            init_PREFIX_2_n[i][j] += init_PREFIX_2_n[i][j - 1];
                            init_PREFIX_2_px[i][j] += init_PREFIX_2_px[i][j - 1];
                            init_PREFIX_2_py[i][j] += init_PREFIX_2_py[i][j - 1];
                            init_PREFIX_2_pxy[i][j] += init_PREFIX_2_pxy[i][j - 1];
                        } else {
                            init_PREFIX_2_n[i][j] += init_PREFIX_2_n[i - 1][j] + init_PREFIX_2_n[i][j - 1] -
                                                     init_PREFIX_2_n[i - 1][j - 1];
                            init_PREFIX_2_px[i][j] += init_PREFIX_2_px[i - 1][j] + init_PREFIX_2_px[i][j - 1] -
                                                      init_PREFIX_2_px[i - 1][j - 1];
                            init_PREFIX_2_py[i][j] += init_PREFIX_2_py[i - 1][j] + init_PREFIX_2_py[i][j - 1] -
                                                      init_PREFIX_2_py[i - 1][j - 1];
                            init_PREFIX_2_pxy[i][j] += init_PREFIX_2_pxy[i - 1][j] + init_PREFIX_2_pxy[i][j - 1] -
                                                       init_PREFIX_2_pxy[i - 1][j - 1];
                        }
                    }
                }
                for (int x_index = 0; x_index < stat.n_x; x_index++) {
                    stat.q[0] = stat.x_L + x_index * stat.incr_x;
                    double q_x_L = stat.q[0] - stat.b_x;
                    double q_x_U = stat.q[0] + stat.b_x;
                    for (int y_index = 0; y_index < stat.n_y; y_index++) {
                        stat.q[1] = stat.y_L + y_index * stat.incr_y;
                        double q_y_L = stat.q[1] - stat.b_y;
                        double q_y_U = stat.q[1] + stat.b_y;
                        Point point = {q_x_U, q_y_U, 0};
                        auto grid_cell = findGridCell_2(point, boundary_line_array_x, boundary_line_array_y);
                        Point point_x = {q_x_L, q_y_U, 0};
                        auto grid_cell_x = findGridCell_2(point_x, boundary_line_array_x, boundary_line_array_y);
                        Point point_y = {q_x_U, q_y_L, 0};
                        auto grid_cell_y = findGridCell_2(point_y, boundary_line_array_x, boundary_line_array_y);
                        Point point_xy = {q_x_L, q_y_L, 0};
                        auto grid_cell_xy = findGridCell_2(point_xy, boundary_line_array_x, boundary_line_array_y);
                        Point point_qx = {stat.q[0], q_y_U, 0};
                        auto grid_cell_qx = findGridCell_2(point_qx, boundary_line_array_x, boundary_line_array_y);
                        Point point_qxy = {stat.q[0], q_y_L, 0};
                        auto grid_cell_qxy = findGridCell_2(point_qxy, boundary_line_array_x, boundary_line_array_y);
                        Point point_qy = {q_x_U, stat.q[1], 0};
                        auto grid_cell_qy = findGridCell_2(point_qy, boundary_line_array_x, boundary_line_array_y);
                        Point point_qyx = {q_x_L, stat.q[1], 0};
                        auto grid_cell_qyx = findGridCell_2(point_qyx, boundary_line_array_x, boundary_line_array_y);
                        Point point_qxqy = {stat.q[0], stat.q[1], 0};
                        auto grid_cell_qxqy = findGridCell_2(point_qxqy, boundary_line_array_x, boundary_line_array_y);
                        double PREFIX_2_n_1 = init_PREFIX_2_n[get<0>(grid_cell)][get<1>(grid_cell)] -
                                            init_PREFIX_2_n[get<0>(grid_cell_qx)][get<1>(grid_cell_qx)] -
                                            init_PREFIX_2_n[get<0>(grid_cell_qy)][get<1>(grid_cell_qy)] +
                                            init_PREFIX_2_n[get<0>(grid_cell_qxqy)][get<1>(grid_cell_qxqy)];
                        double PREFIX_2_px_1 = init_PREFIX_2_px[get<0>(grid_cell)][get<1>(grid_cell)] -
                                             init_PREFIX_2_px[get<0>(grid_cell_qx)][get<1>(grid_cell_qx)] -
                                             init_PREFIX_2_px[get<0>(grid_cell_qy)][get<1>(grid_cell_qy)] +
                                             init_PREFIX_2_px[get<0>(grid_cell_qxqy)][get<1>(grid_cell_qxqy)];
                        double PREFIX_2_py_1 = init_PREFIX_2_py[get<0>(grid_cell)][get<1>(grid_cell)] -
                                             init_PREFIX_2_py[get<0>(grid_cell_qx)][get<1>(grid_cell_qx)] -
                                             init_PREFIX_2_py[get<0>(grid_cell_qy)][get<1>(grid_cell_qy)] +
                                             init_PREFIX_2_py[get<0>(grid_cell_qxqy)][get<1>(grid_cell_qxqy)];
                        double PREFIX_2_pxy_1 = init_PREFIX_2_pxy[get<0>(grid_cell)][get<1>(grid_cell)] -
                                              init_PREFIX_2_pxy[get<0>(grid_cell_qx)][get<1>(grid_cell_qx)] -
                                              init_PREFIX_2_pxy[get<0>(grid_cell_qy)][get<1>(grid_cell_qy)] +
                                              init_PREFIX_2_pxy[get<0>(grid_cell_qxqy)][get<1>(grid_cell_qxqy)];
                        double PREFIX_2_n_2 = init_PREFIX_2_n[get<0>(grid_cell_qx)][get<1>(grid_cell_qx)] -
                                              init_PREFIX_2_n[get<0>(grid_cell_x)][get<1>(grid_cell_x)] -
                                              init_PREFIX_2_n[get<0>(grid_cell_qxqy)][get<1>(grid_cell_qxqy)] +
                                              init_PREFIX_2_n[get<0>(grid_cell_qyx)][get<1>(grid_cell_qyx)];
                        double PREFIX_2_px_2 = init_PREFIX_2_px[get<0>(grid_cell_qx)][get<1>(grid_cell_qx)] -
                                               init_PREFIX_2_px[get<0>(grid_cell_x)][get<1>(grid_cell_x)] -
                                               init_PREFIX_2_px[get<0>(grid_cell_qxqy)][get<1>(grid_cell_qxqy)] +
                                               init_PREFIX_2_px[get<0>(grid_cell_qyx)][get<1>(grid_cell_qyx)];
                        double PREFIX_2_py_2 = init_PREFIX_2_py[get<0>(grid_cell_qx)][get<1>(grid_cell_qx)] -
                                               init_PREFIX_2_py[get<0>(grid_cell_x)][get<1>(grid_cell_x)] -
                                               init_PREFIX_2_py[get<0>(grid_cell_qxqy)][get<1>(grid_cell_qxqy)] +
                                               init_PREFIX_2_py[get<0>(grid_cell_qyx)][get<1>(grid_cell_qyx)];
                        double PREFIX_2_pxy_2 = init_PREFIX_2_pxy[get<0>(grid_cell_qx)][get<1>(grid_cell_qx)] -
                                                init_PREFIX_2_pxy[get<0>(grid_cell_x)][get<1>(grid_cell_x)] -
                                                init_PREFIX_2_pxy[get<0>(grid_cell_qxqy)][get<1>(grid_cell_qxqy)] +
                                                init_PREFIX_2_pxy[get<0>(grid_cell_qyx)][get<1>(grid_cell_qyx)];
                        double PREFIX_2_n_3 = init_PREFIX_2_n[get<0>(grid_cell_qxqy)][get<1>(grid_cell_qxqy)] -
                                              init_PREFIX_2_n[get<0>(grid_cell_qyx)][get<1>(grid_cell_qyx)] -
                                              init_PREFIX_2_n[get<0>(grid_cell_qxy)][get<1>(grid_cell_qxy)] +
                                              init_PREFIX_2_n[get<0>(grid_cell_xy)][get<1>(grid_cell_xy)];
                        double PREFIX_2_px_3 = init_PREFIX_2_px[get<0>(grid_cell_qxqy)][get<1>(grid_cell_qxqy)] -
                                               init_PREFIX_2_px[get<0>(grid_cell_qyx)][get<1>(grid_cell_qyx)] -
                                               init_PREFIX_2_px[get<0>(grid_cell_qxy)][get<1>(grid_cell_qxy)] +
                                               init_PREFIX_2_px[get<0>(grid_cell_xy)][get<1>(grid_cell_xy)];
                        double PREFIX_2_py_3 = init_PREFIX_2_py[get<0>(grid_cell_qxqy)][get<1>(grid_cell_qxqy)] -
                                               init_PREFIX_2_py[get<0>(grid_cell_qyx)][get<1>(grid_cell_qyx)] -
                                               init_PREFIX_2_py[get<0>(grid_cell_qxy)][get<1>(grid_cell_qxy)] +
                                               init_PREFIX_2_py[get<0>(grid_cell_xy)][get<1>(grid_cell_xy)];
                        double PREFIX_2_pxy_3 = init_PREFIX_2_pxy[get<0>(grid_cell_qxqy)][get<1>(grid_cell_qxqy)] -
                                                init_PREFIX_2_pxy[get<0>(grid_cell_qyx)][get<1>(grid_cell_qyx)] -
                                                init_PREFIX_2_pxy[get<0>(grid_cell_qxy)][get<1>(grid_cell_qxy)] +
                                                init_PREFIX_2_pxy[get<0>(grid_cell_xy)][get<1>(grid_cell_xy)];
                        double PREFIX_2_n_4 = init_PREFIX_2_n[get<0>(grid_cell_qy)][get<1>(grid_cell_qy)] -
                                              init_PREFIX_2_n[get<0>(grid_cell_qxqy)][get<1>(grid_cell_qxqy)] -
                                              init_PREFIX_2_n[get<0>(grid_cell_y)][get<1>(grid_cell_y)] +
                                              init_PREFIX_2_n[get<0>(grid_cell_qxy)][get<1>(grid_cell_qxy)];
                        double PREFIX_2_px_4 = init_PREFIX_2_px[get<0>(grid_cell_qy)][get<1>(grid_cell_qy)] -
                                               init_PREFIX_2_px[get<0>(grid_cell_qxqy)][get<1>(grid_cell_qxqy)] -
                                               init_PREFIX_2_px[get<0>(grid_cell_y)][get<1>(grid_cell_y)] +
                                               init_PREFIX_2_px[get<0>(grid_cell_qxy)][get<1>(grid_cell_qxy)];
                        double PREFIX_2_py_4 = init_PREFIX_2_py[get<0>(grid_cell_qy)][get<1>(grid_cell_qy)] -
                                               init_PREFIX_2_py[get<0>(grid_cell_qxqy)][get<1>(grid_cell_qxqy)] -
                                               init_PREFIX_2_py[get<0>(grid_cell_y)][get<1>(grid_cell_y)] +
                                               init_PREFIX_2_py[get<0>(grid_cell_qxy)][get<1>(grid_cell_qxy)];
                        double PREFIX_2_pxy_4 = init_PREFIX_2_pxy[get<0>(grid_cell_qy)][get<1>(grid_cell_qy)] -
                                                init_PREFIX_2_pxy[get<0>(grid_cell_qxqy)][get<1>(grid_cell_qxqy)] -
                                                init_PREFIX_2_pxy[get<0>(grid_cell_y)][get<1>(grid_cell_y)] +
                                                init_PREFIX_2_pxy[get<0>(grid_cell_qxy)][get<1>(grid_cell_qxy)];
                        stat.out_tensor_2[x_index][y_index] =
                                PREFIX_2_n_1 * (1+stat.q[0]/stat.b_x+stat.q[1]/stat.b_y+stat.q[0]/stat.b_x*stat.q[1]/stat.b_y) -
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

}