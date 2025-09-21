#pragma once
#ifndef INIT_VISUAL_H
#define INIT_VISUAL_H

#include "Library.h"

const double inf = 999999999999;
const double small_epsilon = 0.01;
const double pi = 3.14159265358979323846;

//#define OBTAIN_BANDWIDTHS

struct statistics
{
	int dim; //Data Dimensionality: 2 
    int n; //number of points in datasets
	int n_x; //number of discrete region in the x-axis, e.g., 1280
	int n_y; //number of discrete region in the y-axis, e.g., 960
	double x_L, x_U; //lower and upper region of data points in x-axis
	double y_L, y_U; //lower and upper region of data points in y-axis
    double incr_x; //incremental value in the x-axis
	double incr_y; //incremental value in the y-axis
	int method; //chosen method 0: scan 1: SLAM 2:MASS_CR 3:MASS_OPT 4:RQS_kd 5:RQS_range
	int k_type_x; //kernel type 1: Epanechnikov kernel, 2: Triangular kernel, 3: Uniform kernel
	int k_type_y; //kernel type 1: Epanechnikov kernel, 2: Triangular kernel, 3: Uniform kernel
 	double b_x_ratio; //control the bandwidth of x kernel
	double b_y_ratio; //control the bandwidth of y kernel
	double b_x; //The bandwidth of x kernel
	double b_y; //The bandwidth of y kernel
	char*dataFileName;
	char*outputFileName;
    double density;



	double*q;
	double**featureVector; //feature vector of all data points (d = 0, 1, and 2 represent x, y, and timestamp, respectively)
    double**out_tensor_2; //output 2d-visualization
    double runtime;

	//prefix_nosort
	double qx_lower_bound_L;
	double qy_lower_bound_L;
	double qz_lower_bound_L;

    const int leafCapacity = 40;

	
};


void initTensor(statistics& stat);
void initStat(int argc, char**argv, statistics& stat);
void extract_FeatureVector(statistics& stat);
void updateRegion_default(statistics& stat);
void update_incr_value(statistics& stat);
void obtain_bandwidth_values(statistics& stat);
void init_visual(int argc, char**argv, statistics& stat);
void save_result_to_file(statistics& stat);

#endif