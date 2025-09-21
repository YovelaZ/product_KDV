#pragma once
#ifndef SLAM_H
#define SLAM_H

#include "init_visual.h"
#include "SCAN.h"

struct index_x_pair
{
	int index;
	double time;
	bool operator<(const index_x_pair& pair) { return time < pair.time; }
};


void SLAM(statistics& stat);
pair<int, double> find_k_and_c(double b, double d);
void sort_FeatureVector(statistics& stat);
double compute_SLWmatrix_2(int m, double p_x, double p_y, double q_y, int k_type_y, double b_y);
double compute_SLWmatrix_3(int m, double p_x, double p_y, double p_z, double q_y, double q_z, int k_type_y, int k_type_z, double b_y, double b_z);

#endif