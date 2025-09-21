#pragma once
#ifndef MASS_OPT_H
#define MASS_OPT_H

#include "init_visual.h"
#include "MASS_CR.h"



void MASS_OPT(statistics& stat);
vector<double> create_array_k1(double x_L, double b, double d, int n_x, int k);
vector<double> create_array_k0(double x_L, double b, double d, int n_x, int k);
vector<double> create_array_c0(double x_L, double b, double d, int n_x, int k);
pair<int, double> find_k_and_c(double b, double d);
int get_index(double p, double L, double b, double d, int k, double c);


#endif