#pragma once
#ifndef MASS_CR_H
#define MASS_CR_H

#include "init_visual.h"


struct Point {
    double x;
    double y;
    double z;
};

void MASS_CR(statistics& stat);
int findGridIndex(const vector<double>& boundary_array, double value);
pair<int, int> findGridCell_2(const Point& point, const vector<double>& boundary_line_array_x, const vector<double>& boundary_line_array_y);
tuple<int, int, int> findGridCell_3(const Point& point,const vector<double>& boundary_line_array_x,const vector<double>& boundary_line_array_y,const vector<double>& boundary_line_array_z);




#endif