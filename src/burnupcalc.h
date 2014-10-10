#ifndef BURNUPCALC_H_INCLUDED
#define BURNUPCALC_H_INCLUDED

#include <iostream>
#include <vector>
#include <regex>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <map>
#include "structures.h"
#include "origenBuilder.h"
#include "cyclus.h"
//#include <eigen3/Eigen/Eigen>

#include <stdio.h>
#include <math.h>
#include <cblas.h>

//#include "lapacke.h"

//#include <armadillo>

double intpol(double y0, double y1, double x0, double x1, double x);
fuelBundle regionCollapse(fuelBundle fuel);
std::map<int, double> tomass (int ti, double time, isoInformation isoinfo);
fuelBundle phicalc_simple(fuelBundle core);
fuelBundle phicalc_cylindrical(fuelBundle core);
double nusigf_finder(batch_info batch);
double siga_finder(batch_info batch);
double Siga_finder(batch_info batch);
double kcalc(fuelBundle core);
fuelBundle burnupcalc(fuelBundle core, int mode, int DA_mode, double tolerance);
fuelBundle DA_calc(fuelBundle fuel);
fuelBundle lib_interpol(fuelBundle input_fuel);
void mass_check(fuelBundle fuel);
double SS_burnupcalc(isoInformation fuel, int N, double delta, double PNL);

#endif // BURNUPCALC_H_INCLUDED

