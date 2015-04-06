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

#include <ctime>
//#include <eigen3/Eigen/Eigen>

#include <sys/time.h>
#include <ctime>


#include <stdio.h>
#include <math.h>
//#include <cblas.h>

//#include <gsl/gsl_math.h>
//#include <gsl/gsl_eigen.h>

//#include "lapacke.h"

//#include <armadillo>

double intpol(double y0, double y1, double x0, double x1, double x);
batch_info BatchCollapse(batch_info &batch);
fuelBundle CoreCollapse(fuelBundle &fuel);
fuelBundle fast_region_collapse(fuelBundle &fuel);
fuelBundle BatchCollapse_old(fuelBundle &fuel);
std::map<int, double> tomass (int ti, double time, isoInformation &isoinfo);
fuelBundle phicalc_simple(fuelBundle &core);
fuelBundle phicalc_eqpow(fuelBundle &core);
fuelBundle phicalc_cylindrical(fuelBundle &core);
double nusigf_finder(batch_info &batch);
double siga_finder(batch_info &batch);
double kcalc(fuelBundle &core);
void burnupcalc(fuelBundle &core, int mode, int DA_mode, double delta);
void burnupcalc_CR(fuelBundle &core, int mode, int DA_mode, double delta);
double burnupcalc_BU(fuelBundle &core, int mode, int DA_mode, double delta);
fuelBundle DA_calc(fuelBundle &fuel);
fuelBundle lib_interpol(fuelBundle &input_fuel);
void mass_check(fuelBundle &fuel);
double SS_burnupcalc(fuelBundle &core, int mode, int DA_mode, double delta, int N, double ss_fluence, double target_burnup);
double SS_burnupcalc_depricated(fuelBundle &core, int mode, int DA_mode, double delta, int N, double ss_fluence);
double SS_burnupcalc_CR(fuelBundle &core, int mode, int DA_mode, double delta, int N, double ss_fluence, double target_burnup);
std::pair<double, std::pair<double, std::map<int, double> > > blending_calc(fuelBundle &fuel, double BU_end, int mode, int da_mode, double time_step);
double CR_finder(fuelBundle &core);
double CR_batch(fuelBundle &core, int i);
void print_library(std::string name, fuelBundle &core);
#endif // BURNUPCALC_H_INCLUDED

