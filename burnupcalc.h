#ifndef BURNUPCALC_H_INCLUDED
#define BURNUPCALC_H_INCLUDED

#include<iostream>
#include<vector>
#include<regex>
#include<iterator>
#include<fstream>
#include<algorithm>
#include<map>
#include "structures.h"
#include "origenBuilder.h"
#include "nucname.h"

double intpol(double y0, double y1, double x0, double x1, double x);
map<int, double> tomass (int ti, double time, isoInformation isoinfo);
pair<double, map<int, double> > burnupcalc(isoInformation tempone, int N, double tolerance);
double enrichcalc(double BU_end, int N, double tolerance, string type, vector<isoInformation> input_stream);

#endif // BURNUPCALC_H_INCLUDED

