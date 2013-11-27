#ifndef ORIGENBUILDER_H_INCLUDED
#define ORIGENBUILDER_H_INCLUDED

#include<iostream>
#include<vector>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "structures.h"


isoInformation BuildIsotope(std::ifstream &input);
isoInformation DataReader(isoInformation test1, int type, vector<isoInformation> input_stream);
isoInformation FuelBuilder(vector<isoInformation> fuel_values, double u235_mass);

#endif // ORIGENBUILDER_H_INCLUDED
