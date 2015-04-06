#ifndef ORIGENBUILDER_H_INCLUDED
#define ORIGENBUILDER_H_INCLUDED

#include<iostream>
#include<vector>
#include<fstream>
#include<sstream>
#include<algorithm>
#include <cstdio>
#include <ctime>
#include <boost/timer.hpp>

#include "structures.h"

isoInformation BuildIsotope(std::ifstream &input);
fuelBundle StructReader(fuelBundle &core);
isoInformation DataReader(isoInformation test1, int type, std::vector<isoInformation> input_stream);
isoInformation BuildIsotope2(std::ifstream &input, isoInformation &iso);
std::vector<isoInformation> DataReader2(std::string type, std::vector<isoInformation> &input_stream);
double flux_finder(std::string type);
isoInformation BurnupBuilder(std::vector<isoInformation> &fuel_values);
isoInformation FuelBuilder(std::vector<isoInformation> &fuel_values);
std::vector<nonActinide> NonActinideReader(std::string file_name);

#endif // ORIGENBUILDER_H_INCLUDED
