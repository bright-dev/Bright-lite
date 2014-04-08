#ifndef ORIGENBUILDER_H_INCLUDED
#define ORIGENBUILDER_H_INCLUDED

#include<iostream>
#include<vector>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "structures.h"

isoInformation BuildIsotope(std::ifstream &input);
isoInformation DataReader(isoInformation test1, int type, std::vector<isoInformation> input_stream);

isoInformation BuildIsotope2(std::ifstream &input, isoInformation &iso);
std::vector<isoInformation> DataReader2(std::string type, std::vector<isoInformation> &input_stream);

isoInformation FuelBuilder(std::vector<isoInformation> fuel_values, double u235_mass);
std::vector<nonActinide> NonActinideReader(std::string file_name);

#endif // ORIGENBUILDER_H_INCLUDED
