#ifndef PARSER_H_INCLUDED
#define PARSER_H_INCLUDED
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <ostream>
#include <algorithm>
#include "structures.h"

isoInformation ParseOriginFile(std::string file_location);
std::vector<double> FindNeutronProduction(std::string line);

#endif // PARSER_H_INCLUDED
