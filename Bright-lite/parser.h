#ifndef PARSER_H_INCLUDED
#define PARSER_H_INCLUDED
#include "structures.h"


isoInformation ParseOriginFile(std::string file_location);
vector<double> FindNeutronProduction(std::string line);

#endif // PARSER_H_INCLUDED
