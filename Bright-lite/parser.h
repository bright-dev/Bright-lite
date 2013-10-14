#ifndef PARSER_H_INCLUDED
#define PARSER_H_INCLUDED

#include "origenBuilder.h"
#include "burnupcalc.h"

isoInformation ParseOriginFile(string file_location);
vector<double> FindNeutronProduction(string line);

#endif // PARSER_H_INCLUDED
