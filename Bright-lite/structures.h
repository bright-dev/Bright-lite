#ifndef STRUCTURES_H_INCLUDED
#define STRUCTURES_H_INCLUDED

#include <utility>

using namespace std;

struct daughter {
    string name;
    vector<double> mass;
    double fraction;
};

struct isoInformation {
    string name;
    vector<double> neutron_prod;
    vector<double> neutron_dest;
    vector<double> k_inf;
    vector<double> BUd;
    vector<double> time;
    vector<daughter> iso_vector;
    double fraction;
};

#endif // STRUCTURES_H_INCLUDED
