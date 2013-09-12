#ifndef STRUCTURES_H_INCLUDED
#define STRUCTURES_H_INCLUDED

using namespace std;

struct daughter {
    string name;
    vector<double> mass;
};

struct isoInformation {
    vector<double> neutron_prod;
    vector<double> neutron_dest;
    vector<double> k_inf;
    vector<double> BUd;
    vector<double> time;
    vector<daughter> iso_vector;
};

#endif // STRUCTURES_H_INCLUDED
