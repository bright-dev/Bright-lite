#ifndef STRUCTURES_H_INCLUDED
#define STRUCTURES_H_INCLUDED

#include <utility>



struct daughter {
    std::string name;
    std::vector<double> mass;
    double fraction;
};

struct isoInformation {
    std::string name;
    std::vector<double> neutron_prod;
    std::vector<double> neutron_dest;
    double struct_prod;
    double struct_dest;
    std::vector<double> k_inf;
    std::vector<double> BUd;
    std::vector<double> time;
    std::vector<daughter> iso_vector;
    double fraction;
};

#endif // STRUCTURES_H_INCLUDED
