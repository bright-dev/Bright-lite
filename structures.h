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
    std::vector<double> k_inf;
    std::vector<double> BUd;
    std::vector<double> fluence;
    std::vector<daughter> iso_vector;
    int region;
    double fraction;
};

struct nonActinide {
    std::string name;
    double sng;
    double sn2n;
    double snp;
    double sngx;
    double sn2nx;
    double yyn;
    double total_prod;
    double total_dest;
};

#endif // STRUCTURES_H_INCLUDED
