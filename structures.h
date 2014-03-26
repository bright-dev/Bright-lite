#ifndef STRUCTURES_H_INCLUDED
#define STRUCTURES_H_INCLUDED

#include <utility>



struct daughter {
    int name;
    std::vector<double> mass;
    double fraction;
};

struct isoInformation {
    int name; //changed this to int so it contains nucid
    int region;
    char type;
    double fraction;
    std::vector<double> neutron_prod;
    std::vector<double> neutron_dest;
    std::vector<double> k_inf;
    std::vector<double> BUd;
    std::vector<double> fluence;
    std::vector<daughter> iso_vector;
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

struct fuelBundle {
    std::string name;
    std::vector<isoInformation> iso;
};

struct fluxVector {
    std::string name;
    std::vector<double> flux;
};


#endif // STRUCTURES_H_INCLUDED
