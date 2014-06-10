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
    std::string type;
    bool blending;
    //double fraction;
    bool fuel;
    double batch_fluence;
    std::vector<double> fraction;
    double sigs;
    double siga;
    std::vector<double> neutron_prod;
    std::vector<double> neutron_dest;
    std::vector<double> k_inf;
    std::vector<double> BUd;
    std::vector<double> fluence;
    std::vector<daughter> iso_vector;
};

struct nonActinide {
    int name;
    double sng;
    double scattering;
    double sn2n;
    double snp;
    double sngx;
    double sn2nx;
    double yyn;
    double total_prod;
    double total_dest;
};

struct interpol_pair {
    std::string metric;
    double value;
    double scaled_value;
};

struct fuelBundle {
    std::string name;
    std::string operation_type;
    int batch;
    bool libcheck;
    double pnl; //leakage
    double tres;
    double batch_fluence;
    double target_BUd;
    double fuel_radius;
    double moderator_radius;
    double moderator_sigs;
    double moderator_siga;
    std::vector<isoInformation> iso;
    std::vector<interpol_pair> interpol_pairs;
    std::vector<std::string> interpol_libs;
    std::vector<double> stream_fraction;
};




#endif // STRUCTURES_H_INCLUDED
