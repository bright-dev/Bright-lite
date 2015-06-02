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
    double batch_BU;
    double base_power;
    double base_mass;
    double base_flux;
    double fraction; //fraction of this isotope in the batch/stream
    double sigs;
    double siga;
    std::vector<double> neutron_prod;
    std::vector<double> neutron_dest;
    std::vector<double> k_inf;
    std::vector<double> BUd;
    std::vector<double> BU; //total BU
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

struct batch_info {
    std::string name;
    double batch_fluence; //fluence at the end of cycle, not used during burnup/composition calc
    std::vector<double> fraction; //blending fraction
    std::vector<isoInformation> iso;
    double batch_area; //[cm2] the total area of the batch for cylindrical flux calc
    isoInformation collapsed_iso;
    double BUg; //burnup guess, used during burnup/composition calc
    double Fg; //fluence guess, used during burnup/composition calc
    double rflux; //relative flux of batch
    double DA; //thermal disadvantage, phi_M/phi_F
    double discharge_BU; //the discharge burnup of the batch
    double delta_BU; // used to calculate the cycle length change in BU between cycles
    double discharge_CR; //the discharge conversion ratio
    std::map<int, double> comp; //current composition of batch at this batch_fluence
    double return_BU(){
        int ii;
        for(ii = 0; collapsed_iso.fluence[ii] <= Fg; ii++){}
        if(ii == 0){return 0;}
        return collapsed_iso.BU[ii-1] + ((collapsed_iso.BU[ii]-collapsed_iso.BU[ii-1])*
            (Fg - collapsed_iso.fluence[ii-1])/(collapsed_iso.fluence[ii] - collapsed_iso.fluence[ii-1]));
    }
};

struct fuelBundle {
    std::string name;
    std::string operation_type;
    int tot_batch;
    bool libcheck;
    double pnl; //leakage
    double tres; //residence time
    double base_flux;
    double base_power;
    double base_mass;
    double target_BU; //target discharge burnup, used for first guess in burnupcalc
    double fuel_area; //[cm2] the total area of the fuel for cylindrical flux calc
    double cylindrical_delta; //the increment used in cylindrical flux calc
    double mod_Sig_a; //Macroscopic absorption  cross section of the moderator
    double mod_Sig_tr; //Macroscopic transport cross section of the moderator
    double mod_Sig_f; //Macroscopic fission cross section of the moderator
    double mod_thickness; //radial thickness of the moderator [cm]
    double fuel_Sig_tr; //transport cs of the fuel
    double fuel_radius;
    double moderator_radius;
    double moderator_sigs;
    double moderator_siga;
    double disadv_a; //disadvantage calc fuel rad
    double disadv_b; //disadvantage calc mod rad
    double disadv_mod_siga; //disadvantage calc moderator Sig a
    double disadv_mod_sigs; //disadvantage calc mod Sig s
    double disadv_fuel_sigs; //disadvantage calc fuel Sig s
    double struct_prod; //neutron production rate of structural materials
    double struct_dest; //neutron destruction rate of structural materials
    double CR_upper;
    double CR_lower;
    double CR_target; //target CR value
    double CR; //current conversion ratio
    double SS_tolerance; //convergence tolerance for the code
    std::vector<int> CR_fissile; //list of fissile isotopes that are tracked for CR calc
    std::vector<batch_info> batch;
    std::vector<isoInformation> all_iso; //change to manifest
    std::vector<interpol_pair> interpol_pairs;
    std::vector<std::string> interpol_libs;
    std::vector<double> stream_fraction;
};

struct fuelInfo{
    double fluence;
    double burnup;
    std::map<int, double> burnup_info;
};


#endif // STRUCTURES_H_INCLUDED
