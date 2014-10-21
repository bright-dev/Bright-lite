#include "burnupcalc.h"
#include <boost/math/special_functions/bessel.hpp>


using namespace std;

double intpol(double y0, double y1, double x0, double x1, double x) {
    // linear interpolation function
    double y = y0 + (y1 - y0)*(x - x0)/(x1 - x0);
    return y;
}
/*
fuelBundle NBuilder(fuelBundle fuel, vector<nonActinide> nona){
    //assumes the actinides are built in fuel, and that fuel has correct mass fractions
    int name;
    double total_prod = 0;
    double total_dest = 0;
    for(int i=0; i < fuel.iso.size(); i++){
        for(int j=0; j < nona.size(); j++){
            if(fuel.iso[i].name == nona[j].name){
                name = fuel.iso[i].name;
                name = name % 10000;
                name = name / 10;
                fuel.iso[i].neutron_prod.push_back(nona[j].total_prod*1000*0.602/name);
                total_prod += nona[j].total_prod*1000*0.602/name*fuel.iso[i].fraction[3];
                fuel.iso[i].neutron_dest.push_back(nona[j].total_dest*1000*0.602/name);
                total_dest += nona[j].total_dest*1000*0.602/name*fuel.iso[i].fraction[3];
            }
        }
    }
    cout << total_dest << "     " << total_prod << endl;
    int datasize;
    for(int i=0; i<fuel.iso.size(); i++){
        if(fuel.iso[i].type == "A"){
            datasize = fuel.iso[i].neutron_prod.size();
            break;
        }
    }

    for(int i =0; i<fuel.iso.size(); i++){
        if(fuel.iso[i].type == "N"){
            for(int j=0; j< datasize; j++){
                fuel.iso[i].neutron_prod.push_back(fuel.iso[i].neutron_prod[0]);
                fuel.iso[i].neutron_dest.push_back(fuel.iso[i].neutron_dest[0]);
            }
        }
    }
    return fuel;
}*/


fuelBundle regionCollapse(fuelBundle fuel){
///add micro region flux effects
//struct effects accounted here
    //cout << "Begin regionCollapse" << endl;
    for(int i = 0; i < fuel.batch.size(); i++){
        //for(int j = 0; j < fuel.batch[i].iso.size(); j ++){
        fuel.batch[i].collapsed_iso = FuelBuilder(fuel.batch[i].iso);
        //}
        //builds total BU from BUd
        fuel.batch[i].collapsed_iso.BU.push_back(fuel.batch[i].collapsed_iso.BUd[0]);

        for(int j = 1; j < fuel.batch[i].collapsed_iso.BUd.size(); j++){
        //cout << "    1.5tst" << i+1 << "  " << j << endl;
            fuel.batch[i].collapsed_iso.BU.push_back(fuel.batch[i].collapsed_iso.BU[j-1]+fuel.batch[i].collapsed_iso.BUd[j]);
        }
        //test to see if the prod/dest vectors are the same length
        if(fuel.batch[i].collapsed_iso.neutron_prod.size() != fuel.batch[i].collapsed_iso.neutron_dest.size()){
            cout << "Error. Neutron production/destruction rate vector length mismatch." << endl;
        }
        if(fuel.batch[i].collapsed_iso.neutron_prod[0] == 0 || fuel.batch[i].collapsed_iso.neutron_dest[0] == 0){
            fuel.batch[i].collapsed_iso.neutron_prod.erase(fuel.batch[i].collapsed_iso.neutron_prod.begin());
            fuel.batch[i].collapsed_iso.neutron_dest.erase(fuel.batch[i].collapsed_iso.neutron_dest.begin());
            fuel.batch[i].collapsed_iso.BUd.erase(fuel.batch[i].collapsed_iso.BUd.begin());
            fuel.batch[i].collapsed_iso.BU.erase(fuel.batch[i].collapsed_iso.BU.begin());
            fuel.batch[i].collapsed_iso.fluence.erase(fuel.batch[i].collapsed_iso.fluence.begin());
            for(int j = 0; j < fuel.batch[i].collapsed_iso.iso_vector.size(); j++){
                fuel.batch[i].collapsed_iso.iso_vector[j].mass.erase(fuel.batch[i].collapsed_iso.iso_vector[j].mass.begin());
            }
        }
    }
    return fuel;
};


map<int, double> tomass (int ti, double fluence, isoInformation isoinfo) {
    map<int, double> out = map<int, double>();
    double mass_i;
    int name_i;
    int nucid;
    for (int i = 0; i < isoinfo.iso_vector.size(); i++){
        name_i = isoinfo.iso_vector[i].name;
        nucid = pyne::nucname::zzaaam(name_i);
        mass_i = intpol(isoinfo.iso_vector[i].mass[ti-1],
                        isoinfo.iso_vector[i].mass[ti],
                        isoinfo.fluence[ti-1],
                        isoinfo.fluence[ti],
                        fluence);
        out[nucid] = mass_i;
    }
    return out;
}


fuelBundle phicalc_simple(fuelBundle core){
    //updates the rflux of each batch in core.batch
    //assumes the flux of a batch is proportional to the inverse neutron prod rate
    double maxphi = 0;

    //finds the inverse of neutron production at the batch_fluence
    //stores it in core.batch.rflux
    for(int i = 0; i < core.batch.size(); i++){
        int ii;
        for(ii = 0; core.batch[i].collapsed_iso.fluence[ii] < core.batch[i].Fg; ii++){}

        if(core.batch[i].collapsed_iso.fluence.back() < core.batch[i].Fg){
            cout << endl << "Maximum fluence error! Batch fluence exceeded max library fluence. (method2)" << endl;
            cout << "  Values on max fluence will be used. Do not trust results." << endl;
            ii = core.batch[i].collapsed_iso.fluence.size() - 1;
        }

        if(ii == 0){
            core.batch[i].rflux = 1/core.batch[i].collapsed_iso.neutron_prod[0];
            if(maxphi < core.batch[i].rflux){
                maxphi = core.batch[i].rflux;
            }

        } else{

            core.batch[i].rflux =
            1/intpol(core.batch[i].collapsed_iso.neutron_prod[ii-1],core.batch[i].collapsed_iso.neutron_prod[ii], core.batch[i].collapsed_iso.fluence[ii-1], core.batch[i].collapsed_iso.fluence[ii], core.batch[i].Fg+1);

        }

        if(maxphi < core.batch[i].rflux){
            maxphi = core.batch[i].rflux;
        }
    }

    //normalizes all the flux values
    for(int i = 0; i < core.batch.size(); i++){
        core.batch[i].rflux /= maxphi;
    }

    return core;
}

double nusigf_finder(batch_info batch){
    double Nusig_f;
    int ii;
    double prod = 0;
    int A = 235;

    //find the index to interpolate on
    for(ii = 0; batch.collapsed_iso.fluence[ii] < batch.Fg; ii++);
    if(batch.collapsed_iso.fluence.back() < batch.Fg){
        ii = batch.collapsed_iso.fluence.size() - 1;
    }

    prod = intpol(batch.collapsed_iso.neutron_prod[ii-1], batch.collapsed_iso.neutron_prod[ii], batch.collapsed_iso.fluence[ii-1], batch.collapsed_iso.fluence[ii], batch.Fg);

    //Nusig_f[barn] = prod rate / Avagadros number * mass * barn
    Nusig_f = prod/(6.0221413E+23)*A*1E24;

    return Nusig_f;
}

double siga_finder(batch_info batch){
//returns the microscopic cross section [barn] at fluence Fg
    double sig_a;
    int ii;
    double dest = 0;
    int A = 235;

    //find the index to interpolate on
    for(ii = 0; batch.collapsed_iso.fluence[ii] < batch.Fg; ii++);
    if(batch.collapsed_iso.fluence.back() < batch.Fg){
        ii = batch.collapsed_iso.fluence.size() - 1;
    }


    dest = intpol(batch.collapsed_iso.neutron_dest[ii-1], batch.collapsed_iso.neutron_dest[ii], batch.collapsed_iso.fluence[ii-1], batch.collapsed_iso.fluence[ii], batch.Fg);

    //sig_a[barn] = dest rate [cm2/kg] / Avagadros number * mass * barn
    sig_a = dest/(6.0221413E+23)*A*1E24/100;

    return sig_a;
}

double Siga_finder(batch_info batch){
//returns the macroscopic cross section [cm-1] at fluence Fg
    double Sig_a;
    int ii;
    double dest = 0;
    double rho = 0.01097; // [g/cm3]

    //find the index to interpolate on
    for(ii = 0; batch.collapsed_iso.fluence[ii] < batch.Fg; ii++);
    if(batch.collapsed_iso.fluence.back() < batch.Fg){
        ii = batch.collapsed_iso.fluence.size() - 1;
    }

    dest = intpol(batch.collapsed_iso.neutron_dest[ii-1], batch.collapsed_iso.neutron_dest[ii], batch.collapsed_iso.fluence[ii-1], batch.collapsed_iso.fluence[ii], batch.Fg);

    if(ii == 0){
        dest = batch.collapsed_iso.neutron_dest[0];
    }

    //Sig_a[cm-1] = dest rate * density
    Sig_a = dest*rho;

    return Sig_a;
}


////////////////////////////////////////////////////////////////////////////////////////////
fuelBundle phicalc_cylindrical(fuelBundle core){
cout << "cylindrical" << endl<<endl<<endl<<endl<<endl;

//the word 'region' is interchanged with the word 'batch' in comments, water layer is the outermost region

    int region = core.batch.size(); // number of regions (batches) not counting the outer region
    double delta = core.cylindrical_delta;
    double R[region+1]; //radial thickness of each region
    int N[region+1]; //number of mesh points in each region
    int NTotal; //total number of mesh points
    double dd2[region+1]; // D/delta^2
    double Sigma_a[region+1]; //mac. abs. cs of each region
    double NuSigma_f[region+1]; //nu sigma f of each region
    double Sigma_tr[region+1]; //mac. transport cs of each region
    double D[region+1]; //diff coef. for each region
    double LSquared[region+1];

    //set the radial thickness of each region
    R[0] = sqrt(core.fuel_area/region/3.141592);
    for(int i = 1; i < region; i++){
        R[i] = sqrt(core.fuel_area/region/3.141592*(i+1));
    }
    R[region] = R[region-1] + core.mod_thickness; //this is the moderator region

    //assign fuel cross sections
    for(int i = 0; i < region; i++){
        int ii;
        double prod, dest;
        //find the index to interpolate on
        for(ii = 0; core.batch[i].collapsed_iso.fluence[ii] < core.batch[i].Fg; ii++);
        if(core.batch[i].collapsed_iso.fluence.back() < core.batch[i].Fg){
            ii = core.batch[i].collapsed_iso.fluence.size() - 1;
        }

        prod = intpol(core.batch[i].collapsed_iso.neutron_prod[ii-1], core.batch[i].collapsed_iso.neutron_prod[ii], core.batch[i].collapsed_iso.fluence[ii-1], core.batch[i].collapsed_iso.fluence[ii], core.batch[i].Fg);
        NuSigma_f[i] = 0.0203;

        Sigma_a[i] = 1;

        Sigma_tr[i] = core.fuel_Sig_tr;
        D[i] = 1/(Sigma_tr[i]*3.);
        LSquared[i] = D[i]/Sigma_a[i];

    }

    //assign moderator cross sections
    Sigma_a[region] = core.mod_Sig_a;
    Sigma_tr[region] = core.mod_Sig_tr;
    D[region] = 1/(Sigma_tr[region]*3.);
    LSquared[region] = D[region]/Sigma_a[region];
    NuSigma_f[region] = core.mod_Sig_f;

    //populate dd2
    for(int i = 0; i < region+1; i++){
        dd2[i] = D[i]/delta/delta;
    }

    //populate N, number of mesh points in each region
    N[0] = R[0]/delta;
    NTotal = N[0];

    for(int i = 1; i < region+1; i++){
        N[i] = round((R[i] - R[i-1]) / delta);
        NTotal += N[i];
    }
    NTotal += 1;

    /*
    Eigen::MatrixXd A(NTotal, NTotal);
    Eigen::MatrixXd F(NTotal, NTotal);

    A.setZero();
    F.setZero();

    int jprev = 0;
    int j;

    //build the matrix A and F
    for(int i = 0; i <= region; i++){
        for(j = jprev+1; j < jprev + N[i]; j++){
            A(j,j-1) = -dd2[i]*(2.*(j)-1)/(2.*j);
            A(j,j) = 2*dd2[i] + Sigma_a[i];
            A(j,j+1) = -dd2[i]*(2.*(j)+1)/(2.*j);
            F(j,j) = NuSigma_f[i];
        }

        jprev = j;
        if(i != region){
            A(jprev,jprev-1) = D[i];
            A(jprev,jprev) = -D[i]-D[i+1];
            A(jprev,jprev+1) = D[i+1];
        }

    }

    A(0,0) = 8;//1;
    A(0,1) = 9;//-1;
    A(NTotal-1,NTotal-1) = 10;//1;

    // eigen.tuxfamily.org/dox/classEigen_1_1GeneralizedEigenSolver.html
    Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges;

    ges.compute(A, F);
    cout << ges.eigenvalues() << endl;
   */
   /*
    arma::mat  A = arma::zeros<arma::mat>(5,5);
    arma::mat F = arma::zeros<arma::mat>(5,5);

    for(int i = 0; i < 5; i++){
        A(i,i) = i+1;
        F(i,i) = i+2;
    }

    arma::cx_vec eigval;
    arma::cx_mat eigvec;

    arma::eig_pair(eigval, eigvec, A, F);

    cout << A << endl;
    cout << F << endl;

    cout << eigval << endl;
    cout << eigvec << endl;
    */


    for(int i = 0; i < core.batch.size(); i++){
        core.batch[i].rflux = 1;
    }

    return core;
}
/////////////////////////////////////////////////////////////////////////////////////////////

double kcalc(fuelBundle core){
    //uses the fluence values in core.batch.Fg to calculate the core criticality
    int N = core.batch.size();
    double prod_tot = 0;
    double dest_tot = 0;
    double pnl = core.pnl;

    //finds the index j where fluence is just under target, then interpolates
    //  to find the prod and dest values for each batch
    int j;
    for(int i = 0; i < N; i++){

        //find dicrete point j to interpolate on
        for(j = 0; core.batch[i].collapsed_iso.fluence[j] < core.batch[i].Fg; j++){}
        if(core.batch[i].collapsed_iso.fluence.back() < core.batch[i].Fg){
            //cout << endl << "Maximum fluence error! Batch fluence exceeded max library fluence. (kcalc)" << endl;
            //cout << "  Values on max fluence will be used. Do not trust results." << endl;
            j = core.batch[i].collapsed_iso.fluence.size() - 1;
        }

        //add the production rate of batch i to total production
        prod_tot += intpol(core.batch[i].collapsed_iso.neutron_prod[j-1], core.batch[i].collapsed_iso.neutron_prod[j], core.batch[i].collapsed_iso.fluence[j-1], core.batch[i].collapsed_iso.fluence[j], core.batch[i].Fg);

        //add structural material production of this batch, scaled up using disadvantage factor
        prod_tot += core.struct_prod * core.batch[i].DA;

        //add the destruction rate of batch i to total destruction
        dest_tot += intpol(core.batch[i].collapsed_iso.neutron_dest[j-1], core.batch[i].collapsed_iso.neutron_dest[j], core.batch[i].collapsed_iso.fluence[j-1], core.batch[i].collapsed_iso.fluence[j], core.batch[i].Fg);

        //add structural material destruction of this batch, scaled up using disadvantage factor
        dest_tot += core.struct_dest * core.batch[i].DA;

    }
    if(pnl < 0 || pnl > 1.5){
        cout << endl << "Error in core nonleakage! Assumed new value: 1.000" << endl;
        pnl = 1.000;
    }

    return prod_tot * pnl / dest_tot;
}

fuelBundle burnupcalc(fuelBundle core, int mode, int DA_mode, double delta) {
    //this function only uses the COLLAPSED_ISO of each BATCH in the structure CORE
    //all factors that contribute to a change in neutron prod/dest rates have to be factored
    //      before calling this function
    //cout << endl << "Burnupcalc" << endl;

    double BUg = 40; //burnup guess, gets updated if a guess in passed. unused
    int N = core.batch.size(); //number of batches
    double dt = delta*24*60*60; //days to [s]
    double kcore, kcore_prev;
    double y0, y1, x0, x1;
    double burnup = 0;

    //assign the batch_fluence to Fg
    for(int i = 0; i < N; i++){
        core.batch[i].Fg = core.batch[i].batch_fluence;
    }

    //turn zero fluences to one to avoid zero rates
    for(int i = 0; i < N; i++){
        if(core.batch[i].Fg == 0){core.batch[i].Fg = 1;}
    }

    for(int i = 0; i < core.batch[0].collapsed_iso.fluence.size(); i++){
        //cout << core.batch[0].collapsed_iso.fluence[i] << "  " << core.batch[0].collapsed_iso.BU[i] << "  " << core.batch[0].collapsed_iso.neutron_prod[i]/core.batch[0].collapsed_iso.neutron_dest[i] << endl;
    }
    //cout << "Batches: " << core.batch.size() << endl <<endl;
    for(int i = 0; i < core.batch.size(); i++){
        //cout << "fluence: " << core.batch[i].batch_fluence << endl;
    }

    kcore = 3.141592;
    kcore_prev = kcalc(core);


    //<--------------------------------------------------
    double cburnup[N];
    for(int i = 0; i < N; i++){cburnup[i] = 0;}
    //-------------------------------------------------->


    //more forward in time until kcore drops under 1
    while(kcore > 1){
        kcore_prev = kcore;
        //if(kcore != 3.141592){cout << "kcore: " << kcore << endl;}

        //find the normalized relative flux of each batch
        if(mode == 1){
            //simplest case, all batches get the same flux
            for(int i = 0; i < N; i++){
                core.batch[i].rflux = 1;
            }
        }else if(mode == 2){
            //inverse-production flux calculation
            core = phicalc_simple(core);
        }else if(mode == 3){
            core = phicalc_cylindrical(core);
        }else{
            cout << endl << "Error in mode input for batch-level flux calculation." << endl;
            return core;
        }

        //disadvantage calculation
        if(DA_mode == 1){
            core = DA_calc(core);
        }

        //update fluences
        for(int i = 0; i < N; i++){
            //cout << "  Added fluence: " << core.batch[i].rflux * core.base_flux * dt << endl;
            core.batch[i].Fg += core.batch[i].rflux * core.base_flux * dt;
        }

        //<--------------------------------------------------
        double added[N];
        double fnow[N];
        double bburnup[N];
        double totburnup = 0;

        for(int i = 0; i < N; i++){
            added[i] = core.batch[i].rflux * core.base_flux * dt;
            fnow[i] = core.batch[i].Fg;

            int ii;
            for(ii = 0; core.batch[i].collapsed_iso.fluence[ii] < fnow[i]-added[i]; ii++){}
            bburnup[i] = intpol(core.batch[i].collapsed_iso.BU[ii-1],core.batch[i].collapsed_iso.BU[ii], core.batch[i].collapsed_iso.fluence[ii-1], core.batch[i].collapsed_iso.fluence[ii], fnow[i]-added[i]);

            for(ii = 0; core.batch[i].collapsed_iso.fluence[ii] < fnow[i]; ii++){}

            bburnup[i] = intpol(core.batch[i].collapsed_iso.BU[ii-1],core.batch[i].collapsed_iso.BU[ii], core.batch[i].collapsed_iso.fluence[ii-1], core.batch[i].collapsed_iso.fluence[ii], fnow[i]) - bburnup[i];

            totburnup += bburnup[i];

            cburnup[i] += bburnup[i];
            //cout << "burnup this cycle of batch " << i+1 << ": " << bburnup[i] << endl;



        }
        //------------------------------------------------->

        kcore = kcalc(core);

        //cout << "k=" << kcore << "  ";
    }
    //cout << endl;

    //<--------------------------------------------------
    for(int i = 0; i < N; i++){
        //cout << "total burnup of batch " << i+1 << ": " << cburnup[i] << endl;
    }
    //------------------------------------------------->


    //update core fluences
    for(int i = 0; i < N; i++){
        //y0 is the fluence value before the last interation
        y0 = core.batch[i].Fg - core.batch[i].rflux * core.base_flux * dt;
        y1 = core.batch[i].Fg;
        //cout << y0 << " and " << y1 << endl; //<-------
        //cout << "  " << kcore_prev << " " << kcore << endl; //<-------
        core.batch[i].batch_fluence = intpol(y0, y1, kcore_prev, kcore, 1);
        //cout << "  fluence end of burnupcalc: " << core.batch[i].batch_fluence << endl;

    }

    //update current composition of batches
    for(int i = 0; i < N; i++){
        core.batch[i].comp.clear();
        int ii;
        for(ii = 0; core.batch[i].collapsed_iso.fluence[ii] < core.batch[i].batch_fluence; ii++){}
        if(core.batch[i].collapsed_iso.fluence.back() < core.batch[i].batch_fluence){
            cout << endl << "Maximum fluence error! Batch fluence exceeded max library fluence. (burnupcalc1)" << endl;
            cout << "  Values on max fluence will be used. Do not trust results." << endl;
            ii = core.batch[i].collapsed_iso.fluence.size() - 1;
        }

        for(int j = 0; j < core.batch[i].collapsed_iso.iso_vector.size(); j++){
            core.batch[i].comp[core.batch[i].collapsed_iso.iso_vector[j].name] =
            intpol(core.batch[i].collapsed_iso.iso_vector[j].mass[ii-1],core.batch[i].collapsed_iso.iso_vector[j].mass[ii], core.batch[i].collapsed_iso.fluence[ii-1], core.batch[i].collapsed_iso.fluence[ii], core.batch[i].batch_fluence)/1000;
        }
    }

    //the oldest batch is always index=0
    int ii;
    for(ii = 0; core.batch[0].collapsed_iso.fluence[ii] < core.batch[0].batch_fluence; ii++){}
    if(core.batch[0].collapsed_iso.fluence.back() < core.batch[0].batch_fluence){
        cout << endl << "Maximum fluence error! Batch fluence exceeded max library fluence. (burnupcalc2)" << endl;
        cout << "  Values on max fluence will be used. Do not trust results." << endl;
        ii = core.batch[0].collapsed_iso.fluence.size() - 1;
    }
    burnup = intpol(core.batch[0].collapsed_iso.BU[ii-1], core.batch[0].collapsed_iso.BU[ii], core.batch[0].collapsed_iso.fluence[ii-1], core.batch[0].collapsed_iso.fluence[ii], core.batch[0].batch_fluence);

    core.batch[0].discharge_BU = burnup;

    //cout << endl << "Discharge burnup: " << burnup << endl << endl;

    /************************output file*********************************/
    std::ofstream outfile;
    outfile.open("../output_cyclus_recent.txt", std::ios::app);

    outfile << "Discharge burnup: " << burnup;
    if(core.batch[0].collapsed_iso.fluence.back() < core.batch[0].batch_fluence){
        outfile << "\r\n\r\nBatch fluence exceeded max fluence in library! Do not trust results.\r\n     Discharge values extrapolated from last two discrete points.\r\n\r\n";
    }

    outfile << "\r\n Fluences at the end of this burnup cycle:";

    for(int i = 0; i < N; i++){
        outfile << "\r\n  Batch " << i+1 << ": " << core.batch[i].batch_fluence;
    }

    outfile << "\r\n    Discharge fuel fractions '[ISOTOPE]  [FRACTION]':";
    //depends on variable ii, which was found during burnup calculation
    for(std::map<int,double>::iterator it = core.batch[0].comp.begin(); it!=core.batch[0].comp.end(); ++it){
        outfile << "\r\n       " << it->first << "  " << it->second;
     }

    outfile << "\r\n";

    outfile.close();
    /************************End of output file**************************/

    return core;
}



fuelBundle DA_calc(fuelBundle fuel){
// calculates the flux within each batch in fuelBundle
// DA = phi_thermal_Mod / phi_thermal_Fuel

    double a = fuel.disadv_a; // radius of the fuel rod
    double b = fuel.disadv_b; // radius of the equivalent cell
    double Sig_sF = fuel.disadv_fuel_sigs; // macroscopic scatter CS of fuel
    double Sig_sM = fuel.disadv_mod_sigs; //macroscopic scatter CS of moderator
    double Sig_aM = fuel.disadv_mod_siga; // macroscopic abs. CS of moderator

    double L_F; // diffusion length of fuel
    double L_M; // diffusion length of moderator
    double Sig_aF; // macroscopic abs. CS of fuel
    double V_F; // volume of fuel
    double V_M; // volume of moderator
    double Sig_trF; // macroscopic transport CS of fuel
    double Sig_trM; // macroscopic transport CS of moderator
    double Sig_tF; // macroscopic total CS of fuel
    double Sig_tM; //macroscopic total CS of moderator
    double D_F; // diffusion coef. of fuel
    double D_M; // diffusion coef. of moderator
    double A_F; // A number of fuel
    double A_M; // A number of moderator
    double x, y, z; // calculated equivalent dimensions
    double F, E; // lattice functions
    double f; // flux of fuel divided by total flux(fuel+moderator)


    //moderator
    Sig_tM = Sig_aM + Sig_sM;
    A_F = 235;
    A_M = 18;
    Sig_trM = Sig_tM - 2/3/A_M*Sig_sM;
    D_M = 1 / (3 * Sig_trM);
    L_M = sqrt(D_M/Sig_aM);
    y = a/L_M;
    z = b/L_M;
    V_M = pow(b,2)*3.141592 - pow(a,2)*3.141592;
    V_F = pow(a,2)*3.141592;

    //cout << endl << "thermal disadvtg cals: " << endl;
    for(int i = 0; i < fuel.batch.size(); i++){

        Sig_aF = Siga_finder(fuel.batch[i]);
        //cout << "Siga: " << Sig_aF << endl;

        Sig_tF = Sig_aF+Sig_sF;
        Sig_trF = Sig_tF - 2/3/A_F*Sig_sF;
        D_F = 1 / (3 * Sig_trF);
        L_F = sqrt(D_F/Sig_aF);
        x = a/L_F;

        /*****book example***
        //should get f = 0.8272 with the values below
        //Lamarsh pg.316
        a=1.02;
        b=14.3;
        x=0.658;
        y=0.0173;
        z=0.242;
        V_M=195.6;
        V_F=1;
        Sig_aM=0.0002728;
        Sig_aF=0.3668;
        ******************/

        F = x * boost::math::cyl_bessel_i(0,x) / (2 * boost::math::cyl_bessel_i(1, x));

        E = (z*z - y*y) / (2 * y) * ( (boost::math::cyl_bessel_i(0, y) * boost::math::cyl_bessel_k(1, z)+ boost::math::cyl_bessel_k(0, y) * boost::math::cyl_bessel_i(1, z)) / (boost::math::cyl_bessel_i(1, z) * boost::math::cyl_bessel_k(1, y) - boost::math::cyl_bessel_k(1, z) * boost::math::cyl_bessel_i(1, y)));

        f = pow((((Sig_aM * V_M)/(Sig_aF * V_F)) * F + E), (-1.));
        //cout << f << "  Disadvtg: " << (Sig_aF*V_F - f*Sig_aF*V_F)/(f*Sig_aM*V_M)<<endl;

        fuel.batch[i].DA = (Sig_aF*V_F - f*Sig_aF*V_F)/(f*Sig_aM*V_M);
    }

    return fuel;
}


double SS_burnupcalc(isoInformation fuel, int N, double delta, double PNL, double base_flux){
    double burnup = 0;
    double dt = delta*24*60*60; //days to [s]
    fuelBundle core;
    core.pnl = PNL;
    core.base_flux = base_flux;
    double BU, BU_est, F_est;
    double y0, y1;
    double kcore_prev;

    //find when k drops under 1 for single batch
    int ii = 0;
    for(ii = 0; fuel.neutron_prod[ii]*PNL/fuel.neutron_dest[ii] >= 1; ii++){
        if(ii == fuel.neutron_prod.size()-1){
            cout << endl << "SS_burnupcalc error! Fuel criticality doesn't drop below 1." << endl;
            ii -= 1;
            return fuel.BU[ii+1];
        }
    }

    //find the just critical fluence for single batch
    F_est = intpol(fuel.fluence[ii-1], fuel.fluence[ii], fuel.neutron_prod[ii-1]*PNL/fuel.neutron_dest[ii-1], fuel.neutron_prod[ii]*PNL/fuel.neutron_dest[ii], 1);

    //find the corresponding BU
    BU_est = intpol(fuel.BU[ii-1], fuel.BU[ii], fuel.fluence[ii-1], fuel.fluence[ii], F_est);

    //estimate the N batch BU
    BU_est = 2. * N * BU_est / (N + 1.);

    //assign the linearly divided burnup and fluence to each batch
    for(int i = 0; i < N; i++){
        batch_info temp_batch;
        temp_batch.collapsed_iso = fuel;
        temp_batch.BUg = BU_est * i / N; //first batch gets zero BU


        for(ii = 0; temp_batch.collapsed_iso.BU[ii] < temp_batch.BUg; ii++){}
        temp_batch.Fg = intpol(temp_batch.collapsed_iso.fluence[ii-1], temp_batch.collapsed_iso.fluence[ii], temp_batch.collapsed_iso.BU[ii-1], temp_batch.collapsed_iso.BU[ii], temp_batch.BUg);

        if(temp_batch.Fg < 0){temp_batch.Fg = 0;}

        core.batch.push_back(temp_batch);
    }

    //this function ignores structural material effects (for now)
    core.struct_prod = 0;
    core.struct_dest = 0;
    for(int i = 0; i < core.batch.size(); i++){
        core.batch[i].DA = 0;
    }


    double kcore = 3.141592;

    while(kcore > 1){
        kcore_prev = kcore;

        //calculate necessary parameters
        core = phicalc_simple(core);

        for(int i = 0; i < N; i++){
            core.batch[i].Fg += core.batch[i].rflux * core.base_flux * dt;
            cout << "  Fg: " << core.batch[i].Fg << endl;
        }
        kcore = kcalc(core);
        cout << "kcalc:" << kcore  << endl;
    }

    //update core fluences
    for(int i = 0; i < N; i++){
        //y0 is the fluence value before the last interation
        y0 = core.batch[i].Fg - core.batch[i].rflux * core.base_flux * dt;
        y1 = core.batch[i].Fg;
        //cout << y0 << " and " << y1 << endl; //<-------
        //cout << "  " << kcore_prev << " " << kcore << endl; //<-------
        core.batch[i].batch_fluence = intpol(y0, y1, kcore_prev, kcore, 1);
        //cout << "  fluence end of burnupcalc: " << core.batch[i].batch_fluence << endl;

    }

    for(ii = 0; core.batch[N-1].collapsed_iso.fluence[ii] < core.batch[N-1].batch_fluence; ii++){}
    if(core.batch[N-1].collapsed_iso.fluence.back() < core.batch[N-1].batch_fluence){
        cout << endl << "Maximum fluence error!(SS_burnupcalc) Batch fluence exceeded max library fluence. (burnupcalc2)" << endl;
        cout << "  Values on max fluence will be used. Do not trust results." << endl;
        ii = core.batch[N-1].collapsed_iso.fluence.size() - 1;
    }
    burnup = intpol(core.batch[N-1].collapsed_iso.BU[ii-1], core.batch[N-1].collapsed_iso.BU[ii], core.batch[N-1].collapsed_iso.fluence[ii-1], core.batch[N-1].collapsed_iso.fluence[ii], core.batch[N-1].batch_fluence);

    return burnup;
}


/*std::pair<double, std::pair<double, map<int,double> > > blending_calc(fuelBundle fuel, double BU_end, int mode, int da_mode, double time_step) {
    std::pair<double, std::pair<double, std::map<int, double> > > rtn;
    double X;
    double BU_guess, enrich_guess;
    double enrich_previous, BU_previous;
    double enrich_lower, enrich_upper;
    double BU_lower, BU_upper;
    std::vector <int> blend_vector;

    // This is a super quick hack to benchmark enrichcalc
    fuelBundle fuel1;
    fuelBundle fuel2;
    for (int i = 0; i < fuel.stream_fraction.size(); i++){
        if (fuel.stream_fraction[i] == -1){
            blend_vector.push_back(i);
        }
    }
    for (int i = 1; i < 100; i++){
        X = i / 100.;
        fuel.stream_fraction[blend_vector[0]] = X;
        fuel.stream_fraction[blend_vector[1]] = 1 - X;
        fuel = enrich_collapse(fuel);
        BU_lower = burnupcalc(fuel, mode, da_mode, time_step);
        //cout << X << "  BU_LOWER   " << BU_lower << endl;
        if (i > 99){
            cout << "IT'S ALL BROKEN" << endl;
            rtn.first = 0;
            return rtn;
        }
        if (BU_lower > 0){
            enrich_lower = X;
            break;
        }
    }

    for (int i = enrich_lower*100; i < enrich_lower*100+10; i++){
        X = i / 100.;
        fuel.stream_fraction[blend_vector[0]] = X;
        fuel.stream_fraction[blend_vector[1]] = 1 - X;
        fuel1 = enrich_collapse(fuel);
        BU_upper = burnupcalc(fuel1, mode, da_mode, time_step);
        //cout << X << "  BU_UPPER   " << BU_upper << endl;
        if (BU_upper > 2*BU_lower){
            enrich_upper = X;
            break;
        }
    }
    X = enrich_lower + (BU_end - BU_lower)*(enrich_upper - enrich_lower)/(BU_upper - BU_lower);
    if (BU_upper == 0){
        cout << "Burn up code failed" << endl;
        rtn.first = 0;
        return rtn;
    }
    fuel.stream_fraction[blend_vector[0]] = X;
    fuel.stream_fraction[blend_vector[1]] = 1 - X;
    fuel2 = enrich_collapse(fuel);
    BU_guess = burnupcalc(fuel2, mode, da_mode, time_step);
    if (BU_guess == 0){
        cout << "Burn up code failed" << endl;
        rtn.first = 0;
        return rtn;
    }
    // enrichment iteration
    if (abs(enrich_guess - enrich_lower) > abs(enrich_guess - enrich_upper)){
        enrich_previous = enrich_lower;
        BU_previous = BU_lower;
    } else {
        enrich_previous = enrich_upper;
        BU_previous = BU_upper;
    }
    enrich_guess = X;
    int i = 0;
    while ((abs(BU_end - BU_guess)/BU_end) > 0.001){
        //cout << enrich_previous << "   " << BU_previous << endl;
        //cout << enrich_guess << "   " << BU_guess << endl << endl;
        fuelBundle fuel3;
        X = enrich_previous + (BU_end - BU_previous)*(enrich_guess - enrich_previous)/(BU_guess - BU_previous);
        BU_previous = BU_guess;
        enrich_previous = enrich_guess;
        fuel.stream_fraction[blend_vector[0]] = X;
        fuel.stream_fraction[blend_vector[1]] = 1 - X;
        fuel3 = enrich_collapse(fuel);
        BU_guess = burnupcalc(fuel3, mode, da_mode, time_step);
        enrich_guess = X;
        /// TODO FIX THIS QUICK HACK
        if (i > 20) {
            X = (enrich_guess + enrich_previous)/2;
            break;
        }
        i++;
    }
    rtn.first = X;
    fuelBundle fuel3;
    fuel.stream_fraction[blend_vector[0]] = X;
    fuel.stream_fraction[blend_vector[1]] = 1 - X;
    fuel3 = enrich_collapse(fuel);
    rtn.second = burnupcalc(fuel3, mode, da_mode, time_step);
    return rtn;
}*/


/*
fuelBundle InputReader(){
    std::string name, fraction;
    int region, N, t_res;
    char type, word[8];
    int nucid;
    double pnl, intpol_val, target_BUd, mass;
    fuelBundle fuel;

    string line;
    ifstream fin("inputFile2.txt");

    int i=0;
	while(getline(fin, line))
	{
        if(line.find("REACTOR") == 0){
            istringstream iss(line);
            iss >> name >> name;
            fuel.name = name;
        }
        if (line.find("BURNUP") == 0){
            fuel.operation_type = "BURNUP";
        }
        if (line.find("BLENDING") == 0){
            fuel.operation_type = "BLENDING";
            istringstream iss(line);
            iss >> name >> target_BUd;
            fuel.target_BUd = target_BUd;
        }
        if(line.find("STREAMS") == 0){
            istringstream iss(line);
            iss >> name;
            while (iss >> fraction){
                if (fraction == "X"){
                    fuel.stream_fraction.push_back(-1.);
                } else {
                    fuel.stream_fraction.push_back(atof(fraction.c_str()));
                }
            }
        }
        if(line.find("REGIONS") == 0){
            while(getline(fin, line)){
                if(line.find("END") == 0) break;
                isoInformation temp;
                istringstream iss(line);
                iss >> region >> type >> nucid;
                temp.name = nucid;
                temp.region = region;
                temp.type = type;
                temp.fraction.push_back(0.0);
                while (iss >> mass){
                    temp.fraction.push_back(mass);
                }
                temp.fuel = false;
                fuel.all_iso.push_back(temp);
            }
        }
        if(line.find("BATCH") == 0){
            istringstream iss(line);
            iss >> word >> N;
        }
        if(line.find("FUELRES") == 0){
            istringstream iss(line);
            iss >> word >> t_res;
        }
        if(line.find("LEAK") == 0){
            istringstream iss(line);
            iss >> word >> pnl;
        }
        if(line.find("INTERPOLATE") == 0){
            fuel.libcheck = true;
            while(getline(fin, line)){
                if(line.find("INTERPOLEND") == 0) break;
                if(line.find("INTERPOL") == 0){
                    istringstream iss(line);
                    iss >> name;
                    while (iss >> name >> intpol_val){
                        interpol_pair lib_pair;
                        lib_pair.metric = name;
                        lib_pair.value = intpol_val;
                        fuel.interpol_pairs.push_back(lib_pair);
                    }
                }
                if(line.find("INTLIBS") == 0){
                     istringstream iss(line);
                     iss >> name;
                     while (iss >> name){
                        fuel.interpol_libs.push_back(name);
                     }
                }
            }

        }
	}
    fuel.batch = N;
    fuel.pnl = pnl;
    fuel.tres = t_res;
    fin.close();

    return fuel;
}
*/

fuelBundle lib_interpol(fuelBundle input_fuel){
    vector<fuelBundle> fuel_pairs;
    for (int i = 0; i < input_fuel.interpol_libs.size(); i++){
        fuelBundle lib_bundle;
        for(int j = 0; j < input_fuel.all_iso.size(); j++){
            isoInformation iso;
            iso.name = input_fuel.all_iso[j].name;
            iso.fraction = input_fuel.all_iso[j].fraction;
            iso.type = input_fuel.all_iso[j].type;
            iso.blending = input_fuel.all_iso[j].blending;
            iso.region = input_fuel.all_iso[j].region;
            lib_bundle.all_iso.push_back(iso);
        }
        lib_bundle.name = input_fuel.interpol_libs[i];
        lib_bundle.batch = input_fuel.batch;
        lib_bundle.tres = input_fuel.tres;
        lib_bundle.pnl = input_fuel.pnl;
        lib_bundle.target_BU = input_fuel.target_BU;
        fuel_pairs.push_back(lib_bundle);
    }
    for (int i = 0; i < fuel_pairs.size(); i++){
        DataReader2(fuel_pairs[i].name, fuel_pairs[i].all_iso);
    }
    // Reading Metrics //
    vector<vector<double> > metrics;
    std::string metric_name;
    double metric_value;
    for (int i = 0; i < input_fuel.interpol_pairs.size(); i++){
        vector<double> metric_vect;
        metrics.push_back(metric_vect);
        for (int j = 0; j < fuel_pairs.size(); j++){
            std::string line;
            ifstream inf(fuel_pairs[j].name + "/params.txt");
            while(getline(inf, line)){
                if (line.find(input_fuel.interpol_pairs[i].metric) == 0){
                    istringstream iss(line);
                    iss >> metric_name >> metric_value;
                    metrics[i].push_back(metric_value);
                }
            }
            inf.close();
        }
    }

    double alpha = 1;
    for (int i = 0; i < metrics.size(); i++){
        double max = 0;
        double min = metrics[i][0];
        for(int j = 0; j < metrics[i].size(); j++){
            if (metrics[i][j] > max){
                max = metrics[i][j];
            }
            if (metrics[i][j] < min){
                min = metrics[i][j];
            }
        }
        if (input_fuel.interpol_pairs[i].value > max || input_fuel.interpol_pairs[i].value < min){
            alpha = 40;
        }
        for(int j = 0; j < metrics[i].size(); j++){
            metrics[i][j] -= min;
            metrics[i][j] /= max;
        }
        input_fuel.interpol_pairs[i].scaled_value = input_fuel.interpol_pairs[i].value - min;
        input_fuel.interpol_pairs[i].scaled_value /= max;
    }
    // distances
    vector<double> metric_distances;
    for(int i = 0; i < input_fuel.interpol_libs.size(); i++){
        double distance_measure = 1;
        for (int j = 0; j < metrics.size(); j++){
            distance_measure *= pow(input_fuel.interpol_pairs[j].scaled_value - metrics[j][i], 2);
        }
        if(distance_measure == 0){
            return fuel_pairs[i];
        }
        metric_distances.push_back(pow(distance_measure, alpha/2));
        std::cout << distance_measure << std::endl;
    }
    double met_dist_sum;
    for (int i = 0; i < metric_distances.size(); i++){
        met_dist_sum += metric_distances[i];
    }
    std::cout << "Sum : " << met_dist_sum << std::endl;
    // Fuel Bundle instead of iso //
    fuelBundle new_fuel;
    new_fuel.tres = input_fuel.tres;
    new_fuel.pnl = input_fuel.pnl;
    new_fuel.batch = input_fuel.batch;
    new_fuel.operation_type = input_fuel.operation_type;
    new_fuel.target_BU = input_fuel.target_BU;
    for(int i = 0; i < fuel_pairs.size(); i++){
        if(i==0){
            for(int j = 0; j < fuel_pairs[i].all_iso.size(); j++){
                isoInformation new_iso;
                new_iso.name = fuel_pairs[i].all_iso[j].name;
                new_iso.region = fuel_pairs[i].all_iso[j].region;
                new_iso.blending = fuel_pairs[i].all_iso[j].blending;
                new_iso.fraction = fuel_pairs[i].all_iso[j].fraction;
                new_iso.type = fuel_pairs[i].all_iso[j].type;
                for (int k = 0; k < fuel_pairs[i].all_iso[j].fluence.size(); k++){
                    new_iso.fluence.push_back(fuel_pairs[i].all_iso[j].fluence[k]);
                }
                for (int k = 0; k < fuel_pairs[i].all_iso[j].neutron_prod.size(); k++){
                    new_iso.neutron_prod.push_back(fuel_pairs[i].all_iso[j].neutron_prod[k]*metric_distances[i]/met_dist_sum);
                }
                for (int k = 0; k < fuel_pairs[i].all_iso[j].neutron_dest.size(); k++){
                    new_iso.neutron_dest.push_back(fuel_pairs[i].all_iso[j].neutron_dest[k]*metric_distances[i]/met_dist_sum);
                }
                for (int k = 0; k < fuel_pairs[i].all_iso[j].BUd.size(); k++){
                    new_iso.BUd.push_back(fuel_pairs[i].all_iso[j].BUd[k]*metric_distances[i]/met_dist_sum);
                }
                for (int k = 0; k < fuel_pairs[i].all_iso[j].iso_vector.size(); k++){
                    daughter new_daughter;
                    new_daughter.name = fuel_pairs[i].all_iso[j].iso_vector[k].name;
                    for(int ii = 0; ii < fuel_pairs[i].all_iso[j].iso_vector[k].mass.size(); ii++){
                        new_daughter.mass.push_back(fuel_pairs[i].all_iso[j].iso_vector[k].mass[ii]*metric_distances[i]/met_dist_sum);
                    }
                    new_iso.iso_vector.push_back(new_daughter);
                }
                new_fuel.all_iso.push_back(new_iso);
            }
        } else {
            for(int j = 0; j < fuel_pairs[i].all_iso.size(); j++){
                bool iso_check1 = false;
                for(int m = 0; m < new_fuel.all_iso.size(); m++){
                    if(fuel_pairs[i].all_iso[j].name == new_fuel.all_iso[m].name){
                        iso_check1 = true;
                        for (int k = 0; k < fuel_pairs[i].all_iso[j].neutron_prod.size(); k++){
                            new_fuel.all_iso[m].neutron_prod[k] += fuel_pairs[i].all_iso[j].neutron_prod[k]*metric_distances[i]/met_dist_sum;
                        }
                        for (int k = 0; k < fuel_pairs[i].all_iso[j].neutron_dest.size(); k++){
                            new_fuel.all_iso[m].neutron_dest[k] += fuel_pairs[i].all_iso[j].neutron_dest[k]*metric_distances[i]/met_dist_sum;
                        }
                        for (int k = 0; k < fuel_pairs[i].all_iso[j].BUd.size(); k++){
                            new_fuel.all_iso[m].BUd[k] += fuel_pairs[i].all_iso[j].BUd[k]*metric_distances[i]/met_dist_sum;
                        }
                        for (int k = 0; k < fuel_pairs[i].all_iso[j].iso_vector.size(); k++){
                            bool iso_check2 = false;
                            for (int mm = 0; mm < new_fuel.all_iso[m].iso_vector.size(); mm++){
                                if(fuel_pairs[i].all_iso[j].iso_vector[k].name == new_fuel.all_iso[m].iso_vector[mm].name){
                                    iso_check2 = true;
                                    for(int n = 0; n < new_fuel.all_iso[m].iso_vector[mm].mass.size(); n++){
                                        new_fuel.all_iso[m].iso_vector[mm].mass[n] += fuel_pairs[i].all_iso[j].iso_vector[k].mass[n]*metric_distances[i]/met_dist_sum;
                                    }
                                }
                            }
                            if (iso_check2 == false){
                                daughter new_daught;
                                new_daught.name = fuel_pairs[i].all_iso[j].iso_vector[k].name;
                                for(int n = 0; n < fuel_pairs[i].all_iso[j].iso_vector[k].mass.size(); n++){
                                    new_daught.mass.push_back(fuel_pairs[i].all_iso[j].iso_vector[k].mass[n]*metric_distances[i]/met_dist_sum);
                                }
                                new_fuel.all_iso[m].iso_vector.push_back(new_daught);
                            }
                        }
                    }
                } if (iso_check1 == false){
                    //std::cout << fuel_pairs[i].all_iso[j].name << endl;
                    isoInformation new_iso;
                    new_iso.name = fuel_pairs[i].all_iso[j].name;
                    new_iso.region = fuel_pairs[i].all_iso[j].region;
                    new_iso.blending = fuel_pairs[i].all_iso[j].blending;
                    new_iso.fraction = fuel_pairs[i].all_iso[j].fraction;
                    new_iso.type = fuel_pairs[i].all_iso[j].type;
                    for (int k = 0; k < fuel_pairs[i].all_iso[j].fluence.size(); k++){
                        new_iso.fluence.push_back(fuel_pairs[i].all_iso[j].fluence[k]);
                    }
                    for (int k = 0; k < fuel_pairs[i].all_iso[j].neutron_prod.size(); k++){
                        new_iso.neutron_prod.push_back(fuel_pairs[i].all_iso[j].neutron_prod[k]*metric_distances[i]/met_dist_sum);
                    }
                    for (int k = 0; k < fuel_pairs[i].all_iso[j].neutron_dest.size(); k++){
                        new_iso.neutron_dest.push_back(fuel_pairs[i].all_iso[j].neutron_dest[k]*metric_distances[i]/met_dist_sum);
                    }
                    for (int k = 0; k < fuel_pairs[i].all_iso[j].BUd.size(); k++){
                        new_iso.BUd.push_back(fuel_pairs[i].all_iso[j].BUd[k]*metric_distances[i]/met_dist_sum);
                    }
                    for (int k = 0; k < fuel_pairs[i].all_iso[j].iso_vector.size(); k++){
                        daughter new_daughter;
                        new_daughter.name = fuel_pairs[i].all_iso[j].iso_vector[k].name;
                        for(int ii = 0; ii < fuel_pairs[i].all_iso[j].iso_vector[k].mass.size(); ii++){
                            new_daughter.mass.push_back(fuel_pairs[i].all_iso[j].iso_vector[k].mass[ii]*metric_distances[i]/met_dist_sum);
                        }
                    }
                    new_fuel.all_iso.push_back(new_iso);
                }
            }
        }
    }
    return new_fuel;
}

void mass_check(fuelBundle fuel){
    for(int i = 0; i < fuel.all_iso.size(); i++){
        for(int k = 0; k < fuel.all_iso[i].fluence.size(); k++){
            double mass = 0;
            for(int j = 0; j < fuel.all_iso[i].iso_vector.size(); j++){
                mass += fuel.all_iso[i].iso_vector[j].mass[k];
            }
            std::cout << mass << std::endl;
        }
    }
}
/*
fuelBundle burnup_collapse(fuelBundle fuel){
    for(int i = 0; i < fuel.iso.size(); i++){
        fuel.iso[i].fraction[0] = fuel.iso[i].fraction[1];
    }
    return fuel;
}

void iso_output(pair<double, map <int, double> > iso_vector){
    typedef std::map<int, double>::iterator test_map;
    for (test_map iterator = iso_vector.second.begin(); iterator != iso_vector.second.end(); iterator++){
        cout << iterator->first << "       " << iterator->second << endl;
    }
}
*/
/*
int main(){

    fuelBundle fuel;

    vector<int> iso_index;
    vector<fuelBundle> batches;


    fuel = InputReader();

    fuel = fluxcalc_reader(fuel, "fluxCalcinput");
    DataReader2(fuel.name, fuel.iso);

    vector<nonActinide> nona; //"NONA"ctinide ;)
    nona = NonActinideReader(fuel.name + "/TAPE9.INP");
    fuel = NBuilder(fuel, nona);


    //this should probably be its own function
    //builds bundle-vector and assigns fluence
    batches.push_back(fuel);
    batches[0].batch_fluence = 0;
    for(int i = 0; i < fuel.batch-1; i++){
        batches.push_back(batches[0]);
        batches[i+1].batch_fluence = batches[i].batch_fluence + 600.2; //123.2 is a random choice
    }

    burnupcalc(batches, fuel.pnl, 0.001, 1);

    if (fuel.libcheck == true){
        fuel = lib_interpol(fuel);
    }

    if (fuel.operation_type == "BURNUP"){
        fuel = burnup_collapse(fuel);
        pair< double, map < int, double> > test = burnupcalc(fuel, fuel.batch, fuel.pnl, 0.001);
        cout <<"burnup: "<< test.first << endl;
        iso_output(test);
    } else {
        pair< double, pair<double, map<int, double> > > enrichment = enrichcalc(fuel.target_BUd, fuel.batch, 0.001, fuel);
        if (enrichment.first == 0){
            cout << "Reactor could not reach desired burnup." << endl;
        } else {
            cout << "Enrichment:    " << enrichment.first << endl;
        }
        pair< double, map < int, double> > test = enrichment.second;
        iso_output(test);
    }

    return 0;
}

*/
