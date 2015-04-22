#include "burnupcalc.h"
#include <boost/math/special_functions/bessel.hpp>


using namespace std;

typedef unsigned long long timestamp_t;

static timestamp_t
get_timestamp ()
    {
      struct timeval now;
      gettimeofday (&now, NULL);
      return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
    }

double intpol(double y0, double y1, double x0, double x1, double x) {
    // linear interpolation function
    double y = y0 + (y1 - y0)*(x - x0)/(x1 - x0);
    return y;
}


fuelBundle CoreCollapse(fuelBundle &fuel){
    for(int i = 0; i < fuel.batch.size(); i++){
        fuel.batch[i] = BatchCollapse(fuel.batch[i]);
    }

    return fuel;
}

batch_info BatchCollapse(batch_info &batch){
    //cout << " 1  " << batch.iso.size() << endl;
    batch.collapsed_iso = FuelBuilder(batch.iso);

    //builds total BU from BUd
    //cout << " 2" << endl;
    batch.collapsed_iso.BU.push_back(batch.collapsed_iso.BUd[0]);

    for(int j = 1; j < batch.collapsed_iso.BUd.size(); j++){
    //cout << "    1.5tst" << i+1 << "  " << j << endl;
        batch.collapsed_iso.BU.push_back(batch.collapsed_iso.BU[j-1]+batch.collapsed_iso.BUd[j]);
    }
    //test to see if the prod/dest vectors are the same length
    if(batch.collapsed_iso.neutron_prod.size() != batch.collapsed_iso.neutron_dest.size()){
        cout << "Error. Neutron production/destruction rate vector length mismatch." << endl;
    }
    if(batch.collapsed_iso.neutron_prod[0] == 0 || batch.collapsed_iso.neutron_dest[0] == 0){
        batch.collapsed_iso.neutron_prod.erase(batch.collapsed_iso.neutron_prod.begin());
        batch.collapsed_iso.neutron_dest.erase(batch.collapsed_iso.neutron_dest.begin());
        batch.collapsed_iso.BUd.erase(batch.collapsed_iso.BUd.begin());
        batch.collapsed_iso.BU.erase(batch.collapsed_iso.BU.begin());
        batch.collapsed_iso.fluence.erase(batch.collapsed_iso.fluence.begin());
        for(int j = 0; j < batch.collapsed_iso.iso_vector.size(); j++){
            batch.collapsed_iso.iso_vector[j].mass.erase(batch.collapsed_iso.iso_vector[j].mass.begin());
        }
    }
    //cout << " 3" << endl;
    return batch;
}

fuelBundle BatchCollapse_old(fuelBundle &fuel){
///add micro region flux effects
//struct effects accounted here
    cout << "Begin regionCollapse_old" << endl;
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


fuelBundle fast_region_collapse(fuelBundle &fuel){
///add micro region flux effects
//struct effects accounted here
    //cout << "Begin regionCollapse" << endl;
    for(int i = 0; i < fuel.batch.size(); i++){
        //for(int j = 0; j < fuel.batch[i].iso.size(); j ++){
        fuel.batch[i].collapsed_iso = BurnupBuilder(fuel.batch[i].iso);
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
        }
    }
    return fuel;
};


map<int, double> tomass (int ti, double fluence, isoInformation &isoinfo) {
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


fuelBundle phicalc_simple(fuelBundle &core){
//updates the rflux of each batch in core.batch
//assumes the flux of a batch is proportional to the inverse neutron prod rate
    double maxphi = 0;
    //finds the inverse of neutron production at the batch_fluence
    //stores it in core.batch.rflux


    for(int i = 0; i < core.batch.size(); i++){
        int ii;
        if(core.batch[i].collapsed_iso.fluence.back() < core.batch[i].Fg){
            //cout << endl << "Maximum fluence error! Batch fluence exceeded max library fluence. (method2)" << endl;
            //cout << " Values on max fluence will be used. Do not trust results." << endl;
            for(int j = 0; j < core.batch.size(); j++){
                core.batch[j].rflux = 1;
            }
            return core;
        } else {
            for(ii = 0; core.batch[i].collapsed_iso.fluence[ii] < core.batch[i].Fg; ii++){}

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
        //cout << core.batch[i].rflux << endl;
    }
    //cout << endl;
    return core;
}

fuelBundle phicalc_eqpow(fuelBundle &core){
    //updates the rflux of each batch in core.batch
    //assumes the flux of a batch is proportional to the inverse neutron prod rate
    double maxphi = 0;
    //finds the inverse of neutron production at the batch_fluence
    //stores it in core.batch.rflux
    for(int i = 0; i < core.batch.size(); i++){
        int ii;


        if(core.batch[i].collapsed_iso.fluence.back() < core.batch[i].Fg){
            //cout << endl << "Maximum fluence error! Batch fluence exceeded max library fluence. (method2)" << endl;
            //cout << "  Values on max fluence will be used. Do not trust results." << endl;
            for(int j = 0; j < core.batch.size(); j++){
                core.batch[j].rflux = 1;
            }
            return core;

        } else {
            for(ii = 0; core.batch[i].collapsed_iso.fluence[ii] < core.batch[i].Fg; ii++){}
        }

        if(ii == 0){
            //core.batch[i].rflux = 1/core.batch[i].collapsed_iso.neutron_prod[0];
            core.batch[i].rflux = (core.batch[i].collapsed_iso.fluence[ii+1]-core.batch[i].collapsed_iso.fluence[ii])
                /(core.batch[i].collapsed_iso.BU[ii+1]-core.batch[i].collapsed_iso.BU[ii]);

        } else{

            core.batch[i].rflux = (core.batch[i].collapsed_iso.fluence[ii]-core.batch[i].collapsed_iso.fluence[ii-1])
                /(core.batch[i].collapsed_iso.BU[ii]-core.batch[i].collapsed_iso.BU[ii-1]);
        }

        if(maxphi < core.batch[i].rflux){
            maxphi = core.batch[i].rflux;
        }
    }

    //normalizes all the flux values
    //cout << "EqPow Flux: ";
    for(int i = 0; i < core.batch.size(); i++){
        core.batch[i].rflux /= maxphi;
        //cout << core.batch[i].rflux << "   ";
    }
    //cout << endl;

    return core;
}

double nusigf_finder(batch_info &batch){
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


    if(ii == 0){
        prod = batch.collapsed_iso.neutron_prod[0];
    }

    //Nusig_f[cm-1] = prod rate / Avagadros number * mass * barn
    //Nusig_f = prod/(6.0221413E+23)*A*1E24/1000;
    Nusig_f = prod*0.01097;

    return Nusig_f;
}

double siga_finder(batch_info &batch){
//returns the macroscopic cross section at fluence Fg
    double sig_a;
    int ii;
    double dest = 0;
    int A = 235;

    //find the index to interpolate on
    for(ii = 0; batch.collapsed_iso.fluence[ii] < batch.Fg; ii++);
    if(batch.collapsed_iso.fluence.back() < batch.Fg){
        ii = batch.collapsed_iso.fluence.size() - 1;
    }

    if(ii == 0){
        dest = batch.collapsed_iso.neutron_dest[0];
    }


    dest = intpol(batch.collapsed_iso.neutron_dest[ii-1], batch.collapsed_iso.neutron_dest[ii], batch.collapsed_iso.fluence[ii-1], batch.collapsed_iso.fluence[ii], batch.Fg);

    //sig_a[barn] = dest rate [cm2/kg] / Avagadros number * mass * barn
    //sig_a = dest/(6.0221413E+23)*A*1E24/1000;
    sig_a = dest*0.01097;

    return sig_a;
}


////////////////////////////////////////////////////////////////////////////////////////////
fuelBundle phicalc_cylindrical(fuelBundle &core){
//cout << "cylindrical" << endl<<endl<<endl<<endl<<endl;

//the word 'region' is interchanged with the word 'batch' in comments, water layer is the outermost region

    int region = core.batch.size(); // number of regions (batches) not counting the outer region
    double delta = core.cylindrical_delta;
    double R[region+1]; //radial thickness of each region
    int N[region+1]; //number of mesh points in each region
    int NC[region+1]; //cumulative N
    int NTotal; //total number of mesh points
    double dd2[region+1]; // D/delta^2
    double Sigma_a[region+1]; //mac. abs. cs of each region
    double NuSigma_f[region+1]; //nu sigma f of each region
    double Sigma_tr[region+1]; //mac. transport cs of each region
    double D[region+1]; //diff coef. for each region
    double LSquared[region+1];
    double k = 1;
    double k_prev = 0.9;
    double prod, prod_prev;
    double flux[region+1], maxflux=0;
    double sum = 0;

    //set the radial thickness of each region
    R[0] = sqrt(core.fuel_area/region/3.141592);
    for(int i = 1; i < region; i++){
        R[i] = sqrt(core.fuel_area/region/3.141592*(i+1));
    }
    R[region] = R[region-1] + core.mod_thickness; //this is the moderator region


/*
    Sigma_a[0] = 0.0230;
    Sigma_a[1] = 0.0246;
    Sigma_a[2] = 0.0324;
    NuSigma_f[0] = 0.0184;
    NuSigma_f[1] = 0.0217;
    NuSigma_f[2] = 0.0382;*/
    //assign fuel cross sections
    for(int i = 0; i < region; i++){

        NuSigma_f[i] = nusigf_finder(core.batch[i]);

        Sigma_a[i] = siga_finder(core.batch[i]);

        Sigma_tr[i] = core.fuel_Sig_tr;
        D[i] = 1/(Sigma_tr[i]*3.);
        LSquared[i] = D[i]/Sigma_a[i];

    }

    //assign moderator cross sections
    Sigma_a[region] = 0.0066;// core.mod_Sig_a;
    Sigma_tr[region] = core.mod_Sig_tr;
    D[region] = 1/(Sigma_tr[region]*3.);
    LSquared[region] = D[region]/Sigma_a[region];
    NuSigma_f[region] = core.mod_Sig_f;

    //populate dd2
    for(int i = 0; i < region+1; i++){
        dd2[i] = D[i]/delta/delta;
    }

    /*
    cout << "Cross sections in cylindrical calc: " << endl << "Sig_a      Sig_f      R" << endl;
    for(int r = 0; r < region+1; r++){
        cout << Sigma_a[r] << "  " << NuSigma_f[r] << "  " << R[r] << endl;
    }*/

    //populate N, number of mesh points in each region
    N[0] = R[0]/delta;
    NC[0] = N[0];
    NTotal = N[0];

    for(int i = 1; i < region+1; i++){
        N[i] = ceil((R[i] - R[i-1]) / delta);
        NC[i] = NC[i-1] + N[i];
        NTotal += N[i];
    }
    NC[region] += 1;
    NTotal += 1;/*
    Eigen::MatrixXf A(NTotal, NTotal);
    Eigen::MatrixXf F(NTotal, 1);
    Eigen::MatrixXf phi(NTotal, 1);
    Eigen::MatrixXf phi_prev(NTotal, 1);
    Eigen::MatrixXf S(NTotal, 1);
    Eigen::MatrixXf S_prev(NTotal, 1);

    A.setZero();
    F.setZero();

    int jprev = 0;
    int j;

    int r = 0; //region index
    for(int i = 1; i < NTotal-1; i++){
        A(i, i-1) = (-1.)*dd2[r]*(2*i-1)/(2*i);
        A(i,i) = dd2[r]*2. + Sigma_a[r];
        A(i, i+1) = (-1.)*dd2[r]*(2*i+1)/(2*i);
        if(i == NC[r]){
            r += 1;
        }
        if(core.batch[r].rflux < 1.1 && core.batch[r].rflux > 0){
            phi_prev(i) = core.batch[r].rflux; //uses last runs results if available
        } else{
            phi_prev(i) = 1;
        }

    }
    A(0,0) = 1;
    A(0,1) = -1;
    A(NTotal-1,NTotal-1) = 1;
    phi_prev(NTotal-1) = 0;

    for(r = 0; r < region; r++){
        A(NC[r],NC[r]-1) = D[r];
        A(NC[r],NC[r]) = -D[r]-D[r+1];
        A(NC[r],NC[r]+1) = D[r+1];
    }

    r = 0;
    for(int i = 1; i < NTotal; i++){
        if(i != NC[r]){
            F(i) = NuSigma_f[r];
        }
        if( i == NC[r]){
            r += 1;
        }
        S_prev(i) = F(i)*phi_prev(i);
    }


    for(int iter = 0; iter < 10; iter++){
        phi = A.colPivHouseholderQr().solve(S_prev)/k_prev;


        if(!phi.allFinite()){
            phi = phi_prev;
        }

        S = F.array() * phi.array();

        prod = (0.25)*3.141592*NuSigma_f[0]*phi(0)*delta*delta;
        prod_prev = (0.25)*3.141592*NuSigma_f[0]*phi_prev(0)*delta*delta;
        r = 0;

        for(int i = 0; i < NTotal; i++){
            if(i == NC[r]){
                prod += 3.141592*NuSigma_f[r]*phi(i)*(i-0.25)*delta*delta;
                prod_prev += 3.141592*NuSigma_f[r]*phi_prev(i)*(i-0.25)*delta*delta;
                r += 1;
                prod += 3.141592*NuSigma_f[r]*phi(i)*(i-0.25)*delta*delta;
                prod_prev += 3.141592*NuSigma_f[r]*phi_prev(i)*(i-0.25)*delta*delta;
            } else {
                prod += 2.*3.141592*NuSigma_f[r]*phi(i)*i*delta*delta;
                prod_prev += 2.*3.141592*NuSigma_f[r]*phi_prev(i)*i*delta*delta;
            }
        }

        if(abs((k_prev-k)/k) < 0.01 && iter > 3){
            break;
        }
        //cout << "prod: " << prod << "  prod_prev: " << prod_prev << "  k: " << prod/prod_prev << endl;
        k = prod/prod_prev*k_prev;
        phi_prev = phi;
        k_prev = k;
        S_prev = S;


    }

    //NO NEED TO normalize phi
    //phi = phi.array()/phi.maxCoeff();

    //find area weighted average phi per batch
    r = 0;
    flux[0] = 0;
    for(int i = 0; i < NTotal; i++){
        flux[r] += phi(i)*(2*(i+1)-1);
        sum += (2*(i+1)-1);

        if(i == NC[r] || i == NTotal-1){
            flux[r] /= sum;
            sum = 0;
            r += 1;
            flux[r] = 0;
        }
    }
    for(r = 0; r < region+1; r++){
        if(flux[r] > maxflux){
            maxflux = flux[r];
        }
    }

    //cout << "flux: ";
    //cout << core.fuel_area << " " << delta << " | " << Sigma_a[0] << " " << Sigma_a[1] << " " << Sigma_a[2] << " | "
    //<< NuSigma_f[0] << " " << NuSigma_f[1] << " " << NuSigma_f[2] << " | ";
    //normalize the fluxes
    for(r = 0; r < region+1; r++){
        flux[r] /= maxflux;
        //cout << flux[r] << " ";
    }
    //ĺcout << endl;




    //NO NEED TO normalize phi
    //phi = phi.array()/phi.maxCoeff();

    //find area weighted average phi per batch
    r = 0;
    flux[0] = 0;
    for(int i = 0; i < NTotal; i++){
        flux[r] += phi(i)*(2*(i+1)-1);
        sum += (2*(i+1)-1);

        if(i == NC[r] || i == NTotal-1){
            flux[r] /= sum;
            sum = 0;
            r += 1;
            flux[r] = 0;
        }
    }
    for(r = 0; r < region+1; r++){
        if(flux[r] > maxflux){
            maxflux = flux[r];
        }
    }

    //cout << "flux: ";
    //cout << core.fuel_area << " " << delta << " | " << Sigma_a[0] << " " << Sigma_a[1] << " " << Sigma_a[2] << " | "
    //<< NuSigma_f[0] << " " << NuSigma_f[1] << " " << NuSigma_f[2] << " | ";
    //normalize the fluxes
    for(r = 0; r < region+1; r++){
        flux[r] /= maxflux;
        //cout << flux[r] << " ";
    }
    //cout << endl;


    for(int i = 0; i < core.batch.size(); i++){
        core.batch[i].rflux = flux[i];
    }
*/
    return core;
}
/////////////////////////////////////////////////////////////////////////////////////////////

double kcalc(fuelBundle &core){
    //uses the fluence values in core.batch.Fg to calculate the core criticality
    //boost::timer t;
    int N = core.batch.size();
    double prod_tot = 0;
    double dest_tot = 0;
    double pnl = core.pnl;


    //finds the index j where fluence is just under target, then interpolates
    //  to find the prod and dest values for each batch
    int j;
    for(int i = 0; i < N; i++){

        if(core.batch[i].collapsed_iso.fluence.back() < core.batch[i].Fg){
            //cout << endl << "Maximum fluence error! Batch fluence exceeded max library fluence. (kcalc)" << endl;
            //cout << "  Values on max fluence will be used. Do not trust results." << endl;
            j = core.batch[i].collapsed_iso.fluence.size() - 1;
        } else {
            //find dicrete point j to interpolate on
            for(j = 0; core.batch[i].collapsed_iso.fluence[j] < core.batch[i].Fg; j++){}
        }

        //add the production rate of batch i to total production
        prod_tot += intpol(core.batch[i].collapsed_iso.neutron_prod[j-1], core.batch[i].collapsed_iso.neutron_prod[j], core.batch[i].collapsed_iso.fluence[j-1], core.batch[i].collapsed_iso.fluence[j], core.batch[i].Fg);

        //add structural material production of this batch, scaled up using disadvantage factor
        //prod_tot += core.struct_prod * core.batch[i].DA;

        //add the destruction rate of batch i to total destruction
        dest_tot += intpol(core.batch[i].collapsed_iso.neutron_dest[j-1], core.batch[i].collapsed_iso.neutron_dest[j], core.batch[i].collapsed_iso.fluence[j-1], core.batch[i].collapsed_iso.fluence[j], core.batch[i].Fg);

        //add structural material destruction of this batch, scaled up using disadvantage factor
        //dest_tot += core.struct_dest * core.batch[i].DA;
        //cout << "diff: " << prod_tot - dest_tot << endl;
    }
    if(pnl < 0 || pnl > 1.5){
        cout << endl << "Error in core nonleakage! Assumed new value: 1.000" << endl;
        pnl = 1.000;
        }

    //std::cout << "k calc time " << t.elapsed() << std::endl;
    return prod_tot * pnl / dest_tot;
}


double CR_batch(fuelBundle &core, int i){
    //
    //i is the batch number starting from zero
    //boost::timer t;
    double FP = 0, FP0 = 0, FP1 = 0;
    double fissile = 0, fissile0 = 0, fissile1 = 0;
    double ini_fissile = 0;
    double CR;
    int ii, ZZ;

    //cout << "\nCR_finder starting." << endl << "  Upper: " << core.CR_upper << "  lower: " << core.CR_lower << endl;



    for(ii = 0; core.batch[i].collapsed_iso.fluence[ii] < core.batch[i].Fg; ii++){}
    if(ii == 0){
        ii = 1;
    } else if (ii > core.batch[i].collapsed_iso.fluence.size()-1){
        ii = core.batch[i].collapsed_iso.fluence.size()-1;
    }

    //cout << " fluence: " << core.batch[i].Fg;

    for(int j = 0; j < core.batch[i].collapsed_iso.iso_vector.size(); j++){
        //convert name to mass number
        ZZ = core.batch[i].collapsed_iso.iso_vector[j].name;
        ZZ = ZZ % 10000;
        ZZ /= 10;

        //add up the FP
        if(ZZ < core.CR_upper && ZZ > core.CR_lower){
            //cout << "    Batch" << i+1 << " ZZ: " << ZZ << endl;
            //interpolation will be done at the end
            FP0 += core.batch[i].collapsed_iso.iso_vector[j].mass[ii-1];
            FP1 += core.batch[i].collapsed_iso.iso_vector[j].mass[ii];
            //cout << "ZZ " << ZZ  << " mass " << core.batch[i].collapsed_iso.iso_vector[j].mass[ii]<<endl;
        }

        //add up fissiles
        for(int fis = 0; fis < core.CR_fissile.size(); fis++){
            if(core.batch[i].collapsed_iso.iso_vector[j].name == core.CR_fissile[fis]){
                fissile0 += core.batch[i].collapsed_iso.iso_vector[j].mass[ii-1];
                fissile1 += core.batch[i].collapsed_iso.iso_vector[j].mass[ii];

                ini_fissile += core.batch[i].collapsed_iso.iso_vector[j].mass[0]; //mass at fluence zero
            }
        }
    }


    // recycling variable FP0 here to check greater than zero
    FP0 = intpol(FP0, FP1, core.batch[i].collapsed_iso.fluence[ii-1], core.batch[i].collapsed_iso.fluence[ii], core.batch[i].Fg);
    if(FP0 > 0){
        FP += FP0;
    }
    //cout << fissile0 << "  " << fissile1 << "   " << core.batch[i].collapsed_iso.fluence[ii-1] << " " << core.batch[i].collapsed_iso.fluence[ii] << "   " << core.batch[i].Fg << endl;
    fissile = intpol(fissile0, fissile1, core.batch[i].collapsed_iso.fluence[ii-1], core.batch[i].collapsed_iso.fluence[ii], core.batch[i].Fg);

    //cout << endl << setprecision(4) <<  "  FP: " << FP << "  fiss: " << fissile << "  ini_fiss: " << ini_fissile << "  num: " << FP+fissile-ini_fissile << "      CR: " << (FP+fissile-ini_fissile)/FP << endl;

    if(FP > 0){
        CR = (FP+fissile-ini_fissile)/FP;
    } else {
        //cout << " -- CR finder (batch " << i+1 << ") has fission product mass of zero." << endl;
        CR  = 0;
    }
    //std::cout << "CR Calc " << t.elapsed() << std::endl;
    return CR;
}


double CR_finder(fuelBundle &core){
    //
    //i is the batch number starting from zero
    //boost::timer t;
    double FP = 0, FP0 = 0, FP1 = 0;
    double fissile = 0, fissile0 = 0, fissile1 = 0;
    double ini_fissile = 0;
    double CR;
    int ii, ZZ;

    //cout << "\nCR_finder starting." << endl << "  Upper: " << core.CR_upper << "  lower: " << core.CR_lower << endl;


    for(int i = 0; i < core.batch.size(); i++){
        for(ii = 0; core.batch[i].collapsed_iso.fluence[ii] < core.batch[i].Fg; ii++){}
        if(ii == 0){
            ii = 1;
        } else if (ii > core.batch[i].collapsed_iso.fluence.size()-1){
            ii = core.batch[i].collapsed_iso.fluence.size()-1;
        }

        //cout << " fluence: " << core.batch[i].Fg;

        for(int j = 0; j < core.batch[i].collapsed_iso.iso_vector.size(); j++){
            //convert name to mass number
            ZZ = core.batch[i].collapsed_iso.iso_vector[j].name;
            ZZ = ZZ % 10000;
            ZZ /= 10;

            //add up the FP
            if(ZZ < core.CR_upper && ZZ > core.CR_lower){
                //cout << "    Batch" << i+1 << " ZZ: " << ZZ << endl;
                //interpolation will be done at the end
                FP0 += core.batch[i].collapsed_iso.iso_vector[j].mass[ii-1];
                FP1 += core.batch[i].collapsed_iso.iso_vector[j].mass[ii];
                //cout << "ZZ " << ZZ  << " mass " << core.batch[i].collapsed_iso.iso_vector[j].mass[ii]<<endl;
            }

            //add up fissiles
            for(int fis = 0; fis < core.CR_fissile.size(); fis++){
                if(core.batch[i].collapsed_iso.iso_vector[j].name == core.CR_fissile[fis]){
                    fissile0 += core.batch[i].collapsed_iso.iso_vector[j].mass[ii-1];
                    fissile1 += core.batch[i].collapsed_iso.iso_vector[j].mass[ii];

                    ini_fissile += core.batch[i].collapsed_iso.iso_vector[j].mass[0]; //mass at fluence zero
                }
            }
        }


        // recycling variable FP0 here to check greater than zero
        FP0 = intpol(FP0, FP1, core.batch[i].collapsed_iso.fluence[ii-1], core.batch[i].collapsed_iso.fluence[ii], core.batch[i].Fg);
        if(FP0 > 0){
            FP += FP0;
        }
        //cout << fissile0 << "  " << fissile1 << "   " << core.batch[i].collapsed_iso.fluence[ii-1] << " " << core.batch[i].collapsed_iso.fluence[ii] << "   " << core.batch[i].Fg << endl;
        fissile += intpol(fissile0, fissile1, core.batch[i].collapsed_iso.fluence[ii-1], core.batch[i].collapsed_iso.fluence[ii], core.batch[i].Fg);

        FP0 = 0;
        FP1 = 0;
        fissile0 = 0;
        fissile1 = 0;
    }

    //cout << endl << setprecision(4) <<  "  FP: " << FP << "  fiss: " << fissile << "  ini_fiss: " << ini_fissile << "  num: " << FP+fissile-ini_fissile << "      CR: " << (FP+fissile-ini_fissile)/FP << endl;

    if(FP > 0){
        CR = (FP+fissile-ini_fissile)/FP;
    } else {
        cout << " CR finder has fission product mass of zero." << endl;;
        CR  = 0;
    }
    //std::cout << "CR Calc " << t.elapsed() << std::endl;
    return CR;
}


void burnupcalc(fuelBundle &core, int mode, int DA_mode, double delta) {
timestamp_t t0 = get_timestamp();
    //this function only uses the COLLAPSED_ISO of each BATCH in the structure CORE
    //all factors that contribute to a change in neutron prod/dest rates have to be factored
    //      before calling this function
    //cout << endl << "Burnupcalc" << endl;
    boost::timer t;

    int N = core.batch.size(); //number of batches
    double dt = delta*24*60*60; //days to [s]
    double kcore, kcore_prev;
    double y0, y1;
    double burnup = 0, burnup_1 = 0;

/*
    cout << "burnupcalc starting neutron prods: " << core.batch[0].collapsed_iso.neutron_prod[0] << " " << core.batch[1].collapsed_iso.neutron_prod[0] << " " << core.batch[2].collapsed_iso.neutron_prod[0]  << endl;
    cout << "   " << core.batch[0].collapsed_iso.BU[1] << " " << core.batch[1].collapsed_iso.BU[1] << " " << core.batch[2].collapsed_iso.BU[1]  << endl;
    cout << "   " << core.batch[0].collapsed_iso.fluence[1] << " " << core.batch[1].collapsed_iso.fluence[1] << " " << core.batch[2].collapsed_iso.fluence[1]  << endl;
    cout << "     " << dt << endl;
    cout << " " << core.batch[0].batch_fluence << " " << core.batch[1].batch_fluence << " " << core.batch[2].batch_fluence  << endl;
*/

    //assign the batch_fluence to Fg
    for(int i = 0; i < N; i++){
        core.batch[i].Fg = core.batch[i].batch_fluence;
        //cout << core.batch[i].Fg << "  " << core.batch[i].collapsed_iso.neutron_prod[0] << endl;
    }
    int jk;
    for(jk = 0; core.batch[0].collapsed_iso.fluence[jk] < core.batch[0].Fg; jk++){}
    burnup_1 = intpol(core.batch[0].collapsed_iso.BU[jk-1], core.batch[0].collapsed_iso.BU[jk], core.batch[0].collapsed_iso.fluence[jk-1], core.batch[0].collapsed_iso.fluence[jk], core.batch[0].batch_fluence);

    kcore = 3.141592;

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
        }else if(mode == 0){
            core = phicalc_eqpow(core);
        }else{
            cout << endl << "Error in mode input for batch-level flux calculation." << endl;
            return;
        }

        //disadvantage calculation
        if(DA_mode == 1){
            core = DA_calc(core);
        }

        //CR calculation
        //cout << CR_numerator(core, 5)/CR_denominator(core, 5) << " ";

        //update fluences
        for(int i = 0; i < N; i++){
            //cout << "flux: " << core.batch[i].rflux << endl;
            core.batch[i].Fg += core.batch[i].rflux * core.base_flux * dt;

        }
        //core.CR = CR_finder(core);
        kcore = kcalc(core);

        /*if(core.CR_target != 0 && std::abs(core.CR - core.CR_target)/core.CR < 0.1){
            break;
        }*/
    }
    core.CR = CR_finder(core);
    core.batch[0].discharge_CR = CR_batch(core, 0);
    //cout << " CR: " << core.CR << endl;


    //update core fluences and CR
    for(int i = 0; i < N; i++){
        core.batch[i].batch_fluence = intpol(core.batch[i].Fg - (core.batch[i].rflux * core.base_flux * dt), core.batch[i].Fg, kcore_prev, kcore, 1);
        //cout << "batch " << i+1 << " CR: " << core.batch[i].CR << endl;
        //cout << "  fluence end of burnupcalc: " << core.batch[i].batch_fluence << endl;

    }

    int ii;
    //update current composition of batches
    for(int i = 0; i < N; i++){
        core.batch[i].comp.clear();
        //cout << "   " << core.batch[i].batch_fluence << "  " ;
        for(ii = 0; core.batch[i].collapsed_iso.fluence[ii] < core.batch[i].batch_fluence; ii++){}
        //cout << "ii: " << ii << endl;
        if(core.batch[i].collapsed_iso.fluence.back() < core.batch[i].batch_fluence){
            cout << endl << "Maximum fluence error! Batch fluence exceeded max library fluence. (burnupcalc1)" << endl;
            cout << "  Values on max fluence will be used. Do not trust results." << endl;
            core.batch[i].batch_fluence = core.batch[i].collapsed_iso.fluence.back();
        }

        double slope = (core.batch[i].batch_fluence - core.batch[i].collapsed_iso.fluence[ii-1])
            /(core.batch[i].collapsed_iso.fluence[ii] - core.batch[i].collapsed_iso.fluence[ii-1]);

        for(int j = 0; j < core.batch[i].collapsed_iso.iso_vector.size(); j++){
            core.batch[i].comp[core.batch[i].collapsed_iso.iso_vector[j].name] =
                (core.batch[i].collapsed_iso.iso_vector[j].mass[ii-1] + (core.batch[i].collapsed_iso.iso_vector[j].mass[ii]
                - core.batch[i].collapsed_iso.iso_vector[j].mass[ii-1])*slope)/1000;
        }
    }

    //the oldest batch is index=0
    for(ii = 1; core.batch[0].collapsed_iso.fluence[ii] < core.batch[0].batch_fluence; ii++){}
    //cout << " ii: " << ii << "  BU[ii-1]: " << core.batch[0].collapsed_iso.BU[ii-1] << endl;
    burnup = intpol(core.batch[0].collapsed_iso.BU[ii-1], core.batch[0].collapsed_iso.BU[ii], core.batch[0].collapsed_iso.fluence[ii-1], core.batch[0].collapsed_iso.fluence[ii], core.batch[0].batch_fluence);

    core.batch[0].delta_BU = burnup - burnup_1;
    core.batch[0].discharge_BU = burnup;
    //core.batch[0].CR = CR_numerator(core, 0)/CR_denominator(core, 0);

timestamp_t t1 = get_timestamp();
}

void burnupcalc_CR(fuelBundle &core, int mode, int DA_mode, double delta) {
    //this function only uses the COLLAPSED_ISO of each BATCH in the structure CORE
    //all factors that contribute to a change in neutron prod/dest rates have to be factored
    //      before calling this function
    //cout << endl << "Burnupcalc" << endl;

    int N = core.batch.size(); //number of batches
    double dt = delta*24*60*60; //days to [s]
    double kcore, kcore_prev;
    double y0, y1;
    double burnup = 0, burnup_1 = 0, burnup_prev;
    double target_burnup = core.target_BU;


    //assign the batch_fluence to Fg
    for(int i = 0; i < N; i++){
        core.batch[i].Fg = core.batch[i].batch_fluence;
        //cout << "BATCH " << i << " FLUENCE " << core.batch[i].batch_fluence << endl;
    }
    int jk;
    for(jk = 0; core.batch[0].collapsed_iso.fluence[jk] < core.batch[0].Fg; jk++){}
    burnup_1 = intpol(core.batch[0].collapsed_iso.BU[jk-1], core.batch[0].collapsed_iso.BU[jk], core.batch[0].collapsed_iso.fluence[jk-1], core.batch[0].collapsed_iso.fluence[jk], core.batch[0].Fg);

    kcore = 3.141592;

    //more forward in time until kcore drops under 1
    while(burnup < target_burnup){
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
        }else if(mode == 0){
            core = phicalc_eqpow(core);
        }else{
            cout << endl << "Error in mode input for batch-level flux calculation." << endl;
            return;
        }

        //disadvantage calculation
        if(DA_mode == 1){
            core = DA_calc(core);
        }

        //CR calculation
        //cout << CR_numerator(core, 5)/CR_denominator(core, 5) << " ";

        //update fluences
        for(int i = 0; i < N; i++){
            //cout << "flux: " << core.batch[i].rflux << endl;
            core.batch[i].Fg += core.batch[i].rflux * core.base_flux * dt;

        }
        int ii;
        for(ii = 0; core.batch[0].collapsed_iso.fluence[ii] < core.batch[0].Fg; ii++){}
        burnup_prev = burnup;
        burnup = intpol(core.batch[0].collapsed_iso.BU[ii-1], core.batch[0].collapsed_iso.BU[ii], core.batch[0].collapsed_iso.fluence[ii-1], core.batch[0].collapsed_iso.fluence[ii], core.batch[0].Fg);
        //cout << "BURNUP " << burnup << "  TARGET " << target_burnup << endl;
    }

    //update core fluences and CR
    for(int i = 0; i < N; i++){
        core.batch[i].batch_fluence = intpol(core.batch[i].Fg - (core.batch[i].rflux * core.base_flux * dt), core.batch[i].Fg, burnup_prev, burnup, target_burnup);
    }

    int ii;
    //update current composition of batches
    for(int i = 0; i < N; i++){
        core.batch[i].comp.clear();
        //cout << "   " << core.batch[i].batch_fluence << "  " ;
        for(ii = 0; core.batch[i].collapsed_iso.fluence[ii] < core.batch[i].batch_fluence; ii++){}
        //cout << "ii: " << ii << endl;
        if(core.batch[i].collapsed_iso.fluence.back() < core.batch[i].batch_fluence){
            cout << endl << "Maximum fluence error! Batch fluence exceeded max library fluence. (burnupcalc1)" << endl;
            cout << "  Values on max fluence will be used. Do not trust results." << endl;
            core.batch[i].batch_fluence = core.batch[i].collapsed_iso.fluence.back();
        }

        double slope = (core.batch[i].batch_fluence - core.batch[i].collapsed_iso.fluence[ii-1])
            /(core.batch[i].collapsed_iso.fluence[ii] - core.batch[i].collapsed_iso.fluence[ii-1]);

        for(int j = 0; j < core.batch[i].collapsed_iso.iso_vector.size(); j++){
            core.batch[i].comp[core.batch[i].collapsed_iso.iso_vector[j].name] =
                (core.batch[i].collapsed_iso.iso_vector[j].mass[ii-1] + (core.batch[i].collapsed_iso.iso_vector[j].mass[ii]
                - core.batch[i].collapsed_iso.iso_vector[j].mass[ii-1])*slope)/1000;
        }
    }
    core.CR = CR_finder(core);
    core.batch[0].discharge_CR = CR_batch(core, 0);
    //the oldest batch is index=0
    //cout << "BURNUP = " << burnup << "  BURNUP_1 = " << burnup_1 << endl;
    core.batch[0].delta_BU = burnup - burnup_1;
    core.batch[0].discharge_BU = burnup;

    //std::cout << "K of core at discharge "<< kcalc(core) << std::endl;
}

double burnupcalc_BU(fuelBundle &core, int mode, int DA_mode, double delta) {
    //this function only uses the COLLAPSED_ISO of each BATCH in the structure CORE
    //all factors that contribute to a change in neutron prod/dest rates have to be factored
    //      before calling this function
    cout << endl << "Burnupcalc "    << endl;

    int N = core.batch.size(); //number of batches
    double dt = delta*24*60*60; //days to [s]
    double kcore, kcore_prev;
    double y0, y1;
    double burnup = 0;

    //assign the batch_fluence to Fg
    for(int i = 0; i < N; i++){
        core.batch[i].Fg = core.batch[i].batch_fluence;
    }

    //turn zero fluences to one to avoid zero rates
    for(int i = 0; i < N; i++){
        if(core.batch[i].Fg == 0){core.batch[i].Fg = 1;}
    }

    /*for(int i = 0; i < core.batch[0].collapsed_iso.fluence.size(); i++){
        //cout << core.batch[0].collapsed_iso.fluence[i] << "  " << core.batch[0].collapsed_iso.BU[i] << "  " << core.batch[0].collapsed_iso.neutron_prod[i]/core.batch[0].collapsed_iso.neutron_dest[i] << endl;
    }*/
    //cout << "Batches: " << core.batch.size() << endl <<endl;
    for(int i = 0; i < core.batch.size(); i++){
        //cout << "fluence: " << core.batch[i].batch_fluence << endl;
    }

    kcore = 3.141592;
    kcore_prev = kcalc(core);
    int counter = 0;
//cout << "counter: ";
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
            return 0;
        }

        //disadvantage calculation
        if(DA_mode == 1){
            core = DA_calc(core);
        }

        //update fluences
        for(int i = 0; i < N; i++){
            //cout << "  Added fluence: " << core.batch[i].rflux * core.base_flux * dt << endl;
            core.batch[i].Fg += core.batch[i].rflux * core.base_flux * dt;
            //cout << "i " << i << "  " <<  core.batch[i].Fg << endl;
        }
        //cout << " " <<counter++;
        kcore = kcalc(core);
        //cout << kcore << endl;
    }

    //core.batch[0].Fg = intpol(y0, y1, kcore_prev, kcore, 1);
    //cout << core.batch[0].Fg <<" asdfa" <<endl;

    //the oldest batch is always index=0
    int ii;
    for(ii = 0; core.batch[0].collapsed_iso.fluence[ii] < core.batch[0].Fg; ii++){}
    if(core.batch[0].collapsed_iso.fluence.back() < core.batch[0].Fg){
        cout << endl << "Maximum fluence error! Batch fluence exceeded max library fluence. (burnupcalc3)" << endl;
        cout << "  Values on max fluence will be used. Do not trust results." << endl;
        ii = core.batch[0].collapsed_iso.fluence.size() - 1;
       //return core.batch[0].collapsed_iso.BU[ii];
    }
    burnup = intpol(core.batch[0].collapsed_iso.BU[ii-1], core.batch[0].collapsed_iso.BU[ii], core.batch[0].collapsed_iso.fluence[ii-1], core.batch[0].collapsed_iso.fluence[ii], core.batch[0].Fg);

    //cout<< "BURNUPUPUPUPU: " << burnup<< endl;
    core.batch[0].discharge_BU = burnup;

    return burnup;
}

fuelBundle DA_calc(fuelBundle &fuel){
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

        Sig_aF = siga_finder(fuel.batch[i]);
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


double SS_burnupcalc(fuelBundle &core, int mode, int DA_mode, double delta, int N, double ss_fluence, double target_burnup){
    //used to find the steady state burnup of the given fuel
    //N:number of batches; delta: burnup time advancement in days; PNL: nonleakage; base_flux: flux of library
    //THE FINAL BATCH IS THE OLDEST ONE
    //cout << "SS burnupcalc begin" << endl;
    boost::timer t;

    isoInformation fuel = core.batch[0].collapsed_iso;
    double burnup_prev;
    double burnup = 0;
    double dt = delta*24*60*60; //days to [s]
    double BU, BU_est, F_est;
    double kcore_prev;
    bool notsteady = true;
    int ii = 0;
    double BU_prev = 0;
    int counter = 0;
    double kcore = 1.10;



    batch_info temp_batch;
    temp_batch.collapsed_iso = fuel;
    temp_batch.batch_fluence = 0;
    for(int i = 1; i < N; i++){
        temp_batch.Fg = ss_fluence/i; // use old cycles youngest batch discharge to guess
        core.batch.push_back(temp_batch);
    }

    while(notsteady){
        int iter = 0;
        burnup = 0;
        kcore_prev = kcore;
        //cout << "  New core" << endl << endl << endl;

        while(burnup < target_burnup){


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
            }else if(mode == 0){
                core = phicalc_eqpow(core);
            } else {
                cout << endl << "Error in mode input for batch-level flux calculation(SS_burnupcalc_CR)." << endl;
                return 0;
            }

            //disadvantage calculation
            if(DA_mode == 1){
                core = DA_calc(core);
            }

            for(int i = 0; i < N; i++){
                //!NOTE when updating this, MUST also update interpolation calcs below
                core.batch[i].Fg += core.batch[i].rflux * core.base_flux * dt;
                //cout << "  Fg: " << core.batch[i].Fg << "  rflux: " << core.batch[i].rflux << endl;
            }
            for(ii = 0; core.batch[0].collapsed_iso.fluence[ii] < core.batch[0].Fg; ii++){}
            //cout << core.batch[0].collapsed_iso.BU[ii-1] << "  " << core.batch[0].collapsed_iso.BU[ii] << "  " <<  core.batch[0].collapsed_iso.fluence[ii-1];
            //cout << " " << core.batch[0].collapsed_iso.fluence[ii] << "  " <<  core.batch[0].Fg << endl;
            burnup = intpol(core.batch[0].collapsed_iso.BU[ii-1], core.batch[0].collapsed_iso.BU[ii], core.batch[0].collapsed_iso.fluence[ii-1], core.batch[0].collapsed_iso.fluence[ii], core.batch[0].Fg);
            kcore = kcalc(core);
            //std::cout << "Burnup " << burnup << std::endl;
            //std::cout << "kcore " << kcore << std::endl;
            iter++;
            if(iter > 20){
                //cout << "  iter+20" << endl;
                //mode = 1;
                //dt *= 2;
                if(iter > 100){
                    cout << "SS_burnupcalc_cr exceeds 100 iterations with error " << endl;
                    cout << "  k: " <<  kcore << "  " << endl;
                    return kcore;
                }
            }
        }

        //update core fluences
        for(int i = 0; i < N; i++){
            //cout << "    " << core.batch[i].rflux;
            core.batch[i].batch_fluence = intpol(core.batch[i].Fg - (core.batch[i].rflux * core.base_flux * dt), core.batch[i].Fg, burnup_prev, burnup, target_burnup);
            if(core.batch[i].batch_fluence > core.batch[i].collapsed_iso.fluence.back()){
            }
        }

        if(abs(kcore - kcore_prev)/kcore < 0.01 && counter > N+1){
            notsteady = false;
        }
        counter++;
        if(counter > 50){
            cout << "Interations exceeded in SS_Burnupcalc core may not reach equalibrium" << endl;
            return kcore;
        }

        //move each batch one left, now the last and second from last are the same
        for(int i = 0; i < N-1; i++){
            core.batch[i].batch_fluence = core.batch[i+1].batch_fluence;
            core.batch[i].Fg = core.batch[i].batch_fluence;
        }

        //fix last one
        core.batch[N-1].batch_fluence = 0;
        core.batch[N-1].Fg = 0;



    }
    //cout << "SSrflux: " << core.batch[1].rflux << endl;
    //std::cout << "SSBurnupCalc " << core.name << " BU: " << burnup << "  Time: " << t.elapsed() << std::endl;
    if(core.batch[0].batch_fluence > 5E23){
        std::cout << "OH YOU DONE MESSED UP, FLUENCE OVER 5E23 at " << core.batch[0].batch_fluence << std::endl;
        //CR = CR+.1;
    }
    return kcore;
}


double SS_burnupcalc_depricated(fuelBundle &core, int mode, int DA_mode, double delta, int N, double ss_fluence){
    //used to find the steady state burnup of the given fuel
    //N:number of batches; delta: burnup time advancement in days; PNL: nonleakage; base_flux: flux of library
    //THE FINAL BATCH IS THE OLDEST ONE
    //cout << "SS burnupcalc begin" << endl;
    boost::timer t;

    isoInformation fuel = core.batch[0].collapsed_iso;
    double burnup = 0;
    double dt = delta*24*60*60; //days to [s]
    double BU, BU_est, F_est;
    double kcore_prev;
    bool notsteady = true;
    int ii = 0;
    double BU_prev = 0;
    int counter = 0;



    batch_info temp_batch;
    temp_batch.collapsed_iso = fuel;
    temp_batch.batch_fluence = 0;
    for(int i = 1; i < N; i++){
        temp_batch.Fg = ss_fluence*(i+1)*0.7; // use old cycles youngest batch discharge to guess
        core.batch.push_back(temp_batch);
    }

    while(notsteady){
        double kcore = 1.10000;
        int iter = 0;
        BU_prev = burnup;
        //cout << "  New core" << endl;

        while(kcore > 1){
            kcore_prev = kcore;

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
            }else if(mode == 0){
                core = phicalc_eqpow(core);
            } else {
                cout << endl << "Error in mode input for batch-level flux calculation(SS_burnupcalc)." << endl;
                return 0;
            }

            //disadvantage calculation
            if(DA_mode == 1){
                core = DA_calc(core);
            }

            for(int i = 0; i < N; i++){
                //!NOTE when updating this, MUST also update interpolation calcs below
                core.batch[i].Fg += core.batch[i].rflux * core.base_flux * dt;
                //cout << "  Fg: " << core.batch[i].Fg << "  rflux: " << core.batch[i].rflux << "  k: " << kcore << endl;
            }
            kcore = kcalc(core);
            //std::cout<<"kcore "<<kcore << std::endl;
            iter++;
            if(iter > 20){
                //cout << "  iter+20" << endl;
                mode = 1;
                dt *= 2;
                if(iter > 100){
                    cout << "SS_burnupcalc exceeds 100 iterations with error " << abs(burnup - BU_prev)/burnup << "%" << endl;
                    cout << "  BU: " << burnup << "  " << endl;
                    return burnup;
                }
            }
        }

        //update core fluences
        for(int i = 0; i < N; i++){
            //cout << "    " << core.batch[i].rflux;
            core.batch[i].batch_fluence = intpol(core.batch[i].Fg - (core.batch[i].rflux * core.base_flux * dt), core.batch[i].Fg, kcore_prev, kcore, 1);

            if(core.batch[i].batch_fluence > core.batch[i].collapsed_iso.fluence.back()){
                //core.batch[i].batch_fluence > core.batch[i].collapsed_iso.fluence.back() + 1;
            }
        }
        //cout << endl;
        for(ii = ii/2; core.batch[0].collapsed_iso.fluence[ii] < core.batch[0].batch_fluence; ii++){}
        burnup = intpol(core.batch[0].collapsed_iso.BU[ii-1], core.batch[0].collapsed_iso.BU[ii], core.batch[0].collapsed_iso.fluence[ii-1], core.batch[0].collapsed_iso.fluence[ii], core.batch[0].batch_fluence);
        //cout << ii << " intermed burnup: " << burnup << endl;

        if(abs(burnup - BU_prev)/burnup < core.SS_tolerance && counter > N+1){
            notsteady = false;
        }
        counter++;

        //move each batch one left, now the last and second from last are the same
        for(int i = 0; i < N-1; i++){
            core.batch[i].batch_fluence = core.batch[i+1].batch_fluence;
            core.batch[i].Fg = core.batch[i].batch_fluence;
        }

        //fix last one
        core.batch[N-1].batch_fluence = 0;
        core.batch[N-1].Fg = 0;
    }

    //cout << "SSrflux: " << core.batch[1].rflux << endl;
    //std::cout << "SSBurnupCalc " << core.name << " BU: " << burnup << "  Time: " << t.elapsed() << std::endl;
    return burnup;
}


double SS_burnupcalc_CR(fuelBundle &core, int mode, int DA_mode, double delta, int N, double ss_fluence, double target_burnup){
    //used to find the steady state CR of the given fuel
    //N:number of batches; delta: time advancement in days; PNL: nonleakage; base_flux: flux of library
    //THE FINAL BATCH IS THE OLDEST ONE
    //cout << "SS burnupcalc_CR begin" << endl;
    boost::timer t;

    isoInformation fuel = core.batch[0].collapsed_iso;
    double CR = 0;
    double burnup;
    double burnup_prev;
    double dt = delta*24*60*60; //days to [s]
    double BU, BU_est, F_est;
    double kcore_prev;
    bool notsteady = true;
    int ii = 0;
    double CR_prev = 0;
    int counter = 0;



    batch_info temp_batch;
    temp_batch.collapsed_iso = fuel;
    temp_batch.batch_fluence = 0;
    for(int i = 1; i < N; i++){
        temp_batch.Fg = ss_fluence*1/i; // use old cycles youngest batch discharge to guess
        core.batch.push_back(temp_batch);
    }

    while(notsteady){
        double kcore = 1.10000;
        int iter = 0;
        burnup_prev = burnup;
        burnup = 0;
        CR_prev = CR;
        //cout << "  New core" << endl << endl << endl;

        while(burnup < target_burnup){
            //std::cout << "Burnup " << burnup << std::endl;
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
            }else if(mode == 0){
                core = phicalc_eqpow(core);
            } else {
                cout << endl << "Error in mode input for batch-level flux calculation(SS_burnupcalc_CR)." << endl;
                return 0;
            }

            //disadvantage calculation
            if(DA_mode == 1){
                core = DA_calc(core);
            }

            for(int i = 0; i < N; i++){
                //!NOTE when updating this, MUST also update interpolation calcs below
                core.batch[i].Fg += core.batch[i].rflux * core.base_flux * dt;
                //cout << "  Fg: " << core.batch[i].Fg << "  rflux: " << core.batch[i].rflux << endl;
            }
            //cout << endl;
            for(ii = 0; core.batch[0].collapsed_iso.fluence[ii] < core.batch[0].Fg; ii++){}
            //cout << core.batch[0].collapsed_iso.BU[ii-1] << "  " << core.batch[0].collapsed_iso.BU[ii] << "  " <<  core.batch[0].collapsed_iso.fluence[ii-1];
            //cout << " " << core.batch[0].collapsed_iso.fluence[ii] << "  " <<  core.batch[0].Fg << endl;
            burnup = intpol(core.batch[0].collapsed_iso.BU[ii-1], core.batch[0].collapsed_iso.BU[ii], core.batch[0].collapsed_iso.fluence[ii-1], core.batch[0].collapsed_iso.fluence[ii], core.batch[0].Fg);
            iter++;
            if(iter > 20){
                //cout << "  iter+20" << endl;
                //mode = 1;
                //dt *= 2;
                if(iter > 100){
                    cout << "SS_burnupcalc_cr exceeds 100 iterations with error " << endl;
                    cout << "  CR: " <<  CR << "  " << endl;
                    return CR;
                }
            }
        }

        //update core fluences
        for(int i = 0; i < N; i++){
            //cout << "    " << core.batch[i].rflux;
            //std::cout << burnup_prev << " " <<  burnup << " " << target_burnup <<std::endl;
            //std::cout << core.batch[i].Fg - (core.batch[i].rflux * core.base_flux * dt) << " " << core.batch[i].Fg << std::endl;
            core.batch[i].batch_fluence = intpol(core.batch[i].Fg - (core.batch[i].rflux * core.base_flux * dt), core.batch[i].Fg, burnup_prev, burnup, target_burnup);
            //std::cout << "i:" << i << " fluence: " << core.batch[i].batch_fluence << std::endl;
            if(core.batch[i].batch_fluence > core.batch[i].collapsed_iso.fluence.back()){
                //core.batch[i].batch_fluence > core.batch[i].collapsed_iso.fluence.back() + 1;
            }
        }

        //calculate CR
        CR = CR_batch(core, 0);
        if(abs(CR - CR_prev)/CR < 0.01 && counter > N+1){
            notsteady = false;
        }
        counter++;

        //move each batch one left, now the last and second from last are the same
        for(int i = 0; i < N-1; i++){
            core.batch[i].batch_fluence = core.batch[i+1].batch_fluence;
            core.batch[i].Fg = core.batch[i].batch_fluence;
        }

        //fix last one
        core.batch[N-1].batch_fluence = 0;
        core.batch[N-1].Fg = 0;



    }
    //cout << "SSrflux: " << core.batch[1].rflux << endl;
    //std::cout << "SSBurnupCalc " << core.name << " BU: " << burnup << "  Time: " << t.elapsed() << std::endl;
    /*if(core.batch[0].batch_fluence > 5E23){
        std::cout << "OH YOU DONE MESSED UP, FLUENCE OVER 5E23 at " << core.batch[0].batch_fluence << std::endl;
    }*/
    if(CR < 0 ){CR=0;}
    return CR;

}


fuelBundle lib_interpol(fuelBundle &input_fuel){
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

    double alpha = 0.5;
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
    //std::cout << "Sum : " << met_dist_sum << std::endl;
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

void mass_check(fuelBundle &fuel){
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

double CR_finder_temp(fuelBundle &core, int ii){
    double FP = 0;
    double fissile = 0;
    double ini_fissile = 0;
    double CR;
    int ZZ;
    for(int j = 0; j < core.batch[0].collapsed_iso.iso_vector.size(); j++){
        //convert name to mass number
        ZZ = core.batch[0].collapsed_iso.iso_vector[j].name;
        ZZ = ZZ % 10000;
        ZZ /= 10;

        //add up the FP
        if(ZZ < 160 && ZZ > 40){
            //interpolation will be done at the end
            FP += core.batch[0].collapsed_iso.iso_vector[j].mass[ii];
            //cout << "ZZ " << ZZ  << " mass " << core.batch[i].collapsed_iso.iso_vector[j].mass[ii]<<endl;
        }

        //add up fissiles
        for(int fis = 0; fis < core.CR_fissile.size(); fis++){
            if(core.batch[0].collapsed_iso.iso_vector[j].name == core.CR_fissile[fis]){
                fissile += core.batch[0].collapsed_iso.iso_vector[j].mass[ii];
                ini_fissile += core.batch[0].collapsed_iso.iso_vector[j].mass[0]; //mass at fluence zero
            }
        }
    }

    if(FP > 0){
        CR = (FP+fissile-ini_fissile)/FP;
    } else {
        cout << " CR finder has fission product mass of zero." << endl;;
        CR  = 0;
    }
    //std::cout << "CR Calc " << t.elapsed() << std::endl;
    return CR;
}


void print_library(std::string name, fuelBundle &core){
    isoInformation iso = core.batch[0].collapsed_iso;
    name += ".txt";
    ofstream ofs(name);
    ofs << "Fluence";
    for(int i = 0; i < iso.neutron_prod.size(); i++){
        ofs << " " << iso.fluence[i];
    }
    ofs << std::endl << "k";
    for(int i = 0; i < iso.neutron_prod.size(); i++){
        ofs << " " << iso.neutron_prod[i]/iso.neutron_dest[i];
    }
    ofs << std::endl << "BUd";
    for(int i = 0; i < iso.BUd.size(); i++){
        ofs << " " << iso.BUd[i];
    }
    ofs << std::endl << "CR";
    for(int i = 0; i < iso.BUd.size(); i++){
        double CR_test = CR_finder_temp(core, i);
        ofs << " " << CR_test;
    }
    ofs << std::endl << std::endl;
    ofs.close();
}




