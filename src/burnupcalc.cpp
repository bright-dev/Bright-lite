#include "burnupcalc.h"
#include <boost/math/special_functions/bessel.hpp>


using namespace std;

double intpol(double y0, double y1, double x0, double x1, double x) {
    // linear interpolation function
    double y = y0 + (y1 - y0)*(x - x0)/(x1 - x0);
    return y;
}

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
        if(fuel.iso[i].type == *"A"){
            datasize = fuel.iso[i].neutron_prod.size();
            break;
        }
    }

    for(int i =0; i<fuel.iso.size(); i++){
        if(fuel.iso[i].type == *"N"){
            for(int j=0; j< datasize; j++){
                fuel.iso[i].neutron_prod.push_back(fuel.iso[i].neutron_prod[0]);
                fuel.iso[i].neutron_dest.push_back(fuel.iso[i].neutron_dest[0]);
            }
        }
    }
    return fuel;
}

fuelBundle FuelNormalizer(fuelBundle fuel){
    double actmass = 0; // total mass (or fraction) of all actinides
    int i=0;
    char *type = "A";

    while (i < fuel.iso.size()){
        if(fuel.iso[i].type == *type)
           actmass += fuel.iso[i].fraction[0]; // find the total mass of all actinides
        i++;
    }
    i=0;
    while (i < fuel.iso.size()){
        fuel.iso[i].fraction[0] = fuel.iso[i].fraction[0]/actmass; // normalize every fraction using total mass of actinides
        i++;
    }

    return fuel;

}

isoInformation regioncollapse(fuelBundle fuel, double flux){
    int i = 0, j;
    int r = 1;
    isoInformation singleiso;
    vector<isoInformation> region0, region1, region2;
    vector<isoInformation> regions;


    for(int i=0; i < fuel.iso.size(); i++){ //up to two regions at the moment

        if(fuel.iso[i].region == 0){
            region0.push_back(fuel.iso[i]);
        }
        if(fuel.iso[i].region == 1){
            region1.push_back(fuel.iso[i]);
        }
    }
    for(int i=0; i<region0.size(); i++){ //uses the flux to adjust prod and dest for region0, the fuel
        for(int j =0; j < region0[i].neutron_prod.size(); j++){
            region0[i].neutron_prod[j] *= flux;
            region0[i].neutron_dest[j] *= flux;
        }
    }

    /*for (int i = 0; i < region1.size();i++){
        cout << region1[i].name << "    "<<region1[i].fraction << endl;
    }*/
    regions.push_back(FuelBuilder(region0));
    regions.push_back(FuelBuilder(region1));
    regions[0].fraction.push_back(1.);
    regions[1].fraction.push_back(1.);
    return FuelBuilder(regions);
};
// need to add units

map<int, double> tomass (int ti, double fluence, isoInformation isoinfo) {
    map<int, double> out = map<int, double>();
    double mass_i;
    int name_i;
    int nucid;
    for (int i = 0; i < isoinfo.iso_vector.size(); i++){
        name_i = isoinfo.iso_vector[i].name;
        nucid = pyne::nucname::zzaaam(name_i)/10;
        mass_i = intpol(isoinfo.iso_vector[i].mass[ti-1],
                        isoinfo.iso_vector[i].mass[ti],
                        isoinfo.fluence[ti-1],
                        isoinfo.fluence[ti],
                        fluence);
        out[nucid] = mass_i;
    }
    return out;
}


/**
phicalc calculates the flux of a given batch
inputs: n is the batch number starting from zero, total bathces N, discharge burnup BU_total
        isoInformation of the fuel
**/

double phicalc(int n, int N, double fluence, isoInformation tempone){
// batch n of N, batch n at fluence level fluence, batch n is type tempone
    int i = 0;
    double x0=0, x1, dF_n, dB_n, dt_n;
    double BU_total = 0;

    while (tempone.fluence[i] < fluence) // finds the discrete point i corresponding to the fluence just under fluence
        {
            x0 = tempone.fluence[i];
            i++;
        }

    if(i == 0){
        dF_n = tempone.fluence[1];
        dt_n = N * 180;
    }else{

        x1 = tempone.fluence[i];
        dF_n = fluence - x0; //delta fluence
        for(int j = 0; j < i; j++){
            BU_total += tempone.BUd[j];
        }

        dB_n = intpol(BU_total, BU_total+tempone.BUd[i], x0, x1, fluence) - BU_total;
        dt_n = (N * dB_n*180) / BU_total; //delta time due to the change in fluence

    }

    return dF_n/dt_n; //flux of n'th batch, n indexed from zero
}

double SSphicalc(int n, int N, double BU_total, isoInformation tempone){
    int i = 0;
    double x0=0, x1, F_n, dF_n, dB_n, dt_n;
    double BU_n;
    BU_n = BU_total*(n+1)/N;

    while (x0 + tempone.BUd[i] < BU_n) // finds the discrete point i corresponding to the burnup just under BU_n
        {
            x0 += tempone.BUd[i]; //x0 is the total burnup under BU_n
            i++;
        }

    x1 = x0 + tempone.BUd[i]; // adds on more discrete point for linear interpolation

    F_n = intpol(tempone.fluence[i-1],tempone.fluence[i],x0,x1,BU_n); //fluence at BU_n
    dF_n = (tempone.fluence[i] - F_n); //delta fluence
    dB_n = x1 - BU_n; //delta burnup due to the delta fluence
    dt_n = (N * dB_n*180) / BU_total; //delta time due to the change in fluence
//cout<< "dF_n: "<< dF_n<<" dB_n: "<<dB_n<<" dt_n: "<<dt_n<<endl;
    return dF_n/dt_n; //flux of n'th batch, n indexed from zero

}

double SSkcalc(isoInformation tempone, double BU_total, int N, double pnl){
    double x0 = 0, x1 = 0, BU_n = 0, phi, pbatch[N], dbatch[N], p_total=0, d_total=0;
    int j = 0, i = 0;

            for (j = 0; j < N; j++) {
                BU_n = BU_total*(j+1)/N;
                x0 = 0; //these variables are reset here b/c they're used outside the for-loop as well
                x1= 0;
                i=0;
                while (x0 + tempone.BUd[i] < BU_n) // finds the discrete point i corresponding to the burnup just under BU_n
                    {
                        x0 += tempone.BUd[i];
                        i++;
                    }
                x1 = x0 + tempone.BUd[i]; // adds on more discrete point for linear interpolation

                phi = SSphicalc(j, N, BU_total, tempone);
                pbatch[j] = intpol(tempone.neutron_prod[i-1],tempone.neutron_prod[i], x0, x1, BU_n);
                dbatch[j] = intpol(tempone.neutron_dest[i-1],tempone.neutron_dest[i], x0, x1, BU_n);

                p_total += pbatch[j]*phi;
                d_total += dbatch[j]*phi;
            }
    ofstream outfile("../outputfile.txt");

    outfile<< N << " "<< BU_total << endl << endl;

    j=0; // recycled variable, is also batch index in this function
    while(j < tempone.iso_vector.size()){
        //i is still the discrete point of the last bach
        outfile << tempone.iso_vector[j].name << " ";
        outfile << intpol(tempone.iso_vector[j].mass[i-1], tempone.iso_vector[j].mass[i], x0, x1, BU_n) << endl;
        j++;

    }

    outfile.close();

    return (p_total*pnl)/(d_total);

}

double kcalc(std::vector<isoInformation> isoBatches, std::vector<double> batch_fluence, double fluence, std::vector<double> batch_phi, double pnl){
//takes isoInformation of all batches and returns the k at the fluence level
//starting fluence levels are saved in the isoInformation structure
//the variable fluence passed to this function is how much the bundles will be burned

    double total_prod = 0;
    double total_dest = 0;
    int j;


    for(int i = 0; i < isoBatches.size(); i++){
        j = 0;
        if(batch_fluence[i] + fluence*batch_phi[i] != 0){
            while(isoBatches[i].fluence[j] < batch_fluence[i] + fluence*batch_phi[i]){
                //finds the discrete point where fluence is less than the target fluence
                j++;
            }
            total_prod += intpol(isoBatches[i].neutron_prod[j-1], isoBatches[i].neutron_prod[j],
                                 isoBatches[i].fluence[j-1], isoBatches[i].fluence[j], batch_fluence[i] + fluence*batch_phi[i]);
            total_dest += intpol(isoBatches[i].neutron_dest[j-1], isoBatches[i].neutron_dest[j],
                                 isoBatches[i].fluence[j-1], isoBatches[i].fluence[j], batch_fluence[i] + fluence*batch_phi[i]);
        }
        else{
            total_prod += isoBatches[i].neutron_prod[0];
            total_dest += isoBatches[i].neutron_dest[0];
        }
    }

    if(total_dest == 0){
        cout << "Error in kcalc! Total neutron destruction is zero."<<endl;
    }
    return total_prod*pnl/total_dest;

}

pair<double, map<int, double> > burnupcalc(vector<fuelBundle> batches, double pnl, double tolerance, double flux) {

    pair<double, map<int,double> > rtn(0, map<int, double>());
    double F1, F2, F3, k_total = 0;
    double fluence = 0; //added fluence to all the batches due to burnup
    std::vector<isoInformation> isoBatches;
    std::vector<double> batch_fluence;
    std::vector<double> batch_phi;
    double phimax = 0;
    int N = batches.size();
    double burnup = 0;
    int oldest_batch = N; //index for oldest batch, index starts from zero

    //collapse all the baches to isoInformation and build batch_fluence vector
    for(int i = 0; i < batches.size(); i++){
        isoBatches.push_back(regioncollapse(batches[i], 1));
        batch_fluence.push_back(batches[i].batch_fluence);
        batch_phi.push_back(phicalc(i, N, batch_fluence[i], isoBatches[i]));
        if(batch_phi[i]>phimax){
            phimax = batch_phi[i];
        }
        if(batches[i].batch_fluence > oldest_batch){
            oldest_batch = i;
        }
    }

    for(int i = 0; i < batch_phi.size(); i++){
        //normalize batch_phi using max value
        batch_phi[i] = batch_phi[i]/phimax;
    }



    if(kcalc(isoBatches, batch_fluence, fluence, batch_phi, pnl) < 1){
        cout << "Error! Original core configuration is not critical." << endl;
    }

F1 = 0;
F2 = 100;

int i = 0;
    while(abs(1-k_total)>tolerance){

        F3 = F2 - (kcalc(isoBatches, batch_fluence, F2, batch_phi, pnl)-1)*(F2-F1)/
                        ((kcalc(isoBatches, batch_fluence, F2, batch_phi, pnl))-(kcalc(isoBatches, batch_fluence, F1, batch_phi, pnl)));
        k_total = kcalc(isoBatches, batch_fluence, F3, batch_phi, pnl);
        F1 = F2;
        F2 = F3;
        i++;

        if(i==50){
            cout<< "Warning! Maximum iteration reached in function burnupcalc."<<endl;
            F3 = (F1+F2+F3)/3;
            break;
        }
    }//F3 is the fluence of the cycle


    for(int j = 0; isoBatches[oldest_batch].fluence[j] < F3*batch_phi[oldest_batch]; j++){
        burnup += isoBatches[oldest_batch].BUd[j];
    }
    cout << "Burnup: "<<burnup << "    k at this burnup: "<< k_total << endl;

    rtn.first = burnup;
    //rtn.second = tomass(i, time_f, tempone);
    return rtn;

}

pair<double, map<int, double> > SSburnupcalc(fuelBundle fuel, int N, double pnl, double tolerance, double flux) {
    //Stead State (SS) burnupcalc, takes only the starting bundle
    pair<double, map<int,double> > rtn(0, map<int, double>());
    double time_f; // time when k reaches one, time in days
    double BU_total = 0;
    double k_total = 10;
    double BU1, BU2, BU3;
    int i_stor;

    isoInformation tempone = regioncollapse(fuel, 0.95);

    for (int i = 0; i < tempone.neutron_prod.size(); i++){
        tempone.k_inf.push_back(((tempone.neutron_prod[i]))*pnl/(tempone.neutron_dest[i]));
        cout << tempone.k_inf[i] << endl;
    }

    if (tempone.k_inf[1] < 1){
        cout << "Core setup not critical" << endl;
        return rtn;
    }
    for (int i = 1; i < tempone.k_inf.size(); i++){
        BU_total += tempone.BUd[i];
        if (tempone.k_inf[i] < 1.0){
            BU_total = intpol(BU_total,BU_total+tempone.BUd[i],tempone.k_inf[i-1],tempone.k_inf[i],1.0);
            //cout << BU_total << endl;
            break;
        }
        if (i == tempone.k_inf.size() -1. ){
            cout << "Fluence not high enough to bring core subcritical" << endl;
        }
        i_stor = i;
    }
    time_f = intpol(tempone.fluence[i_stor-1],tempone.fluence[i_stor],tempone.k_inf[i_stor-1],tempone.k_inf[i_stor],1.0);
    BU1 = 2.* N * BU_total/(N+1);
    BU2 = BU1*1.1;
    k_total = 2;
    int i=0;
    while(abs(1-k_total)>tolerance){
        BU3 = BU2 - (SSkcalc(tempone, BU2, N, pnl)-1)*(BU2-BU1)/((SSkcalc(tempone, BU2, N, pnl))-(SSkcalc(tempone, BU1, N, pnl)));
        k_total = SSkcalc(tempone, BU3, N, pnl);
        BU1 = BU2;
        BU2 = BU3;
        i++;

        if(i==50){
            cout<< "Warning! Maximum iteration reached."<<endl;
            BU3 = (BU1+BU2+BU3)/3;
            break;
        }
    }
    BU_total = BU3;

    rtn.first = BU_total;
    rtn.second = tomass(i_stor, time_f, tempone);
    return rtn;

}

fuelBundle enrich_collapse(fuelBundle fuel){
    for(int i = 0; i < fuel.iso.size(); i++){
        fuel.iso[i].fraction[0] = 0;
        for(int j = 1; j < fuel.iso[i].fraction.size(); j++){
            fuel.iso[i].fraction[0] += fuel.iso[i].fraction[j] * fuel.stream_fraction[j-1];
        }
    }
}

fuelBundle fluxcalc_reader(fuelBundle fuel, string file_name){
    //must be called AFTER InputReader, so that fuelbundle is built
    //assumes the format [type of material][space]sigma_scatter[space]sigma_a in the input file
    string line;
    int nucid;
    double sigs, siga;
    double x;
    char name[10];

    ifstream fin(file_name);


	while(getline(fin, line))
	{
        if(line.find("a") == 0){
            istringstream iss(line);
            iss >> name >> x;
            fuel.fuel_radius = x;
            continue;
        }

        if(line.find("b") == 0){
            istringstream iss(line);
            iss >> name >> x;
            fuel.moderator_radius = x;
            continue;
        }

        if(line.find("MODERATOR") == 0){
            istringstream iss(line);
            iss >> name >> sigs >> siga;
            fuel.moderator_sigs = sigs;
            fuel.moderator_siga = siga;
            continue;
        }

        istringstream iss(line);
        iss >> nucid >> sigs >> siga;

        if(nucid){
            for(int i = 0; i < fuel.iso.size(); i++){
                if(fuel.iso[i].name == nucid){
                    fuel.iso[i].sigs = sigs;
                    fuel.iso[i].siga = siga;
                    fuel.iso[i].fuel = true;
                }
            }
        }
	}


    return fuel;
}

double fluxcalc(fuelBundle fuel){
// calculates the flux of each region in fuelBundle
// probably will need to add reactor identifier as input in the future

    int i = 0;
    int r = 1; //number of regions
    double frac35, frac38; // mass fraction of u235 and u238
    double temp;


    double a = fuel.fuel_radius; // radius of the fuel rod
    double b = fuel.moderator_radius; // radius of the equivalent cell
    double L_F; // diffusion length of fuel
    double L_M; // diffusion length of moderator
    double Sig_aF; // macroscopic abs. CS of fuel
    double Sig_aM = fuel.moderator_siga; // macroscopic abs. CS of moderator
    double V_F; // volume of fuel
    double V_M; // volume of moderator
    double Sig_trF; // macroscopic transport CS of fuel
    double Sig_trM; // macroscopic transport CS of moderator
    double Sig_tF; // macroscopic total CS of fuel
    double Sig_tM; //macroscopic total CS of moderator
    double Sig_sF; // macroscopic scatter CS of fuel
    double Sig_sM = fuel.moderator_sigs; //macroscopic scatter CS of moderator
    double D_F; // diffusion coef. of fuel
    double D_M; // diffusion coef. of moderator
    double A_F; // A number of fuel
    double A_M; // A number of moderator
    double x, y, z; // calculated equivalent dimensions
    double F, E; // lattice functions
    double f; // flux of fuel divided by total flux(fuel+moderator)


/**************moderator****************/
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
/****************************************/


    vector<int> fuel_index;
    for(int i = 0; i < fuel.iso.size(); i++){
        if(fuel.iso[i].fuel == true){
            fuel_index.push_back(i);
        }
    }

    for(int fluence = 0; fluence < fuel.iso[0].neutron_dest.size(); fluence++){
        Sig_aF = 0;
        Sig_sF = 0;

        for(int i = 0; i < fuel_index.size(); i++){
            cout << fuel.iso[fuel_index[i]].sigs << endl;
            Sig_aF += (fuel.iso[fuel_index[i]].neutron_dest[fluence] * fuel.iso[fuel_index[i]].fraction[1]/100);
            Sig_sF += fuel.iso[fuel_index[i]].sigs * fuel.iso[fuel_index[i]].fraction[1];
        }

        Sig_tF = Sig_aF+Sig_sF;
        Sig_trF = Sig_tF - 2/3/A_F*Sig_sF;
        D_F = 1 / (3 * Sig_trF);
        L_F = sqrt(D_F/Sig_aF);
        x = a/L_F;

        /*****book example***
        a=1.02;
        b=14.3;
        x=0.658;
        y=0.0173;
        z=0.242;
        V_M=195.6;
        V_F=1;
        Sig_aM=0.0002728;
        Sig_aF=0.3668;
        *******************/

        F = x * boost::math::cyl_bessel_i(0,x) / (2 * boost::math::cyl_bessel_i(1, x));
        E = (z*z - y*y) / (2 * y) * ( (boost::math::cyl_bessel_i(0, y) * boost::math::cyl_bessel_k(1, z)+ boost::math::cyl_bessel_k(0, y) *
                                       boost::math::cyl_bessel_i(1, z)) / (boost::math::cyl_bessel_i(1, z) *
                                        boost::math::cyl_bessel_k(1, y) - boost::math::cyl_bessel_k(1, z) * boost::math::cyl_bessel_i(1, y)));
        f = pow((((Sig_aM * V_M)/(Sig_aF * V_F)) * F + E), (-1.));
        cout << "Disadvtg: " << (Sig_aF*V_F - f*Sig_aF*V_F)/(f*Sig_aM*V_M)<<endl;

    }

    return (Sig_aF*V_F - f*Sig_aF*V_F)/(f*Sig_aM*V_M);
}


/**
enrichcalc takes the burnup goal, batch number, and tolerance;
and returns the enrichment needed to reach the burnup goal.
The code initially estimates the enrichment needed by
calculating the burnup (for a given number of batches and type
of reactor) of a two and seven percent enriched fuel and then
interpolates/extrapolates from these two points. The code then
calculates the burnup for this estimated enrichment level, and
iterates by changing the guess enrichment until the desired
burnup is achieved by the guess value. This enrichment is then
returned.


Inputs

BU_end: The target burnup. The code finds an enrichment which
will result in a burnup equal to this value.

N: number of batches used in the core.

tolerance: Currently not utilized, used to adjust precision.


Outputs:

Double, the enrichment needed to achieve BU_end.

*/

pair<double, pair<double, map<int,double> > > enrichcalc(double BU_end, int N, double tolerance, fuelBundle fuel) {
    pair<double, pair<double, map<int, double> > > rtn;
    double X;
    double BU_guess, enrich_guess;
    double enrich_previous, BU_previous;
    double enrich_lower, enrich_upper;
    double BU_lower, BU_upper;
    vector <int> blend_vector;

    /** This is a super quick hack to benchmark enrichcalc */
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
        BU_lower = SSburnupcalc(fuel, fuel.batch, fuel.pnl, 0.01, 1).first;
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
        BU_upper = SSburnupcalc(fuel1, fuel.batch, fuel.pnl, 0.01, 1).first;
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
    BU_guess = SSburnupcalc(fuel2, fuel.batch, fuel.pnl, 0.01, 1).first;
    BU_guess = SSburnupcalc(fuel2, fuel.batch, fuel.pnl, 0.01, 1).first;
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
        BU_guess = SSburnupcalc(fuel3, fuel.batch, fuel.pnl, 0.01, 1).first;
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
    rtn.second = SSburnupcalc(fuel3, fuel.batch,fuel.pnl, 0.01, 1);
    return rtn;
}

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
                fuel.iso.push_back(temp);
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

fuelBundle lib_interpol(fuelBundle input_fuel){
    /// TODO Add NONA to these libs before combining ///
    vector<fuelBundle> fuel_pairs;
    for (int i = 0; i < input_fuel.interpol_libs.size(); i++){
        fuelBundle lib_bundle;
        for(int j = 0; j < input_fuel.iso.size(); j++){
            isoInformation iso;
            iso.name = input_fuel.iso[j].name;
            iso.fraction = input_fuel.iso[j].fraction;
            iso.type = input_fuel.iso[j].type;
            iso.blending = input_fuel.iso[j].blending;
            iso.region = input_fuel.iso[j].region;
            lib_bundle.iso.push_back(iso);
        }
        lib_bundle.name = input_fuel.interpol_libs[i];
        lib_bundle.batch = input_fuel.batch;
        lib_bundle.tres = input_fuel.tres;
        lib_bundle.pnl = input_fuel.pnl;
        lib_bundle.target_BUd = input_fuel.target_BUd;
        fuel_pairs.push_back(lib_bundle);
    }
    for (int i = 0; i < fuel_pairs.size(); i++){
        DataReader2(fuel_pairs[i].name, fuel_pairs[i].iso);
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
        double distance_measure = 0;
        for (int j = 0; j < metrics.size(); j++){
            distance_measure += pow(input_fuel.interpol_pairs[j].scaled_value - metrics[j][i], 2);
        }
        if(distance_measure == 0){
            return fuel_pairs[i];
        }
        metric_distances.push_back(pow(distance_measure, alpha/2));
    }
    double met_dist_sum;
    for (int i = 0; i < metric_distances.size(); i++){
        met_dist_sum += 1./metric_distances[i];
    }
    // Fuel Bundle instead of iso //
    fuelBundle new_fuel;
    new_fuel.tres = input_fuel.tres;
    new_fuel.pnl = input_fuel.pnl;
    new_fuel.batch = input_fuel.batch;
    new_fuel.operation_type = input_fuel.operation_type;
    new_fuel.target_BUd = input_fuel.target_BUd;
    for(int i = 0; i < fuel_pairs.size(); i++){
        if(i==0){
            for(int j = 0; j < fuel_pairs[i].iso.size(); j++){
                isoInformation new_iso;
                new_iso.name = fuel_pairs[i].iso[j].name;
                new_iso.region = fuel_pairs[i].iso[j].region;
                new_iso.blending = fuel_pairs[i].iso[j].blending;
                new_iso.fraction = fuel_pairs[i].iso[j].fraction;
                new_iso.type = fuel_pairs[i].iso[j].type;
                for (int k = 0; k < fuel_pairs[i].iso[j].fluence.size(); k++){
                    new_iso.fluence.push_back(fuel_pairs[i].iso[j].fluence[k]);
                }
                for (int k = 0; k < fuel_pairs[i].iso[j].neutron_prod.size(); k++){
                    new_iso.neutron_prod.push_back(fuel_pairs[i].iso[j].neutron_prod[k]/metric_distances[i]/met_dist_sum);
                }
                for (int k = 0; k < fuel_pairs[i].iso[j].neutron_dest.size(); k++){
                    new_iso.neutron_dest.push_back(fuel_pairs[i].iso[j].neutron_dest[k]/metric_distances[i]/met_dist_sum);
                }
                for (int k = 0; k < fuel_pairs[i].iso[j].BUd.size(); k++){
                    new_iso.BUd.push_back(fuel_pairs[i].iso[j].BUd[k]/metric_distances[i]/met_dist_sum);
                }
                for (int k = 0; k < fuel_pairs[i].iso[j].iso_vector.size(); k++){
                    daughter new_daughter;
                    new_daughter.name = fuel_pairs[i].iso[j].iso_vector[k].name;
                    for(int ii = 0; ii < fuel_pairs[i].iso[j].iso_vector[k].mass.size(); ii++){
                        new_daughter.mass.push_back(fuel_pairs[i].iso[j].iso_vector[k].mass[ii]/metric_distances[i]/met_dist_sum);
                    }
                }
                new_fuel.iso.push_back(new_iso);
            }
        } else {
            for(int j = 0; j < fuel_pairs[i].iso.size(); j++){
                bool iso_check1 = false;
                for(int m = 0; m < new_fuel.iso.size(); m++){
                    if(fuel_pairs[i].iso[j].name == new_fuel.iso[m].name){
                        iso_check1 = true;
                        for (int k = 0; k < fuel_pairs[i].iso[j].neutron_prod.size(); k++){
                            new_fuel.iso[m].neutron_prod[k] += fuel_pairs[i].iso[j].neutron_prod[k]/metric_distances[i]/met_dist_sum;
                        }
                        for (int k = 0; k < fuel_pairs[i].iso[j].neutron_dest.size(); k++){
                            new_fuel.iso[m].neutron_dest[k] += fuel_pairs[i].iso[j].neutron_dest[k]/metric_distances[i]/met_dist_sum;
                        }
                        for (int k = 0; k < fuel_pairs[i].iso[j].BUd.size(); k++){
                            new_fuel.iso[m].BUd[k] += fuel_pairs[i].iso[j].BUd[k]/metric_distances[i]/met_dist_sum;
                        }
                        for (int k = 0; k < fuel_pairs[i].iso[j].iso_vector.size(); k++){
                            bool iso_check2 = false;
                            for (int mm = 0; mm < new_fuel.iso[m].iso_vector.size(); mm++){
                                if(fuel_pairs[i].iso[j].iso_vector[k].name == new_fuel.iso[m].iso_vector[mm].name){
                                    iso_check2 = true;
                                    for(int n = 0; n < new_fuel.iso[m].iso_vector[mm].mass.size(); n++){
                                        new_fuel.iso[m].iso_vector[mm].mass[n] += fuel_pairs[i].iso[j].iso_vector[k].mass[n]/metric_distances[i]/met_dist_sum;
                                    }
                                }
                            }
                            if (iso_check2 == false){
                                daughter new_daught;
                                new_daught.name = fuel_pairs[i].iso[j].iso_vector[k].name;
                                for(int n = 0; n < fuel_pairs[i].iso[j].iso_vector[k].mass.size(); n++){
                                    new_daught.mass.push_back(fuel_pairs[i].iso[j].iso_vector[k].mass[n]/metric_distances[i]/met_dist_sum);
                                }
                                new_fuel.iso[m].iso_vector.push_back(new_daught);
                            }
                        }
                    }
                } if (iso_check1 == false){
                    isoInformation new_iso;
                    new_iso.name = fuel_pairs[i].iso[j].name;
                    new_iso.region = fuel_pairs[i].iso[j].region;
                    new_iso.blending = fuel_pairs[i].iso[j].blending;
                    new_iso.fraction = fuel_pairs[i].iso[j].fraction;
                    new_iso.type = fuel_pairs[i].iso[j].type;
                    for (int k = 0; k < fuel_pairs[i].iso[j].fluence.size(); k++){
                        new_iso.fluence.push_back(fuel_pairs[i].iso[j].fluence[k]);
                    }
                    for (int k = 0; k < fuel_pairs[i].iso[j].neutron_prod.size(); k++){
                        new_iso.neutron_prod.push_back(fuel_pairs[i].iso[j].neutron_prod[k]/metric_distances[i]/met_dist_sum);
                    }
                    for (int k = 0; k < fuel_pairs[i].iso[j].neutron_dest.size(); k++){
                        new_iso.neutron_dest.push_back(fuel_pairs[i].iso[j].neutron_dest[k]/metric_distances[i]/met_dist_sum);
                    }
                    for (int k = 0; k < fuel_pairs[i].iso[j].BUd.size(); k++){
                        new_iso.BUd.push_back(fuel_pairs[i].iso[j].BUd[k]/metric_distances[i]/met_dist_sum);
                    }
                    for (int k = 0; k < fuel_pairs[i].iso[j].iso_vector.size(); k++){
                        daughter new_daughter;
                        new_daughter.name = fuel_pairs[i].iso[j].iso_vector[k].name;
                        for(int ii = 0; ii < fuel_pairs[i].iso[j].iso_vector[k].mass.size(); ii++){
                            new_daughter.mass.push_back(fuel_pairs[i].iso[j].iso_vector[k].mass[ii]/metric_distances[i]/met_dist_sum);
                        }
                    }
                    new_fuel.iso.push_back(new_iso);
                }
            }
        }
    }
    return new_fuel;
}

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
    */
    return 0;
} 

*/
