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
                total_prod += nona[j].total_prod*1000*0.602/name*fuel.iso[i].fraction;
                fuel.iso[i].neutron_dest.push_back(nona[j].total_dest*1000*0.602/name);
                total_dest += nona[j].total_dest*1000*0.602/name*fuel.iso[i].fraction;
            }
        }
    }
    //cout << total_prod << "     " << total_dest << endl;
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
           actmass += fuel.iso[i].fraction; // find the total mass of all actinides
            i++;
        }

    i=0;
    while (i < fuel.iso.size()){
        fuel.iso[i].fraction = fuel.iso[i].fraction/actmass; // normalize every fraction using total mass of actinides
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
    regions[0].fraction = 1;
    regions[1].fraction = 1;
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

double phicalc(int n, int N, double BU_total, isoInformation tempone){
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

double kcalc(isoInformation tempone, double BU_total, int N, double pnl){
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

                phi = phicalc(j, N, BU_total, tempone);
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

pair<double, map<int, double> > burnupcalc(fuelBundle fuel, int N, double pnl, double tolerance, double flux) {
    isoInformation tempone = regioncollapse(fuel, flux);
    pair<double, map<int,double> > rtn(0, map<int, double>());
    double time_f; // time when k reaches one, time in days
    double BU_total = 0;
    double k_total = 10;
    int i = 0;
    double BU1, BU2, BU3;


    for (int i = 0; i < tempone.neutron_prod.size(); i++){
        tempone.k_inf.push_back(((tempone.neutron_prod[i]))*pnl/(tempone.neutron_dest[i]));
    }

    if (tempone.k_inf[0] < 1){
        cout << "Core setup not critical" << endl;
        return rtn;
    }
    for (int i = 0; i < tempone.k_inf.size(); i++){
        BU_total += tempone.BUd[i];
        if (tempone.k_inf[i] < 1.0){
            BU_total = intpol(BU_total,BU_total+tempone.BUd[i],tempone.k_inf[i-1],tempone.k_inf[i],1.0);
            break;
        }
        if (i == tempone.k_inf.size() -1. ){
            cout << "Fluence not high enough to bring core subcritical" << endl;
        }
    }
    time_f = intpol(tempone.fluence[i-1],tempone.fluence[i],tempone.k_inf[i-1],tempone.k_inf[i],1.0);
    BU1 = 2.* N * BU_total/(N+1);
    BU2 = BU1*1.1;
    k_total = 2;
    i=0;

    while(abs(1-k_total)>tolerance){
        BU3 = BU2 - (kcalc(tempone, BU2, N, pnl)-1)*(BU2-BU1)/((kcalc(tempone, BU2, N, pnl))-(kcalc(tempone, BU1, N, pnl)));
        k_total = kcalc(tempone, BU3, N, pnl);
        BU1 = BU2;
        BU2 = BU3;
        i++;

        if(i==50)
            {
            cout<< "Warning! Maximum iteration reached."<<endl;
            BU3 = (BU1+BU2+BU3)/3;
            break;
            }



    }
    BU_total = BU3;

    rtn.first = BU_total;
    rtn.second = tomass(i, time_f, tempone);
    return rtn;

}

double fluxcalc(fuelBundle fuel){
// calculates the flux of each region in fuelBundle
// probably will need to add reactor identifier as input in the future

    int i = 0;
    int r = 1; //number of regions
    double frac35, frac38; // mass fraction of u235 and u238
    double temp;

    while(i < fuel.iso.size()){ // finds the number of regions and saves mass fraction of fuel
        if(fuel.iso[i].region > r)
            r = fuel.iso[i].region;
        if(fuel.iso[i].name == 922350)
            frac35 = fuel.iso[i].fraction;
        if(fuel.iso[i].name == 922380)
            frac38 = fuel.iso[i].fraction;
        i++;
    }

    //cout << "235: " << frac35 << endl << "238: " << frac38 << endl;


    double flux[r]; // creates a flux vector


    if (r == 1){ //r=1 means two regions
        double a; // radius of the fuel rod
        double b; // radius of the equivalent cell
        double L_F; // diffusion length of fuel
        double L_M; // diffusion length of moderator
        double Sig_aF; // macroscopic abs. CS of fuel
        double Sig_aM; // macroscopic abs. CS of moderator
        double V_F; // volume of fuel
        double V_M; // volume of moderator
        double Sig_trF; // macroscopic transport CS of fuel
        double Sig_trM; // macroscopic transport CS of moderator
        double Sig_tF; // macroscopic total CS of fuel
        double Sig_tM; //macroscopic total CS of moderator
        double Sig_sF; // macroscopic scatter CS of fuel
        double Sig_sM; //macroscopic scatter CS of moderator
        double D_F; // diffusion coef. of fuel
        double D_M; // diffusion coef. of moderator
        double A_F; // A number of fuel
        double A_M; // A number of moderator
        double x, y, z; // calculated equivalent dimensions
        double F, E; // lattice functions
        double f; // flux of fuel divided by total flux(fuel+moderator)

        temp = frac35;
        frac35 = frac35 / (frac35 + frac38);
        frac38 = frac38 / (temp + frac38);

        double abs35, sca35, tot35; //xsecs for u235
        double abs38, sca38, tot38; //xsecs for u238

        //abs35 = 608.4-14.95;
        abs35 = 608.4-14.95;
        sca35 = 14.95;
        tot35 = 608.4;

        //abs38 = 11.77-9.356;
        abs38 = 11.77-9.356;
        sca38 = 9.360;
        tot35 = 11.77;
        Sig_aF = abs35*frac35 + abs38*frac38;
        Sig_aM = 0.000094*pow(10,1);
        a = 0.4095; // [cm]
        b = 0.70749;// [cm]

    // transport CS calculation
        Sig_tF = (tot35*frac35 + tot38*frac38)*pow(10,1); // [cm]
        Sig_tM = 2.75*pow(10,1); // [cm]
        Sig_sF = (sca35*frac35+sca38*frac38)*pow(10,1); // [cm]
        Sig_sM = 2.739*pow(10,1); // [cm]
        A_F = 235;
        A_M = 18;
        Sig_trF = Sig_tF - 2/3/A_F*Sig_sF;
        Sig_trM = Sig_tM - 2/3/A_M*Sig_sM;

    // diffusion calculation
        D_F = 1 / (3 * Sig_trF);
        D_M = 1 / (3 * Sig_trM);

    // diffusion length calculation
        L_F = sqrt(D_F/Sig_aF);
        L_M = sqrt(D_M/Sig_aM);

        x = a/L_F;
        y = a/L_M;
        z = b/L_M;
        V_M = pow(a,2)*3.141592;
        V_F = pow(b,2)*3.141592 - pow(a,2)*3.141592;


        F = x * boost::math::cyl_bessel_i(0,x) / (2 * boost::math::cyl_bessel_i(0, x));
        E = (z*z - y*y) / (2 * y) * ( (boost::math::cyl_bessel_i(0, y) * boost::math::cyl_bessel_k(1, z)+ boost::math::cyl_bessel_k(0, y) *
                                       boost::math::cyl_bessel_i(1, z)) / (boost::math::cyl_bessel_i(1, z) *
                                        boost::math::cyl_bessel_k(0, y) - boost::math::cyl_bessel_k(1, z) * boost::math::cyl_bessel_i(0, y)));
        f = pow((((Sig_aM * V_M)/(Sig_aF * V_F)) * F + E), (-1.));

        flux[1] = 1;
        flux[0] = f / (f - 1);
    } else{
        i=0;
        while(i < r){
            flux[i] = 1;
            i++;
        }

    }


    return flux[0];

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

double enrichcalc(double BU_end, int N, double tolerance, fuelBundle fuel, double flux) {
    //THIS WORKS ONLY FOR TWO ISO'S
    double X;
    double BU_guess, enrich_guess;
    double enrich_previous, BU_previous;
    double enrich_lower, enrich_upper;
    double BU_lower, BU_upper;
    vector <int> blend_vector;

    /** This is a super quick hack to benchmark enrichcalc */
    fuelBundle fuel1;
    fuelBundle fuel2;
    for (int i = 0; i < fuel.iso.size(); i++){
        if (fuel.iso[i].blending == true){
            blend_vector.push_back(i);
        }
    }
    for (int i = 0; i < 100; i++){
        X = i / 100.;
        fuel.iso[blend_vector[0]].fraction = X;
        fuel.iso[blend_vector[1]].fraction = 1 - X;
        fuel = FuelNormalizer(fuel);
        BU_lower = burnupcalc(fuel, fuel.batch, fuel.pnl, 0.01, flux).first;
        if (BU_lower > 0){
            enrich_lower = X;
            break;
        }
    }

    for (int i = enrich_lower*100; i < enrich_lower*100+10; i++){
        X = i / 100.;
        fuel.iso[blend_vector[0]].fraction = X;
        fuel.iso[blend_vector[1]].fraction = 1 - X;
        fuel1 = FuelNormalizer(fuel);
        BU_upper = burnupcalc(fuel1, fuel.batch, fuel.pnl, 0.01, flux).first;
        if (i > 99){
            cout << "IT'S ALL BROKEN" << endl;
            return 0;
        }
        if (BU_upper > 2*BU_lower){
            enrich_upper = X;
            break;
        }
    }
    X = enrich_lower + (BU_end - BU_lower)*(enrich_upper - enrich_lower)/(BU_upper - BU_lower);
    if (BU_upper == 0){
        cout << "Burn up code failed" << endl;
        return 0;
    }

    /** rest of the haxxor */
    fuel.iso[blend_vector[0]].fraction = X;
    fuel.iso[blend_vector[1]].fraction = 1-X;
    fuel2 = FuelNormalizer(fuel);
    BU_guess = burnupcalc(fuel2, fuel.batch, fuel.pnl, 0.01, flux).first;
    if (BU_guess == 0){
        cout << "Burn up code failed" << endl;
        return 0;
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

    while ((abs(BU_end - BU_guess)/BU_end) > 0.01){
        cout << enrich_previous << "   " << BU_previous << endl;
        cout << enrich_guess << "   " << BU_guess << endl;
        cout << "TEST" << endl;
        fuelBundle fuel3;
        X = enrich_previous + (BU_end - BU_previous)*(enrich_guess - enrich_previous)/(BU_guess - BU_previous);
        BU_previous = BU_guess;
        enrich_previous = enrich_guess;
        fuel.iso[blend_vector[0]].fraction = X;
        fuel.iso[blend_vector[1]].fraction = 1. - X;
        fuel3 = FuelNormalizer(fuel);
        BU_guess = burnupcalc(fuel3, fuel.batch, fuel.pnl, 0.01, flux).first;
        enrich_guess = X;
    }
    return X;
}

fuelBundle InputReader(){
    std::string name, mass;
    int region, N, t_res;
    char type, word[8];
    int nucid;
    double pnl, intpol_val, target_BUd;
    fuelBundle fuel;
    isoInformation temp;

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
        if(line.find("REGIONS") == 0){
            while(getline(fin, line)){
                if(line.find("END") == 0) break;
                istringstream iss(line);
                iss >> region >> type >> nucid >> mass;
                temp.name = nucid;
                temp.region = region;
                temp.type = type;
                if (mass == "X"){
                    temp.blending = true;
                    temp.fraction = 0.5;
                } else {
                    temp.blending = false;
                    temp.fraction = atof(mass.c_str());
                }
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
            if(pnl <= 0.92 || pnl > 1)
                cout <<endl << endl << "Warning! Non-leakage value wrong!"<<endl<<" See LEAKG in input file."<<endl<<endl;
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

int main(){
    fuelBundle fuel;
    fuel = InputReader();
    fuel = FuelNormalizer(fuel);
    double flux;
    flux = fluxcalc(fuel);
    DataReader2(fuel.name, fuel.iso);

    vector<nonActinide> nona; //"NONA"ctinide ;)
    nona = NonActinideReader(fuel.name + "/TAPE9.INP");
    fuel = NBuilder(fuel, nona);
    if (fuel.libcheck == true){
        fuel = lib_interpol(fuel);
    }
    if (fuel.operation_type == "BURNUP"){
        pair< double, map < int, double> > test = burnupcalc(fuel, fuel.batch, fuel.pnl, 0.001, flux);
        cout <<"burnup: "<< test.first << endl;
    } else {
        cout << "Enrichment:    " << enrichcalc(fuel.target_BUd, fuel.batch, 0.001, fuel, flux) << endl;
    }
    /*typedef std::map<int, double>::iterator test_map;
    for (test_map iterator = test.second.begin(); iterator != test.second.end(); iterator++){
        cout << iterator->first << "       " << iterator->second << endl;
    }*/

    return 0;
}



