#include "burnupcalc.h"
#include <boost/math/special_functions/bessel.hpp>

using namespace std;

double intpol(double y0, double y1, double x0, double x1, double x) {
    // linear interpolation function
    double y = y0 + (y1 - y0)*(x - x0)/(x1 - x0);
    return y;
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


    for (int i = 0; i < region1.size();i++){
        cout << region1[i].name << "    "<<region1[i].fraction << endl;
    }
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
/*
cout<<"Fn: "<< F_n<<" time:"<<tempone.time[i]<<endl;
cout<<"x1: "<<x1<<" BU_n: "<<BU_n<<endl;
*/
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

pair<double, map<int, double> > burnupcalc(isoInformation tempone, int N, double pnl, double tolerance) {
    pair<double, map<int,double> > rtn(0, map<int, double>());
    double time_f; // time when k reaches one, time in days
    double BU_total = 0;
    double k_total = 10;
    int i = 0;
    double BU1, BU2, BU3;

	i = 0;
    while (i < tempone.neutron_prod.size())
    {
        tempone.k_inf.push_back(((tempone.neutron_prod[i]))*pnl/(tempone.neutron_dest[i]));
        i++;
    }
    i=0;
    while (tempone.k_inf[i] > 1.0){
        BU_total += tempone.BUd[i];
        i++; // finds the number of entry when k drops under 1
    }
    BU_total = intpol(BU_total,BU_total+tempone.BUd[i],tempone.k_inf[i-1],tempone.k_inf[i],1.0);
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
    //rtn.second = tomass(i, time_f, tempone);
    return rtn;

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


double enrichcalc(double BU_end, int N, double tolerance, int type, vector<isoInformation> input_stream) {
//THIS WORKS ONLY FOR TWO ISO'S
double X;
double BU_guess;
double BU2, BU6;
isoInformation test2;


// accurate guess calc, extrapolating from two data points at 2 and 7% enrichment
input_stream[0].fraction = 0.02;
input_stream[1].fraction = 0.98;
BU2 = burnupcalc(DataReader(test2, type, input_stream), N, .98, 0.01).first;
input_stream[0].fraction = 0.06;
input_stream[1].fraction = 0.94;
BU6 = burnupcalc(DataReader(test2, type, input_stream), N, .98, 0.01).first;

X = 0.02 + (BU_end - BU2)*(0.06 - 0.02)/(BU6 - BU2);

input_stream[0].fraction = X;
input_stream[1].fraction = 1-X;
BU_guess = burnupcalc(DataReader(test2, type, input_stream), N, .98, 0.01).first;

// enrichment iteration
while (BU_end < BU_guess)
{
    X = X - 0.001;
    input_stream[0].fraction = X;
    input_stream[1].fraction = 1-X;
    BU_guess = burnupcalc(DataReader(test2, type, input_stream), N, .98, 0.01 ).first;
}


while (BU_end > BU_guess)
{
    X = X + 0.001;
    input_stream[0].fraction = X;
    input_stream[1].fraction = 1-X;
    BU_guess = burnupcalc(DataReader(test2, type, input_stream), N, .98, 0.01 ).first;
}
    return X;
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

        abs35 = 608.4-14.95;
        sca35 = 14.95;
        tot35 = 608.4;

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


fuelBundle InputReader(){
    std::string name;
    int region, N, t_res;
    char type, word[8];
    int nucid;
    double mass, pnl;
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
            cout << name << endl;
        }
        if(line.find("REGIONS") == 0){
            while(getline(fin, line)){
                if(line.find("END") == 0) break;
                istringstream iss(line);
                iss >> region >> type >> nucid >> mass;
                temp.name = nucid;
                temp.region = region;
                temp.type = type;
                temp.fraction = mass;
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
	}
    fuel.batch = N;
    fuel.pnl = pnl;
    fuel.tres = t_res;

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

isoInformation lib_interpol(fuelBundle input_fuel, vector<string> libs, vector<interpol_pair> targets){
    vector<fuelBundle> fuel_pairs;
    for (int i = 0; i < libs.size(); i++){
        fuelBundle lib_bundle;
        for(int j = 0; j < input_fuel.iso.size(); j++){
            if(input_fuel.iso[j].type == *"A"){
                isoInformation iso;
                iso.name = input_fuel.iso[j].name;
                iso.fraction = input_fuel.iso[j].fraction;
                iso.type = input_fuel.iso[j].type;
                iso.region = input_fuel.iso[j].region;
                lib_bundle.iso.push_back(iso);
            }
        }
        lib_bundle.name = libs[i];
        lib_bundle.batch = input_fuel.batch;
        lib_bundle.tres = input_fuel.tres;
        lib_bundle.pnl = input_fuel.pnl;
        fuel_pairs.push_back(lib_bundle);
    }
    for (int i = 0; i < fuel_pairs.size(); i++){
        DataReader2(fuel_pairs[i].name, fuel_pairs[i].iso);
    }
    // Reading Metrics //
    vector<vector<double> > metrics;
    std::string metric_name;
    double metric_value;
    for (int i = 0; i < targets.size(); i++){
        vector<double> metric_vect;
        metrics.push_back(metric_vect);
        for (int j = 0; j < fuel_pairs.size(); j++){
            std::string line;
            ifstream inf(fuel_pairs[j].name + "/params.txt");
            while(getline(inf, line)){
                if (line.find(targets[i].metric) == 0){
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
        if (targets[i].value > max || targets[i].value < min){
            alpha = 40;
        }
        for(int j = 0; j < metrics[i].size(); j++){
            metrics[i][j] -= min;
            metrics[i][j] /= max;
        }
        targets[i].scaled_value = targets[i].value - min;
        targets[i].scaled_value /= max;
    }
    // distances
    vector<double> metric_distances;
    for(int i = 0; i < libs.size(); i++){
        double distance_measure = 0;
        for (int j = 0; j < metrics.size(); j++){
            distance_measure += pow(targets[j].scaled_value - metrics[j][i], 2);
        }
        metric_distances.push_back(pow(distance_measure, alpha/2));
    }
    double met_dist_sum;
    for (int i = 0; i < metric_distances.size(); i++){
        met_dist_sum += 1./metric_distances[i];
    }
    // libs
    vector<isoInformation> test_pairs;
    for(int i = 0; i < fuel_pairs.size(); i++){
        double flux = fluxcalc(fuel_pairs[i]);
        vector<nonActinide> nona; //"NONA"ctinide ;)
        nona = NonActinideReader(fuel_pairs[i].name + "/TAPE9.INP");
        fuel_pairs[i] = NBuilder(fuel_pairs[i], nona);
        isoInformation singleiso;
        singleiso = regioncollapse(fuel_pairs[i], flux);
        test_pairs.push_back(singleiso);
    }
    for (int i = 0; i < test_pairs.size(); i++){
        cout <<"burnup: "<< burnupcalc(test_pairs[i], fuel_pairs[i].batch, fuel_pairs[i].pnl, 0.001).first << endl;
    }
    isoInformation test_iso;
    for (int i = 0; i < test_pairs[0].fluence.size(); i++){
        test_iso.fluence.push_back(test_pairs[0].fluence[i]);
    }
    for (int i = 0; i < test_pairs[0].neutron_prod.size(); i++){
        double temp_prod = 0;
        for (int j = 0; j < test_pairs.size(); j++){
            temp_prod += test_pairs[j].neutron_prod[i]/metric_distances[j];
        }
        test_iso.neutron_prod.push_back(temp_prod/met_dist_sum);
    }
    for (int i = 0; i < test_pairs[0].neutron_dest.size(); i++){
        double temp_dest = 0;
        for (int j = 0; j < test_pairs.size(); j++){
            temp_dest += test_pairs[j].neutron_dest[i]/metric_distances[j];
        }
        test_iso.neutron_dest.push_back(temp_dest/met_dist_sum);
    }
    for (int i = 0; i < test_pairs[0].BUd.size(); i++){
        double temp_BUd = 0;
        for (int j = 0; j < test_pairs.size(); j++){
            temp_BUd += test_pairs[j].BUd[i]/metric_distances[j];
        }
        test_iso.BUd.push_back(temp_BUd/met_dist_sum);
    }
    return test_iso;

}

int main(){
    fuelBundle fuel;
    //test

    fuel = InputReader();

    fuel = FuelNormalizer(fuel);
    double flux;
    flux = fluxcalc(fuel);
    cout << flux << endl;

    DataReader2(fuel.name, fuel.iso);

    vector<nonActinide> nona; //"NONA"ctinide ;)
    nona = NonActinideReader(fuel.name + "/TAPE9.INP");
    fuel = NBuilder(fuel, nona);

    isoInformation singleiso;
    singleiso = regioncollapse(fuel, flux);

    /*vector <string> test_libs;

    test_libs.push_back("E7_100");
    test_libs.push_back("E9_100");
    interpol_pair test_pair;
    test_pair.metric = "ENRICHMENT";
    test_pair.value = 8;
    vector<interpol_pair> test_inter;
    test_inter.push_back(test_pair);
    isoInformation singleiso2 = lib_interpol(fuel, test_libs, test_inter);

    cout <<"burnup: "<< burnupcalc(singleiso2, fuel.batch, fuel.pnl, 0.001).first << endl;*/

    cout <<"burnup: "<< burnupcalc(singleiso, fuel.batch, fuel.pnl, 0.001).first << endl;
    /*bool test_check = false;
    double old_burnup = 1;
    double burnup_test = burnupcalc(singleiso2, fuel.batch, fuel.pnl, 0.001).first
    while (!test_check){
        isoInformation singleiso2 = lib_interpol(fuel, test_libs, test_inter)
        old_burnup = burnup_test;
        burnup_test = burnupcalc(singleiso2, fuel.batch, fuel.pnl, 0.001).first
        if (std::abs(burnup_test - old_burnup)/burnup_test < 0.001) test_check = true;
    }
/*
    while(1){
        cin >> BU_end;
        cout<< enrichcalc(BU_end, 100, 0.001, 1, input_stream)<<endl;
    }

cout << burnupcalc(DataReader(test1, 1, input_stream), 1, .01).first << endl;
cout << burnupcalc(DataReader(test1, 1, input_stream), 2, .01).first << endl;
cout << burnupcalc(DataReader(test1, 1, input_stream), 3, .01).first << endl;
cout << burnupcalc(DataReader(test1, 1, input_stream), 5, .01).first << endl;
cout << burnupcalc(DataReader(test1, 1, input_stream), 3, .01).first << endl;
*/

    return 0;
}



