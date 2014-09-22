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
        
        //add the structural material effects
        for(int j = 0; j < fuel.batch[i].collapsed_iso.neutron_prod.size(); j++){
            //IF theres no rate from actinide, theres no rate from struct
            if(fuel.batch[i].collapsed_iso.neutron_prod[j] != 0){
                fuel.batch[i].collapsed_iso.neutron_prod[j] += fuel.struct_prod;
            }
            if(fuel.batch[i].collapsed_iso.neutron_dest[j] != 0){
                fuel.batch[i].collapsed_iso.neutron_dest[j] += fuel.struct_dest;
            }
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
// need to add units

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
        cout << "ii: " << ii << endl;
        if(ii == 0){
            core.batch[i].rflux = 1/core.batch[i].collapsed_iso.neutron_prod[0];
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

double prod_to_nusigf(double prod){
    double NuSig_f;
    
    NuSig_f = prod*19.1;
    
    NuSig_f = 0.0203;
    
    return NuSig_f;
}

double dest_to_siga(double dest){
    double Sig_a;
    
    Sig_a = dest*19.1;
    
    Sig_a = 0.0200;
    
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
        
        prod = intpol(core.batch[i].collapsed_iso.neutron_prod[ii-1], core.batch[i].collapsed_iso.neutron_prod[ii], core.batch[i].collapsed_iso.fluence[ii-1], core.batch[i].collapsed_iso.fluence[ii], core.batch[i].Fg);
        NuSigma_f[i] = prod_to_nusigf(prod);
        
        dest = intpol(core.batch[i].collapsed_iso.neutron_dest[ii-1], core.batch[i].collapsed_iso.neutron_dest[ii], core.batch[i].collapsed_iso.fluence[ii-1], core.batch[i].collapsed_iso.fluence[ii], core.batch[i].Fg);
        Sigma_a[i] = dest_to_siga(dest);    
        
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
        for(j = 0; core.batch[i].collapsed_iso.fluence[j] < core.batch[i].Fg; j++){}

        prod_tot += intpol(core.batch[i].collapsed_iso.neutron_prod[j-1], core.batch[i].collapsed_iso.neutron_prod[j], core.batch[i].collapsed_iso.fluence[j-1], core.batch[i].collapsed_iso.fluence[j], core.batch[i].Fg);
      
        dest_tot += intpol(core.batch[i].collapsed_iso.neutron_dest[j-1], core.batch[i].collapsed_iso.neutron_dest[j], core.batch[i].collapsed_iso.fluence[j-1], core.batch[i].collapsed_iso.fluence[j], core.batch[i].Fg);   
             
    }
    if(pnl < 0 || pnl > 1){
        cout << endl << "Error in core nonleakage! Assumed new value: 1" << endl;
        pnl = 1;
    }
    
    return prod_tot * pnl / dest_tot;    
}

fuelBundle burnupcalc(fuelBundle core, int mode, double delta) {
    //this function only uses the COLLAPSED_ISO of each BATCH in the structure CORE
    //all factors that contribute to a change in neutron prod/dest rates have to be factored
    //      before calling this function
    cout << endl << "Burnupcalc" << endl;
    
    double BUg = 40; //burnup guess, gets updated if a guess in passed
    int N = core.batch.size(); //number of batches
    double dt = delta*24*60*60; //ten days [s]
    double kcore, kcore_prev;
    double y0, y1, x0, x1;
    double burnup = 0;
    
    //assign the batch_fluence to Fg
    for(int i = 0; i < N; i++){
        core.batch[i].Fg = core.batch[i].batch_fluence;
        cout << "batch flunce in burnupcalc: " << core.batch[i].batch_fluence <<endl;
    }  
      
    //turn zero fluences to one to avoid zero rates
    for(int i = 0; i < N; i++){
        if(core.batch[i].Fg == 0){core.batch[i].Fg = 1;}
    }    
    
    
    
    kcore = kcalc(core);
    if(kcore < 1){
        cout << endl << "Error! Core is not critical." << endl;
    }
    
    //more forward in time until kcore drops under 1
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
        }else{
            cout << endl << "Error in mode input for batch-level flux calculation." << endl;
            return core;
        }
        //update fluences
        for(int i = 0; i < N; i++){
            core.batch[i].Fg += core.batch[i].rflux * core.base_flux * dt;
            cout << "    batch: " << i+1 << "  rflux:" << core.batch[i].rflux << "  Fg:" << core.batch[i].Fg << endl;
        }
        
        kcore = kcalc(core);
        cout << "k: " << kcalc(core) << endl;
    }
    
    //update core fluences

    for(int i = 0; i < N; i++){
        //y0 is the fluence value before the last interation
        y0 = core.batch[i].Fg - core.batch[i].rflux * core.base_flux * dt;
        y1 = core.batch[i].Fg;
        core.batch[i].batch_fluence = intpol(y0, y1, kcore_prev, kcore, 1);
        cout << "  fluence end of burnupcalc: " << core.batch[i].batch_fluence << endl;
        
    }    
    
    //update current composition of batches
    for(int i = 0; i < N; i++){
        core.batch[i].comp.clear();
        int ii;
        for(ii = 0; core.batch[i].collapsed_iso.fluence[ii] < core.batch[i].batch_fluence; ii++){}
        
        for(int j = 0; j < core.batch[i].collapsed_iso.iso_vector.size(); j++){
            core.batch[i].comp[core.batch[i].collapsed_iso.iso_vector[j].name] = 
            intpol(core.batch[i].collapsed_iso.iso_vector[j].mass[ii-1],core.batch[i].collapsed_iso.iso_vector[j].mass[ii], core.batch[i].collapsed_iso.fluence[ii-1], core.batch[i].collapsed_iso.fluence[ii], core.batch[i].batch_fluence)/1000;
        }
    }
    
    //the oldest batch is always index=0
    int ii;
    for(ii = 0; core.batch[0].collapsed_iso.fluence[ii] < core.batch[0].batch_fluence; ii++){}
    burnup = intpol(core.batch[0].collapsed_iso.BU[ii-1], core.batch[0].collapsed_iso.BU[ii], core.batch[0].collapsed_iso.fluence[ii-1], core.batch[0].collapsed_iso.fluence[ii], core.batch[0].batch_fluence);
    
    cout << "Discharge burnup: " << burnup << endl;
        
    /************************output file*********************************/
    std::ofstream outfile;
    outfile.open("../output_cyclus_recent.txt", std::ios::app);
    
    outfile << "Discharge burnup: " << burnup; 
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

/*
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
            for(int i = 0; i < fuel.all_iso.size(); i++){
                if(fuel.all_iso[i].name == nucid){
                    fuel.all_iso[i].sigs = sigs;
                    fuel.all_iso[i].siga = siga;
                    fuel.all_iso[i].fuel = true;
                }
            }
        }
	}


    return fuel;
}
*/






/*
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


////moderator
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
*/
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
/*
        F = x * boost::math::cyl_bessel_i(0,x) / (2 * boost::math::cyl_bessel_i(1, x));
        E = (z*z - y*y) / (2 * y) * ( (boost::math::cyl_bessel_i(0, y) * boost::math::cyl_bessel_k(1, z)+ boost::math::cyl_bessel_k(0, y) *
                                       boost::math::cyl_bessel_i(1, z)) / (boost::math::cyl_bessel_i(1, z) *
                                        boost::math::cyl_bessel_k(1, y) - boost::math::cyl_bessel_k(1, z) * boost::math::cyl_bessel_i(1, y)));
        f = pow((((Sig_aM * V_M)/(Sig_aF * V_F)) * F + E), (-1.));
        cout << "Disadvtg: " << (Sig_aF*V_F - f*Sig_aF*V_F)/(f*Sig_aM*V_M)<<endl;

    }

    return (Sig_aF*V_F - f*Sig_aF*V_F)/(f*Sig_aM*V_M);
}
*/

/*
pair<double, pair<double, map<int,double> > > enrichcalc(double BU_end, int N, double tolerance, fuelBundle fuel) {
    pair<double, pair<double, map<int, double> > > rtn;
    double X;
    double BU_guess, enrich_guess;
    double enrich_previous, BU_previous;
    double enrich_lower, enrich_upper;
    double BU_lower, BU_upper;
    vector <int> blend_vector;

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
*/

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
/*
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
