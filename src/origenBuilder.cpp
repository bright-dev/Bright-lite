#include<iostream>
#include<vector>
#include<regex>
#include<iterator>
#include<fstream>
#include<algorithm>
#include<iostream>
#include<vector>
#include<regex>
#include<iterator>
#include<fstream>
#include<algorithm>
#include<map>
#include "structures.h"
#include "origenBuilder.h"
#include "cyclus.h"
#include "structures.h"

using namespace std;

fuelBundle StructReader(fuelBundle &core){
//reads in structural material fractions from input file and adds the total
//values to core.struct_prod and core.struct_dest
//only needs to be called once to populate core.struct_prod and dest
    //cout << "Struct reader called \n";
    int nucid;
    double fraction;
    double tot_dest = 0;
    double tot_prod = 0;

    string line;
    ifstream fin(cyclus::Env::GetInstallPath() + "/share/brightlite/" + core.name + "/structural.txt");

    vector<nonActinide> nonas;

    nonas = NonActinideReader(cyclus::Env::GetInstallPath() + "/share/brightlite/" + core.name + "/TAPE9.INP");



	while(getline(fin, line))
	{
        istringstream iss(line);
        iss >> nucid >> fraction;

        for(int i = 0; i < nonas.size(); i++){
            if(nonas[i].name == nucid){
                nucid = nucid % 10000;
                nucid = nucid / 10;

                tot_prod += nonas[i].total_prod * fraction * 1000*0.602/nucid;
                tot_dest += nonas[i].total_dest * fraction * 1000*0.602/nucid;
            }
        }

        //cout << nucid << "  " << fraction << endl;

    }
    //cout << "Structs: " << tot_prod << " " << tot_dest << endl;

    core.struct_prod = tot_prod;
    core.struct_dest = tot_dest;


    return core;
}


isoInformation BurnupBuilder(vector<isoInformation> &fuel_values){
    isoInformation fuel;
    //boost::timer t;
    for(int mm = 0; mm < fuel_values.size(); mm++){
        if(fuel_values[mm].iso_vector.size() > 0){
            if(fuel.fluence.size() < 1){
                for (int i = 0; i < fuel_values[mm].fluence.size(); i++){
                    fuel.fluence.push_back(fuel_values[mm].fluence[i]);
                }
                for (int i = 0; i < fuel_values[mm].neutron_prod.size(); i++){
                    fuel.neutron_prod.push_back(fuel_values[mm].fraction*fuel_values[mm].neutron_prod[i]);
                }
                for (int i = 0; i < fuel_values[mm].neutron_dest.size(); i++){
                    fuel.neutron_dest.push_back(fuel_values[mm].fraction*fuel_values[mm].neutron_dest[i]);

                }
                for (int i = 0; i < fuel_values[mm].BUd.size(); i++){
                    fuel.BUd.push_back(fuel_values[mm].fraction*fuel_values[mm].BUd[i]);
                }
            }else{
                for (int i = 0; i < fuel_values[mm].neutron_prod.size(); i++){
                    fuel.neutron_prod[i] += fuel_values[mm].fraction*fuel_values[mm].neutron_prod[i];
                }
                for (int i = 0; i < fuel_values[mm].neutron_dest.size(); i++){
                    fuel.neutron_dest[i] += fuel_values[mm].fraction*fuel_values[mm].neutron_dest[i];
                }
                for (int i = 0; i < fuel_values[mm].BUd.size(); i++){
                    fuel.BUd[i] += fuel_values[mm].fraction*fuel_values[mm].BUd[i];
                }
            }
        }
        else{
            if(fuel.neutron_prod.size()>1){
                for (int i = 0; i < fuel_values[mm].neutron_prod.size(); i++){
                    fuel.neutron_prod[i] = fuel.neutron_prod[i] + fuel_values[mm].fraction*fuel_values[mm].neutron_prod[i];
                }
                for (int i = 0; i < fuel_values[mm].neutron_dest.size(); i++){
                    fuel.neutron_dest[i] = fuel.neutron_dest[i] + fuel_values[mm].fraction*fuel_values[mm].neutron_dest[i];
                }
            }else{
                for (int i = 0; i < fuel_values[mm].neutron_prod.size(); i++){
                    fuel.neutron_prod.push_back(fuel_values[mm].fraction*fuel_values[mm].neutron_prod[i]+(fuel_values[mm].fraction)*fuel_values[1].neutron_prod[i]);
                }
                for (int i = 0; i < fuel_values[mm].neutron_dest.size(); i++){
                    fuel.neutron_dest.push_back(fuel_values[mm].fraction*fuel_values[mm].neutron_dest[i]+(fuel_values[mm].fraction)*fuel_values[1].neutron_dest[i]);
                }
            }
        }
    }
    //std::cout << "burnupblender " << t.elapsed() << std::endl;
    return fuel;
}

/// make this take the prod and dest as doubles
isoInformation FuelBuilder(vector<isoInformation> &fuel_values){
//takes a vetor of isoinfo and uses .fraction to combine isitopes to create one iso
    isoInformation fuel;
    //boost::timer t;
    for(int mm = 0; mm < fuel_values.size(); mm++){
        if(fuel_values[mm].iso_vector.size() > 0){
            if(fuel.fluence.size() < 1){
                for (int i = 0; i < fuel_values[mm].fluence.size(); i++){
                    fuel.fluence.push_back(fuel_values[mm].fluence[i]);
                }
                for (int i = 0; i < fuel_values[mm].neutron_prod.size(); i++){
                    fuel.neutron_prod.push_back(fuel_values[mm].fraction*fuel_values[mm].neutron_prod[i]);
                }
                for (int i = 0; i < fuel_values[mm].neutron_dest.size(); i++){
                    fuel.neutron_dest.push_back(fuel_values[mm].fraction*fuel_values[mm].neutron_dest[i]);

                }
                for (int i = 0; i < fuel_values[mm].BUd.size(); i++){
                    fuel.BUd.push_back(fuel_values[mm].fraction*fuel_values[mm].BUd[i]);
                }
                for (int i = 0; i < fuel_values[mm].iso_vector.size(); i++){
                    fuel.iso_vector.push_back(fuel_values[mm].iso_vector[i]);
                    for(int k = 0; k < fuel.iso_vector[i].mass.size(); k++){
                        fuel.iso_vector[i].mass[k] = fuel_values[mm].fraction*fuel.iso_vector[i].mass[k];
                    }
                }
            }else{
                for (int i = 0; i < fuel_values[mm].neutron_prod.size(); i++){
                    fuel.neutron_prod[i] += fuel_values[mm].fraction*fuel_values[mm].neutron_prod[i];
                }
                for (int i = 0; i < fuel_values[mm].neutron_dest.size(); i++){
                    fuel.neutron_dest[i] += fuel_values[mm].fraction*fuel_values[mm].neutron_dest[i];
                }
                for (int i = 0; i < fuel_values[mm].BUd.size(); i++){
                    fuel.BUd[i] += fuel_values[mm].fraction*fuel_values[mm].BUd[i];
                }
                for (int i = 0; i < fuel_values[mm].iso_vector.size(); i++){
                    bool iso_check = true;
                    for(int j = 0; j < fuel.iso_vector.size(); j++){
                        if (fuel_values[mm].iso_vector[i].name == fuel.iso_vector[j].name){
                            for(int k = 0; k < fuel.iso_vector[j].mass.size(); k++){
                                for(int ii = 0; ii < fuel_values[mm].iso_vector[i].mass.size(); ii ++){
                                    if ( k ==ii ){
                                        fuel.iso_vector[j].mass[k] += fuel_values[mm].fraction*fuel_values[mm].iso_vector[i].mass[ii];
                                    }
                                }
                            }
                            iso_check = false;
                        }
                    }
                    if (iso_check == true) {
                        fuel.iso_vector.push_back(fuel_values[mm].iso_vector[i]);
                        for(int k = 0; k < fuel.iso_vector[fuel.iso_vector.size()-1].mass.size()-1; k++){
                            fuel.iso_vector[fuel.iso_vector.size()-1].mass[k] = fuel.iso_vector[fuel.iso_vector.size()-1].mass[k]*fuel_values[mm].fraction;
                        }
                    }
                }
            }
        }
        else{
            if(fuel.neutron_prod.size()>1){
                for (int i = 0; i < fuel_values[mm].neutron_prod.size(); i++){
                    fuel.neutron_prod[i] = fuel.neutron_prod[i] + fuel_values[mm].fraction*fuel_values[mm].neutron_prod[i];
                }
                for (int i = 0; i < fuel_values[mm].neutron_dest.size(); i++){
                    fuel.neutron_dest[i] = fuel.neutron_dest[i] + fuel_values[mm].fraction*fuel_values[mm].neutron_dest[i];
                }
            }else{
                for (int i = 0; i < fuel_values[mm].neutron_prod.size(); i++){
                    fuel.neutron_prod.push_back(fuel_values[mm].fraction*fuel_values[mm].neutron_prod[i]+(fuel_values[mm].fraction)*fuel_values[1].neutron_prod[i]);
                }
                for (int i = 0; i < fuel_values[mm].neutron_dest.size(); i++){
                    fuel.neutron_dest.push_back(fuel_values[mm].fraction*fuel_values[mm].neutron_dest[i]+(fuel_values[mm].fraction)*fuel_values[1].neutron_dest[i]);
                }
            }
        }
    }
    //std::cout << "fuelBlending " << t.elapsed() << std::endl;
    return fuel;
}



vector<nonActinide> NonActinideReader(string file_name){
    vector<nonActinide> structural_comps;
    ifstream inf(file_name);
    string e_test = "E";
    string e_space = "*E *";
    if (!inf){
        cout << "Failure to read " << file_name << endl;
    }
    string line, library, iso;
    double sng, sn2n, sna, snp, sngx, sn2nx, yyn;
    while(getline(inf, line)){
        istringstream iss(line);
        iss >> iso >> iso >> iso;
        nonActinide na_iso;
        if (iso == "MATERIAL"){
            while(getline(inf, line)){
                istringstream iss1(line);
                iss1 >> library;
                if (library == "-1"){
                    inf.close();
                    for (int i = 0; i < structural_comps.size(); i++){
                        structural_comps[i].total_prod = 2*structural_comps[i].sn2n + 2*structural_comps[i].sn2nx;
                        structural_comps[i].total_dest = structural_comps[i].snp + structural_comps[i].sng +
                        structural_comps[i].sngx + structural_comps[i].sn2n + structural_comps[i].sn2nx;

                    }
                    return structural_comps;
                }
                iss1 >> iso >> sng >> sn2n >> snp >> sngx >> sn2nx >> yyn;
                na_iso.name = atoi(iso.c_str());
                na_iso.sng = sng;
                na_iso.sn2n = sn2n;
                na_iso.snp = snp;
                na_iso.sngx = sngx;
                na_iso.sn2nx = sn2nx;
                na_iso.yyn = yyn;
                structural_comps.push_back(na_iso);
            }
        }
    }
}


isoInformation BuildIsotope2(ifstream &input, isoInformation &iso, double flux_value){
// builds the iso vector from isotope library database
    int i = 0;
    string buffer;
    double value;
    string line;
    while(getline(input, line)){
        istringstream iss(line);
        iss >> buffer;
        if (i >= 4){
            daughter daughter;
            daughter.name = pyne::nucname::zzaaam(buffer);
            while (iss >> value){
                daughter.mass.push_back(value);
            }
            iso.iso_vector.push_back(daughter);
        } else {
            while (iss >> value){
                if (i == 0){
                    iso.fluence.push_back(value*flux_value*84600);
                } else if (i == 1){
                    iso.neutron_prod.push_back(value);
                } else if (i == 2){
                    iso.neutron_dest.push_back(value);
                } else if (i == 3){
                    iso.BUd.push_back(value);
                }
            }
            i++;
        }
    }
    return (iso);
}

double flux_finder(std::string type){
    //std::cout << "ADESEEE " << type + "/params.txt" << std::endl;
    ifstream inf(type + "/params.txt");
    if(!inf){
        cout << "SADF AD " <<type << "/params.txt" << std::endl;
    }
    string buffer;
    double value;
    string line;
    while(getline(inf, line)){
        istringstream iss(line);
        iss >> buffer >> value;
        if (buffer == "FLUX"){
            return value;
        }
    }
}

vector<isoInformation> DataReader2(string type, vector<isoInformation> &input_stream){
//returns iso for this batch
    //std::cout << type << std::endl;
    type = cyclus::Env::GetInstallPath() + "/share/brightlite/" + type;
    double flux_value = flux_finder(type);
    for (int i = 0; i < input_stream.size(); i++){
        if(true){
            ifstream inf(type + "/" +to_string(input_stream[i].name) + ".txt");
            if(!inf){
                cout << "Failed to read file for " + type + " " +  to_string(input_stream[i].name) << endl;
            }
            BuildIsotope2(inf, input_stream[i], flux_value);
            inf.close();
        }
    }

    return input_stream;
}

