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


isoInformation FuelBuilder(vector<isoInformation> fuel_values){
    isoInformation fuel;
    for(int mm = 0; mm < fuel_values.size(); mm++){
        if(fuel_values[mm].iso_vector.size() > 0){
            if(fuel.fluence.size() < 1){
                for (int i = 0; i < fuel_values[mm].fluence.size(); i++){
                    fuel.fluence.push_back(fuel_values[mm].fluence[i]);
                }
                for (int i = 0; i < fuel_values[mm].neutron_prod.size(); i++){
                    fuel.neutron_prod.push_back(fuel_values[mm].fraction[0]*fuel_values[mm].neutron_prod[i]);
                }
                for (int i = 0; i < fuel_values[mm].neutron_dest.size(); i++){
                    fuel.neutron_dest.push_back(fuel_values[mm].fraction[0]*fuel_values[mm].neutron_dest[i]);

                }
                for (int i = 0; i < fuel_values[mm].BUd.size(); i++){
                    fuel.BUd.push_back(fuel_values[mm].fraction[0]*fuel_values[mm].BUd[i]);
                }
                for (int i = 0; i < fuel_values[mm].iso_vector.size(); i++){
                    fuel.iso_vector.push_back(fuel_values[mm].iso_vector[i]);
                    for(int k = 0; k < fuel.iso_vector[i].mass.size(); k++){
                        fuel.iso_vector[i].mass[k] = fuel_values[mm].fraction[0]*fuel.iso_vector[i].mass[k];
                    }
                }
            }else{
                for (int i = 0; i < fuel_values[mm].neutron_prod.size(); i++){
                    fuel.neutron_prod[i] += fuel_values[mm].fraction[0]*fuel_values[mm].neutron_prod[i];
                }
                for (int i = 0; i < fuel_values[mm].neutron_dest.size(); i++){
                    fuel.neutron_dest[i] += fuel_values[mm].fraction[0]*fuel_values[mm].neutron_dest[i];
                }
                for (int i = 0; i < fuel_values[mm].BUd.size(); i++){
                    fuel.BUd[i] += fuel_values[mm].fraction[0]*fuel_values[mm].BUd[i];
                }
                for (int i = 0; i < fuel_values[mm].iso_vector.size(); i++){
                    bool iso_check = true;
                    for(int j = 0; j < fuel.iso_vector.size(); j++){
                        if (fuel_values[mm].iso_vector[i].name == fuel.iso_vector[j].name){
                            for(int k = 0; k < fuel.iso_vector[j].mass.size(); k++){
                                for(int ii = 0; ii < fuel_values[mm].iso_vector[i].mass.size(); ii ++){
                                    if ( k ==ii ){
                                        fuel.iso_vector[j].mass[k] += fuel_values[mm].fraction[0]*fuel_values[mm].iso_vector[i].mass[ii];
                                    }
                                }
                            }
                            iso_check = false;
                        }
                    }
                    if (iso_check == true) {
                        fuel.iso_vector.push_back(fuel_values[mm].iso_vector[i]);
                        for(int k = 0; k < fuel.iso_vector[fuel.iso_vector.size()-1].mass.size()-1; k++){
                            fuel.iso_vector[fuel.iso_vector.size()-1].mass[k] = fuel.iso_vector[fuel.iso_vector.size()-1].mass[k]*fuel_values[mm].fraction[0];
                        }
                    }
                }
            }
        }
        else{
            if(fuel.neutron_prod.size()>1){
                for (int i = 0; i < fuel_values[mm].neutron_prod.size(); i++){
                    fuel.neutron_prod[i] = fuel.neutron_prod[i] + fuel_values[mm].fraction[0]*fuel_values[mm].neutron_prod[i];
                }
                for (int i = 0; i < fuel_values[mm].neutron_dest.size(); i++){
                    fuel.neutron_dest[i] = fuel.neutron_dest[i] + fuel_values[mm].fraction[0]*fuel_values[mm].neutron_dest[i];
                }
            }else{
                for (int i = 0; i < fuel_values[mm].neutron_prod.size(); i++){
                    fuel.neutron_prod.push_back(fuel_values[mm].fraction[0]*fuel_values[mm].neutron_prod[i]+(fuel_values[mm].fraction[0])*fuel_values[1].neutron_prod[i]);
                }
                for (int i = 0; i < fuel_values[mm].neutron_dest.size(); i++){
                    fuel.neutron_dest.push_back(fuel_values[mm].fraction[0]*fuel_values[mm].neutron_dest[i]+(fuel_values[mm].fraction[0])*fuel_values[1].neutron_dest[i]);
                }
            }
        }
    }
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

isoInformation BuildIsotope2(ifstream &input, isoInformation &iso){
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
                    iso.fluence.push_back(value);
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

vector<isoInformation> DataReader2(string type, vector<isoInformation> &input_stream){
    for (int i = 0; i < input_stream.size(); i++){
        if(input_stream[i].type == *"A"){       //cem added this
            ifstream inf(type + "/" +to_string(input_stream[i].name) + ".txt");
            if(!inf){
                cout << "Failed to read file for " + type + " " +  to_string(input_stream[i].name) << endl;
            }
            BuildIsotope2(inf, input_stream[i]);
            inf.close();
        }
    }
    return input_stream;
}

