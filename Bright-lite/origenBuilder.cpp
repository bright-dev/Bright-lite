#include<iostream>
#include<vector>
#include<regex>
#include<iterator>
#include<fstream>
#include<algorithm>
#include "structures.h"

using namespace std;

isoInformation BuildIsotope(ifstream &input){
    isoInformation isotope;
    int i = 0;
    string buffer;
    double value;
    string line;
    while(getline(input, line)){
        istringstream iss(line);
        iss >> buffer;
        if (i >= 4){
            daughter daughter;
            daughter.name = buffer;
            while (iss >> value){
                daughter.mass.push_back(value);
            }
            isotope.iso_vector.push_back(daughter);
        } else {
            while (iss >> value){
                if (i == 0){
                    isotope.time.push_back(value);
                } else if (i == 1){
                    isotope.neutron_prod.push_back(value);
                } else if (i == 2){
                    isotope.neutron_dest.push_back(value);
                } else if (i == 3){
                    isotope.BUd.push_back(value);
                }
            }
            i++;
        }
    }
    return (isotope);
}

isoInformation FuelBuilder(vector<isoInformation> fuel_values, double u235_mass){
    isoInformation fuel;
    for (int i = 0; i < fuel_values[0].time.size(); i++){
        fuel.time.push_back(u235_mass*fuel_values[0].time[i]+(1.0-u235_mass)*fuel_values[1].time[i]);
    }
    for (int i = 0; i < fuel_values[0].neutron_prod.size(); i++){
        fuel.neutron_prod.push_back(u235_mass*fuel_values[0].neutron_prod[i]+(1.0-u235_mass)*fuel_values[1].neutron_prod[i]);
    }
    for (int i = 0; i < fuel_values[0].neutron_dest.size(); i++){
        fuel.neutron_dest.push_back(u235_mass*fuel_values[0].neutron_dest[i]+(1.0-u235_mass)*fuel_values[1].neutron_dest[i]);
    }
    for (int i = 0; i < fuel_values[0].BUd.size(); i++){
        fuel.BUd.push_back(u235_mass*fuel_values[0].BUd[i]+(1.0-u235_mass)*fuel_values[1].BUd[i]);
    }
    for (int i = 0; i < fuel_values[0].iso_vector.size(); i++){
        fuel.iso_vector.push_back(fuel_values[0].iso_vector[i]);
        for(int k = 0; k < fuel.iso_vector[i].mass.size(); k++){
            fuel.iso_vector[i].mass[k] = u235_mass*fuel.iso_vector[i].mass[k];
        }
    }
    for (int jj = 1; jj < fuel_values.size(); jj ++){
        for (int i = 0; i < fuel_values[jj].iso_vector.size(); i++){
            bool iso_check = true;
            for(int j = 0; j < fuel.iso_vector.size(); j++){
                if (fuel_values[jj].iso_vector[i].name == fuel.iso_vector[j].name){
                    for(int k = 0; k < fuel.iso_vector[j].mass.size(); k++){
                        for(int ii = 0; ii < fuel_values[jj].iso_vector[i].mass.size(); ii ++){
                            if ( k ==ii ){
                                fuel.iso_vector[j].mass[k] = fuel.iso_vector[j].mass[k] + (1.0-u235_mass)*fuel_values[jj].iso_vector[i].mass[ii];
                            }
                        }
                    }
                }
                iso_check = false;
            }
            if (iso_check == true) {
                fuel.iso_vector.push_back(fuel_values[jj].iso_vector[i]);
            }
        }
    }
    return fuel;
}

isoInformation DataReader(isoInformation test1, double X){
    vector<isoInformation> mass_stream;
    isoInformation iso_info;
    isoInformation iso_info1;
    ifstream inf("/home/robert/Bright-lite/u235data.txt");
    ifstream inf1("/home/robert/Bright-lite/u238data.txt");
    ofstream outf("/home/robert/Bright-lite/test.txt");
    if (!inf){
        cerr << "Could not read file yo for U-235\n";
    }
    if (!inf1){
        cerr << "Could not readdf file for U-238\n";
    }
    iso_info = BuildIsotope(inf);
    mass_stream.push_back(iso_info);
    iso_info1 = BuildIsotope(inf1);
    mass_stream.push_back(iso_info1);

    inf.close();
    inf1.close();
    test1 = FuelBuilder(mass_stream, X);
    outf << "TIME" << "    ";
    for(int i = 0; i < test1.time.size(); i++){
        outf << test1.time[i] << "   ";
    }
    outf << "\n" << "NEUT_PROD" << "    ";
    for(int i = 0; i < test1.neutron_prod.size(); i++){
        outf << test1.neutron_prod[i] << "   ";
    }
    outf << "\n" << "NEUT_DEST" << "    ";
    for(int i = 0; i < test1.neutron_dest.size(); i++){
        outf << test1.neutron_dest[i] << "   ";
    }
    outf << "\n" << "BUd" << "    ";
    for(int i = 0; i < test1.BUd.size(); i++){
        outf << test1.BUd[i] << "    ";
    }
    outf << "\n" << "";
    for (int i = 0; i < test1.iso_vector.size(); i++){
        if (test1.iso_vector[i].mass[11] > 0.01){
            outf << test1.iso_vector[i].name << "    ";
            for (int j = 0; j < test1.iso_vector[i].mass.size(); j++){
                outf << test1.iso_vector[i].mass[j] << "    ";
            }
            outf << "\n";
        }
    }
    outf.close();
    return test1;
}
