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

isoInformation FuelBuilder(vector<isoInformation> fuel_values){
    isoInformation fuel;
    for (int i = 0; i < fuel_values[0].time.size(); i++){
        fuel.time.push_back(fuel_values[0].fraction*fuel_values[0].time[i]+(fuel_values[0].fraction)*fuel_values[1].time[i]);
    }
    for (int i = 0; i < fuel_values[0].neutron_prod.size(); i++){
        fuel.neutron_prod.push_back(fuel_values[0].fraction*fuel_values[0].neutron_prod[i]+(fuel_values[0].fraction)*fuel_values[1].neutron_prod[i]);
    }
    for (int i = 0; i < fuel_values[0].neutron_dest.size(); i++){
        fuel.neutron_dest.push_back(fuel_values[0].fraction*fuel_values[0].neutron_dest[i]+(fuel_values[0].fraction)*fuel_values[1].neutron_dest[i]);
    }
    for (int i = 0; i < fuel_values[0].BUd.size(); i++){
        fuel.BUd.push_back(fuel_values[0].fraction*fuel_values[0].BUd[i]+(fuel_values[0].fraction)*fuel_values[1].BUd[i]);
    }

    for(int j = 1; j < fuel_values.size(); j++){
        for (int i = 0; i < fuel_values[0].neutron_prod.size(); i++){
            fuel.neutron_prod[i] = fuel.neutron_prod[i] + fuel_values[j].fraction*fuel_values[j].neutron_prod[i];
        }
        for (int i = 0; i < fuel_values[0].neutron_dest.size(); i++){
            fuel.neutron_dest[i] = fuel.neutron_dest[i] + fuel_values[j].fraction*fuel_values[j].neutron_dest[i];
        }
        for (int i = 0; i < fuel_values[0].BUd.size(); i++){
            fuel.BUd[i] = fuel.BUd[i] + fuel_values[j].fraction*fuel_values[j].BUd[i];
        }
    }

    for (int i = 0; i < fuel_values[0].iso_vector.size(); i++){
        fuel.iso_vector.push_back(fuel_values[0].iso_vector[i]);
        for(int k = 0; k < fuel.iso_vector[i].mass.size(); k++){
            fuel.iso_vector[i].mass[k] = fuel_values[0].fraction*fuel.iso_vector[i].mass[k];
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
                                fuel.iso_vector[j].mass[k] = fuel.iso_vector[j].mass[k] + fuel_values[jj].fraction*fuel_values[jj].iso_vector[i].mass[ii];
                            }
                        }
                    }
                    iso_check = false;
                }
            }
            if (iso_check == true) {
                fuel.iso_vector.push_back(fuel_values[jj].iso_vector[i]);
            }
        }
    }
    return fuel;
}

isoInformation DataReader(isoInformation test1, int type, vector<isoInformation> input_stream){
    vector<isoInformation> mass_stream;
    for (int i = 0; i < input_stream.size(); i++){
        isoInformation iso_info;
        string dir;
        if (type == 1) {
            dir = "../LWR/";
        } else if (type == 2){
            dir = "../DUPIC/";
        }
        ifstream inf(dir + input_stream[i].name + ".txt");
        if(!inf){
            cout << "NOOOOO" << endl;
        }
        iso_info = BuildIsotope(inf);
        mass_stream.push_back(iso_info);
        mass_stream[i].name = input_stream[i].name;
        mass_stream[i].fraction = input_stream[i].fraction;
        inf.close();
    }
    ofstream outf("../test.txt");
    test1 = FuelBuilder(mass_stream);
    outf << "TIME" << " ";
    for(int i = 0; i < test1.time.size(); i++){
        outf << test1.time[i] << " ";
    }
    outf << "\n" << "NEUT_PROD" << " ";
    for(int i = 0; i < test1.neutron_prod.size(); i++){
        outf << test1.neutron_prod[i] << " ";
    }
    outf << "\n" << "NEUT_DEST" << " ";
    for(int i = 0; i < test1.neutron_dest.size(); i++){
        outf << test1.neutron_dest[i] << " ";
    }
    outf << "\n" << "BUd" << " ";
    for(int i = 0; i < test1.BUd.size(); i++){
        outf << test1.BUd[i] << " ";
    }
    outf << "\n" << "";
    for (int i = 0; i < test1.iso_vector.size(); i++){
        if (test1.iso_vector[i].mass[11] > 0.01){
            outf << test1.iso_vector[i].name << " ";
            for (int j = 0; j < test1.iso_vector[i].mass.size(); j++){
                outf << test1.iso_vector[i].mass[j] << " ";
            }
            outf << "\n";
        }
    }

    outf.close();
    return test1;
}
