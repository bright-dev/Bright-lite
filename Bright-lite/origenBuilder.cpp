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

void DataReader(){
    vector<isoInformation> mass_stream;
    isoInformation iso_info;
    ifstream inf("C:/Users/Robert/Documents/U-235/u235data.txt");
    ifstream inf1("C:/Users/Robert/Documents/U-238/u238data.txt");
    if (!inf){
        cerr << "Could not read file for U-235\n";
    }
    if (!inf1){
        cerr << "Could not read file for U-238\n";
    }
    iso_info = BuildIsotope(inf);
    mass_stream.push_back(iso_info);
    iso_info = BuildIsotope(inf1);
    mass_stream.push_back(iso_info);

    inf.close();
    inf1.close();
    // This bit is for testing.
    double neutron_prod;
    double neutron_dest;
    double mass;
    for (int i = 0; i < mass_stream[0].neutron_prod.size(); i++){
        neutron_prod = 0;
        neutron_dest = 0;
        for ( int j = 0; j < mass_stream.size(); j++){
            if ( j == 0) {
                mass = 0.03;
            } else {
                mass = 0.97;
            }
            neutron_prod += mass * mass_stream[j].neutron_prod[i];
            neutron_dest += mass * mass_stream[j].neutron_dest[i];
        }
        cout << mass_stream[0].time[i]<< "     " <<neutron_prod / neutron_dest << "\n";
    }
}
