#include<iostream>
#include<vector>
#include<fstream>
#include<sstream>
#include "origenBuilder.h"

using namespace std;

vector< vector<double> > ParseOriginFile(string file_location){
    string tape = "/TAPE6";
    ifstream inf("C:/Users/Robert/Documents/TAPE6.txt");
    // Read error message.
    if (!inf){
        cerr << "Could not read file " + tape + " in " + file_location;
    }
    vector< vector<double> > iso_info;
    string line;
    vector<double> neutron_prod;
    vector<double> neutron_dest;
    vector<double> k_inf;
    vector<double> BUd;
    while (getline(inf, line)) {
        if (line.find("NEUT PRODN") == 0){
            istringstream iss(line);
            string test;
            iss >> test >> test;
            double value;
            while(iss >> value){
                neutron_prod.push_back(value);
            }
        }
        if (line.find("NEUT DESTN") == 0){
            istringstream iss(line);
            string test;
            iss >> test >> test;
            double value;
            while (iss >> value){
                neutron_dest.push_back(value);
            }
        }
        if (line.find("K INFINITY") == 0){
            istringstream iss(line);
            string test;
            iss >> test >> test;
            double value;
            while (iss >> value){
                k_inf.push_back(value);
            }
        }
        if (line.find("BURNUP,MWD") == 0){
            istringstream iss(line);
            string test;
            iss >> test;
            double value;
            while (iss >> value){
                BUd.push_back(value);
            }
        }
    }
    inf.close();
    iso_info.push_back(neutron_prod);
    iso_info.push_back(neutron_dest);
    iso_info.push_back(k_inf);
    iso_info.push_back(BUd);

    return iso_info;
}

int main(){
    vector< vector<double> > testVector;
    testVector = ParseOriginFile("C:/Users/Robert/Documents");
    for (int i = 1; i < testVector[0].size(); i ++){
        cout << ((testVector[0][i]/testVector[1][i]) - testVector[2][i])/testVector[2][i]*100;
        cout << "     " << testVector[3][i];
        cout << "\n";
    }
    OrigenTemplateBuilder();
    return 0;
}
