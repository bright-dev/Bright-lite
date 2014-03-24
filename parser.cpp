#include "parser.h"

using namespace std;

double IssToDouble(istringstream &iss){
    double num1;
    iss >> num1;
    return (num1);
}

isoInformation ParseOriginFile(string file_location){
    string tape = "/TAPE6";
    ifstream inf("home/Robert/Bright-lite/U-235/TAPE6.out");
    ofstream outf;
    outf.open("C:/Users/Robert/Documents/TAPE61.out");
    // Read error message.
    if (!inf){
        cerr << "Could not read file " + tape + " in " + file_location;
    }
    if (!outf){
        cerr << "Could not write file\n";
    }
    isoInformation iso_info;
    string line;
    while (getline(inf, line)) {
        if (line.find("NEUT PRODN") == 0){
            istringstream iss(line);
            string test;
            iss >> test >> test;
            double value;
            while(iss >> value){
                iso_info.neutron_prod.push_back(value);
            }
        }
        if (line.find("NEUT DESTN") == 0){
            istringstream iss(line);
            string test;
            iss >> test >> test;
            double value;
            while (iss >> value){
                iso_info.neutron_dest.push_back(value);
            }
        }
        if (line.find("K INFINITY") == 0){
            istringstream iss(line);
            string test;
            iss >> test >> test;
            double value;
            while (iss >> value){
                iso_info.k_inf.push_back(value);
            }
        }
        if (line.find("BURNUP,MWD") == 0){
            istringstream iss(line);
            string test;
            iss >> test;
            double value;
            while (iss >> value){
                iso_info.BUd.push_back(value);
            }
        }
        if(line.find("TIME, SEC") == 45){
            istringstream iss(line);
            double time_pass;
            string time;
            string line2;
            double flux;
            getline(inf, line2);
            istringstream iss_flux(line2);
            iss >> time >> time >> flux;
            while(iss >> time_pass){
                iso_info.fluence.push_back(time_pass*flux);
            }
        }
        if(line.find("5 SUMMARY TABLE:  C") == 21){
            getline(inf, line);
            getline(inf, line);
            double num1;
            double num2;
            string expon;
            int j = 0;
            bool test = false;
            while( getline(inf, line)){
                if(line.find("SUMTOT") == 0 || j == 51) {
                    break;
                } else {
                    j += 1;
                    string iso_name = line.substr(0, 11);
                    iso_name.erase(remove(iso_name.begin(), iso_name.end(), ' '), iso_name.end());
                    if(iso_name.find_first_of("0123456789") != std::string::npos){
                        for (int i = 0; i < iso_info.iso_vector.size(); i++){
                            if (iso_info.iso_vector[i].name == iso_name){
                                test = true;
                                for ( int j = 11; j < 130 ; j += 10){
                                    istringstream iss(line.substr(j,9));
                                    iso_info.iso_vector[i].mass.push_back(IssToDouble(iss));
                                }
                            }
                        }
                        if (test == false){
                            daughter daughter;
                            daughter.name = iso_name;
                            for ( int j = 11; j < 130 ; j += 10){
                                istringstream iss(line.substr(j,9));
                                daughter.mass.push_back(IssToDouble(iss));
                            }
                            iso_info.iso_vector.push_back(daughter);
                        }
                    }
                }
            }
        }
    }
    inf.close();
    outf << "TIME" << "    ";
    for(int i = 0; i < iso_info.fluence.size(); i++){
        outf << iso_info.fluence[i] << "   ";
    }
    outf << "\n" << "NEUT_PROD" << "    ";
    for(int i = 0; i < iso_info.neutron_prod.size(); i++){
        outf << iso_info.neutron_prod[i] << "   ";
    }
    outf << "\n" << "NEUT_DEST" << "    ";
    for(int i = 0; i < iso_info.neutron_dest.size(); i++){
        outf << iso_info.neutron_dest[i] << "   ";
    }
    outf << "\n" << "BUd" << "    ";
    for(int i = 0; i < iso_info.BUd.size(); i++){
        outf << iso_info.BUd[i] << "    ";
    }
    outf << "\n" << "";
    for (int i = 0; i < iso_info.iso_vector.size(); i++){
        if (iso_info.iso_vector[i].mass[11] > 0.01){
            outf << iso_info.iso_vector[i].name << "    ";
            for (int j = 0; j < iso_info.iso_vector[i].mass.size(); j++){

                outf << iso_info.iso_vector[i].mass[j] << "    ";
            }
            outf << "\n";
        }
    }
    outf.close();
    return iso_info;
}

















