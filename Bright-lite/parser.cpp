#include<iostream>
#include<fstream>

vector ParseOriginFile(String file_location){
    ifstream inf(file_location + "TAPE9.OUT");
    // Read error message.
    if (!inf){
        cerr << "Could not read file TAPE9.OUT in " + file_location;
        exit(1);
    }

    while (inf) {
        std::string strInput;

    }

}

int FindNeutronProduction (String line) {

}
