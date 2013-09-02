#include<iostream>
#include<vector>
#include<regex>
#include<iterator>
#include<fstream>
#include<algorithm>

using namespace std;

void OrigenTemplateBuilder(){
    ifstream inf("C:/Users/Robert/Documents/tape5.fluence_stripper");
    ofstream outf("C:/Users/Robert/Desktop/test");
    if (!inf){
        cerr << "Could not read file";
    }
    string line;
    regex e ("/\bFLUX\b");
    string fix = "3E14";
    while(getline(inf, line)){
        regex_replace(line, e, fix);
        outf << line + "\n";
    }
    inf.close();
    outf.close();
}
