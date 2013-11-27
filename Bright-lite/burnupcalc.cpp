#include "burnupcalc.h"

using namespace std;

double intpol(double y0, double y1, double x0, double x1, double x)
{
    // linear interpolation function
    double y = y0 + (y1 - y0)*(x - x0)/(x1 - x0);
    return y;
}



// need to add units

map<int, double> tomass (int ti, double time, isoInformation isoinfo) {
    map<int, double> out = map<int, double>();
    double mass_i;
    string name_i;
    int nucid;
    for (int i = 0; i < isoinfo.iso_vector.size(); i++){
        name_i = isoinfo.iso_vector[i].name;
        nucid = pyne::nucname::zzaaam(name_i)/10;
        mass_i = intpol(isoinfo.iso_vector[i].mass[ti-1],
                        isoinfo.iso_vector[i].mass[ti],
                        isoinfo.time[ti-1],
                        isoinfo.time[ti],
                        time);
        out[nucid] = mass_i;
    }
    return out;
}


pair<double, map<int, double> > burnupcalc(isoInformation tempone, int N, double tolerance) {
    pair<double, map<int,double> > rtn(0, map<int, double>());
    double mass = 1; //in kilograms
    double BU_f, time_f; // burnup and time when k reaches one, time in days
    double BU_total = 0;
    double BU_finder = 0;
    double BU_n = 0;    // estimated burnup of the n th batch
    double time_finder = 0;
    double k_total =0;
    int i =0;
    double k_batch [N];
    double t_batch;
    double x0=0;
    double x1=0;
    int m =0;
    while (i< tempone.neutron_prod.size())
    {
        tempone.k_inf.push_back(tempone.neutron_prod[i]/tempone.neutron_dest[i]);
        i++;
    }

    i=0;
    while (tempone.k_inf[i] > 1)
        i++;            // finds the number of entry when k drops under 1


    // below line interpolates the end point for burnup when k is one
    BU_f = intpol(tempone.BUd[i-1],tempone.BUd[i],tempone.k_inf[i-1],tempone.k_inf[i],1.0);

    time_f = intpol(tempone.time[i-1],tempone.time[i],tempone.k_inf[i-1],tempone.k_inf[i],1.0);


    int j=0;
    while (j < i)
    {
        BU_total = tempone.BUd[j] + BU_total;
        j++;
    }

    BU_total = BU_total + BU_f; // adds the last value when k reaches one


    if (N == 1){
        rtn.first = BU_total;
        rtn.second = tomass(i, time_f, tempone);
        return rtn;
    }
    BU_total = 2*N*BU_total/(N+1); //linear approximation of mutli batch burnup, used to find a good initial guess of the max burnup


    while (1) {
        for (j = 0; j<N; j++) {
            BU_finder = 0;
            BU_n = BU_total*(j+1)/N;

            i=0;
            while (BU_finder + tempone.BUd[i] < BU_n)  // finds the discrete point i corresponding to the burnup
                {
                    BU_finder = BU_finder + tempone.BUd[i];
                    i++;
                }

            while (m < i) {
                x0 = x0 + tempone.BUd[m];  // sums the burnup for total burnup
                m++;
            }
            x1 = x0 + tempone.BUd[m];     // adds on more discrete point for linear interpolation

            k_batch[j] = intpol(tempone.k_inf[i], tempone.k_inf[i+1], x0, x1, BU_n);  //finds the k of the batch

            x0 = 0;
            x1= 0;
            m = 0;
        }

        j = 0;
        while (j < N) { //sums the k values of every batch
            k_total = k_total + k_batch[j];
            j++;
        }
        k_total = k_total/N;

        cout << k_total << endl;

        if (abs(1 - k_total) < 0.0001 ) //breaks out of loop if k is close enough, tolerance value passed to the function can be used here
        break;


        BU_total = intpol(0,BU_total,tempone.k_inf[0],k_total,1); // updates the guess using (k(0),0) and (k_total, BU_total)



        k_total = 0;

    }

        rtn.first = BU_total;
        rtn.second = tomass(i, time_f, tempone);
        return rtn;

}



double enrichcalc(double BU_end, int N, double tolerance, int type, vector<isoInformation> input_stream)
{

double X;
double BU_guess;
double BU2, BU7;
isoInformation test2;


// accurate guess calc, extrapolating from two data points at 2 and 7% enrichment
BU2 = 20;
BU7 = 100;

X = 0.02 + (BU_end - BU2)*(0.07 - 0.02)/(BU7 - BU2);

BU_guess = burnupcalc(DataReader(test2, type, input_stream), N, 0.01).first;

// enrichment iteration
while (BU_end < BU_guess)
{
    X = X - 0.001;
    BU_guess = burnupcalc(DataReader(test2, type, input_stream), N, 0.01 ).first;
}


while (BU_end > BU_guess)
{
    X = X + 0.001;
    BU_guess = burnupcalc(DataReader(test2, type, input_stream), N, 0.01 ).first;
}
    return X;


}

int main(){
    isoInformation testVector;
    double BUd_sum = 0;
    int N;
    double X;
    isoInformation test1;
    vector<isoInformation> input_stream;
    ifstream inf("../inputFile.txt");
    string line;
    double mass_total;
    while (getline(inf, line)) {
        isoInformation temp_iso;
        istringstream iss(line);
        iss >> temp_iso.name;
        iss >> temp_iso.fraction;
        mass_total = mass_total + temp_iso.fraction;
        input_stream.push_back(temp_iso);
    }
    for (int i = 0; i < input_stream.size(); i++){
        input_stream[i].fraction = input_stream[i].fraction / mass_total;
    }
    int dips;
    cout << "1. LWR" << endl << "2. DUPIC" << endl;
    cin >> dips;

    double BU_end;
    int ip;
    cout << "1. Enrichment to Burnup" << endl << "2. Burnup to Enrichment" << endl;
    cin >> ip;

    switch (ip)
    {
    case 1:
        cout << "Enter number of batches: ";
        cin >> N;
        cout << "Burnup is  " << burnupcalc(DataReader(test1, dips, input_stream), N, .01).first << endl << endl ;
        break;
    case 2:
        cout << "Enter desired Burnup (0-200): ";
        cin >> BU_end;
        cout << "Enter number of batches: ";
        cin >> N;
        cout << "Desired enrichment is about:  " << enrichcalc(BU_end, N, 1, dips, input_stream)*100 << " %" << endl << endl;
        break;
    default:
        cout << endl<< "yeaah, no" << endl << endl;

    }
  return 0;
}
























