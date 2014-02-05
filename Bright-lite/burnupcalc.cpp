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

/**

burnupcalc takes a reactor-database-vector, batch number,
and tolerance; and finds the achievable burnup of the fuel.
The passed database vector is specific to a reactor and
contains neutron production, destruction, time, burnup,
and isotopic information. burnupcalc initially finds the
k-value of the fuel at time-step, and will return the total
burnup when k=1 in the case of a single batch loading. If
there is more than one batch, burnupcalc initially estimates
the burnup using a linear approximation and checks the
validity of the guess until the core average k-value is
within 0.0001. The function returns this burnup.


Inputs

tempone: This variable is of type isoInformation, and is a
reactor-specific database. It contains neutron production,
neutron destruction, burnup, time, and isotope vectors. The
variable type also contains an empty vector for k-infinity,
which is filled by burnupcalc in accordance with additional
assumptions (such as neutron leakage, which is subtracted
from neutron production).

N: Number of batches used in the core. Used to account for
batch loading and burnup calculations.

tolerance: Currently not utilized, used to adjust precision.


Outputs

Double, total burnup reached by the core.
*/

pair<double, map<int, double> > burnupcalc(isoInformation tempone, int N, double tolerance) {
    pair<double, map<int,double> > rtn(0, map<int, double>());
    double mass = 1; //in kilograms
    double BU_f, time_f; // burnup and time when k reaches one, time in days
    double BU_total = 0;
    double BU_finder = 0;
    double BU_n = 0; // estimated burnup of the n th batch
    double time_finder = 0;
    double k_total = 0;
    int i = 0;
    double k_batch [N];
    double fluxWeight[N], fluxWeight_tot=0;
    double t_batch;
    double x0 = 0;
    double x1 = 0;
    int m = 0;
    double k_old;
    double bud_old = tempone.BUd[0];

    //read the structural material contribution information (production[0],
    //destruction[1], leakage[2]) to s_contr
    double s_contr[3];
    ifstream fin("../LWR/LWRSTRUCT.txt");
    double passer;
    string spasser;
	while(i < 3)
	{
        fin >> spasser >> passer;
        s_contr[i]=passer;
        i++;
	}
	i = 0;
    while (i < tempone.neutron_prod.size())
    {
        tempone.k_inf.push_back((tempone.neutron_prod[i]*s_contr[2]/* + s_contr[0]*/)/(tempone.neutron_dest[i]/* + s_contr[1]*/));
        i++;
    }
    k_old = tempone.k_inf[0];
    i=0;
    while (tempone.k_inf[i] > 1.0){
        //cout << i << "  "<< tempone.k_inf[i] << endl;
        i++; // finds the number of entry when k drops under 1
    }

    // below line interpolates the end point for burnup when k is one
    BU_f = intpol(tempone.BUd[i],tempone.BUd[i+1],tempone.k_inf[i],tempone.k_inf[i+1],1.0);

    time_f = intpol(tempone.time[i],tempone.time[i+1],tempone.k_inf[i],tempone.k_inf[i+1],1.0);


    int j=0;
    while (j < i+1)
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
    BU_total = 2.* N * BU_total/(N+1); //linear approximation of mutli batch burnup, used to find a good initial guess of the max burnup
    int mn = 0;
    while (1) {
        for (j = 0; j < N; j++) {
            BU_finder = 0;
            BU_n = BU_total*(j+1)/N;

            i=0;
            while (BU_finder + tempone.BUd[i] < BU_n) // finds the discrete point i corresponding to the burnup
                {
                    BU_finder = BU_finder + tempone.BUd[i];
                    i++;
                }

            while (m < i) {
                x0 = x0 + tempone.BUd[m]; // sums the burnup for total burnup
                m++;
            }
            x1 = x0 + tempone.BUd[m]; // adds on more discrete point for linear interpolation

            k_batch[j] = intpol(tempone.k_inf[i], tempone.k_inf[i+1], x0, x1, BU_n); //finds the k of the batch
            time_f = intpol(tempone.time[i],tempone.time[i+1],x0, x1, BU_n);
            //flux is inversely proportional to neutron production
            fluxWeight[j] = 1./intpol(tempone.neutron_prod[i],tempone.neutron_prod[i+1], x0, x1, BU_n);
            x0 = 0;
            x1= 0;
            m = 0;
        }

        j = 0;
        while (j < N) { //sums the k values of every batch, weighing for flux, inverse of production
            k_total = k_total + k_batch[j]*fluxWeight[j];
            fluxWeight_tot = fluxWeight_tot + fluxWeight[j];
            j++;
        }

        k_total = k_total/N/fluxWeight_tot;
        cout << k_total << endl;

        if (abs(1 - k_total) < 0.00001 ) //breaks out of loop if k is close enough, tolerance value passed to the function can be used here
            break;

        BU_total = intpol(tempone.BUd[0],BU_total,tempone.k_inf[0],k_total,1);
        //BU_total = intpol(bud_old,BU_total,k_old,k_total,1); // updates the guess using (k(0),0) and (k_total, BU_total)
        k_total = 0;
        fluxWeight_tot = 0;
        mn++;
    }
    rtn.first = BU_total;
    rtn.second = tomass(i, time_f, tempone);
    return rtn;

}


/**
enrichcalc takes the burnup goal, batch number, and tolerance;
and returns the enrichment needed to reach the burnup goal.
The code initially estimates the enrichment needed by
calculating the burnup (for a given number of batches and type
of reactor) of a two and seven percent enriched fuel and then
interpolates/extrapolates from these two points. The code then
calculates the burnup for this estimated enrichment level, and
iterates by changing the guess enrichment until the desired
burnup is achieved by the guess value. This enrichment is then
returned.


Inputs

BU_end: The target burnup. The code finds an enrichment which
will result in a burnup equal to this value.

N: number of batches used in the core.

tolerance: Currently not utilized, used to adjust precision.


Outputs:

Double, the enrichment needed to achieve BU_end.

*/


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
    ifstream inf("inputFile.txt");
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
        cout << input_stream[i].name << " " << input_stream[i].fraction << endl;
    }
    double enrichment = input_stream[0].fraction;
    inf.close();
    int dips;
    /*cout << "1. LWR" << endl << "2. DUPIC" << endl;
    cin >> dips;*/
    map<int, double> test_mass;
    map<int, double>::iterator Iter;
    double BU_end;
    /*int ip;
    cout << "1. Enrichment to Burnup" << endl << "2. Burnup to Enrichment" << endl;
    cin >> ip;*/
    /*switch (ip)
    {
    case 1:
        cout << "Enter number of batches: ";
        cin >> N;
        double BU_d;
        BU_d = burnupcalc(DataReader(test1, dips, input_stream), N, .01).first;
        test_mass = burnupcalc(DataReader(test1, dips, input_stream), N, .01).second;
        cout << "Burnup is " << BU_d << endl;
        for (Iter = test_mass.begin(); Iter != test_mass.end(); ++Iter){
            string m = pyne::nucname::name((*Iter).first);
            if ((*Iter).second > 0.01){
                outf2 << m << " " << (*Iter).second << endl;
                /** STUPID UGLY UGLY CODE*/
                /*if (m == "Am241" || m == "Am243" || m == "Cm242" || m == "Cm244" || m == "Np237" || m == "Np238" || m == "Np239"){
                    outf1 << m << " " << (*Iter).second << endl;
                }
                if (m == "Pu238" || m == "Pu239" || m == "Pu240" || m == "Pu241" || m == "Pu242" || m == "U234" || m == "U235"){
                    outf1 << m << " " << (*Iter).second << endl;
                }
                if (m == "U236" || m == "U237" || m == "U238"){
                    outf1 << m << " " << (*Iter).second << endl;
                }
            }
        }
        outf1.close();
        outf2.close();
        break;
    case 2:
        cout << "Enter desired Burnup (0-200): ";
        cin >> BU_end;
        cout << "Enter number of batches: ";
        cin >> N;
        cout << "Desired enrichment is about: " << enrichcalc(BU_end, N, 1, dips, input_stream)*100 << " %" << endl << endl;
        break;
    default:
        cout << endl<< "yeaah, no" << endl << endl;
    }
    outf1.close();*/
    double BU_d;
    ofstream outf2("outputIsosLWR.txt");
    BU_d = burnupcalc(DataReader(test1, 1, input_stream), 3, .01).first;
    test_mass = burnupcalc(DataReader(test1, 1, input_stream), 3, .01).second;
    cout << "TEST" << endl;
    cout << "Burnup is  " << BU_d << endl;
    /*for (Iter = test_mass.begin(); Iter != test_mass.end(); ++Iter){
        string m = pyne::nucname::name((*Iter).first);
        if ((*Iter).second > 0.01){
            outf2 << m << "   " << (*Iter).second << endl;
            /** STUPID UGLY UGLY CODE
            if (m == "Am241" || m == "Am243" || m == "Cm242" || m == "Cm244" || m == "Np237" || m == "Np238" || m == "Np239"){
                outf1 << m << "  " << (*Iter).second << endl;
            }
            if (m == "Pu238" || m == "Pu239" || m == "Pu240" || m == "Pu241" || m == "Pu242" || m == "U234" || m == "U235"){
                outf1 << m << "  " << (*Iter).second << endl;
            }
            if (m == "U236" || m == "U237" || m == "U238"){
                outf1 << m << "  " << (*Iter).second << endl;
            }
        }
    }
    outf1.close();
    outf2.close();
    ofstream outf3("outputIsosDUPIC.txt");
    ifstream inf1("inputFile.txt");
    string line1;
    isoInformation test2;
    vector<isoInformation> input_stream1;
    double mass_total1;
    map<int, double> test_mass1;
    double BU_d1;
    while (getline(inf1, line1)) {
        isoInformation temp_iso;
        istringstream iss(line1);
        iss >> temp_iso.name;
        iss >> temp_iso.fraction;
        mass_total1 = mass_total1 + temp_iso.fraction;
        input_stream1.push_back(temp_iso);
    }
    for (int i = 0; i < input_stream1.size(); i++){
        input_stream1[i].fraction = input_stream1[i].fraction / mass_total1;
        cout << input_stream1[i].name << "    " << input_stream1[i].fraction << endl;
    }
    inf.close();
    BU_d1 = burnupcalc(DataReader(test2, 2, input_stream1), 300, .01).first;
    test_mass1 = burnupcalc(DataReader(test2, 2, input_stream1), 300, .01).second;
    cout << "Burnup is  " << BU_d1 << endl;
    for (Iter = test_mass1.begin(); Iter != test_mass1.end(); ++Iter){
        string m = pyne::nucname::name((*Iter).first);
        if ((*Iter).second > 0.01){
            outf3 << m << "   " << (*Iter).second << endl;
            /** STUPID UGLY UGLY CODE
            if (m == "Am241" || m == "Am243" || m == "Cm242" || m == "Cm244" || m == "Np237" || m == "Np238" || m == "Np239"){
                outf1 << m << "  " << (*Iter).second << endl;
            }
            if (m == "Pu238" || m == "Pu239" || m == "Pu240" || m == "Pu241" || m == "Pu242" || m == "U234" || m == "U235"){
                outf1 << m << "  " << (*Iter).second << endl;
            }
            if (m == "U236" || m == "U237" || m == "U238"){
                outf1 << m << "  " << (*Iter).second << endl;
            }
        }
    }
    fstream outf4("outBUD", fstream::in | fstream::out | fstream::app);
    outf4 << enrichment << "    " << BU_d << "   " << BU_d1 << endl;
    outf4.close();*/
    return 0;
}




