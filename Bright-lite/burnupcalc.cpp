#include<iostream>
#include<vector>
#include<regex>
#include<iterator>
#include<fstream>
#include<algorithm>
#include "structures.h"
#include "origenBuilder.h"

using namespace std;


// need to add units

double burnupcalc(isoInformation tempone, int N, double tolerance)
{

double mass = 1; //in kilograms
double BU_f, time_f; // burnup and time when k reaches one, time in days
double BU_total = 0;

int i =0;
while (i< tempone.neutron_prod.size())
{
    tempone.k_inf.push_back(tempone.neutron_prod[i]/tempone.neutron_dest[i]);
    i++;
}

i=0;

// at 4th entry time is at 460 days, one cycle time




while (tempone.k_inf[i] > 1)
    i++;            // finds the number of entry when k drops under 1


BU_f = tempone.BUd[i-1] + (1 - tempone.k_inf[i-1])*(tempone.BUd[i] - tempone.BUd[i-1])/(tempone.k_inf[i] - tempone.k_inf[i-1]);
// above line extrapolates the end point for burnup when k is one

time_f = tempone.time[i-1] + (1 - tempone.k_inf[i-1])*(tempone.time[i] - tempone.time[i-1])/(tempone.k_inf[i] - tempone.k_inf[i-1]);
// same as BU_f but for time (days)


int j=0;
while (j < i)
{
    BU_total = tempone.BUd[j] + BU_total;
    j++;
}

BU_total = BU_total + BU_f; // adds the last value when k reaches one

BU_total = 2*N*BU_total/(N+1);



return BU_total;

}



double enrichcalc(double BU_end, int N, double tolerance)
{

double X;
double BU_guess;
double BU2, BU7;
isoInformation test2;


// accurate guess calc, extrapolating from two data points at 2 and 7% enrichment
BU2 = burnupcalc(DataReader(test2, 0.02), N, 0.01);
BU7 = burnupcalc(DataReader(test2, 0.07), N, 0.01);

X = 0.02 + (BU_end - BU2)*(0.07 - 0.02)/(BU7 - BU2);

BU_guess = burnupcalc(DataReader(test2, X), N, 0.01);

// enrichment iteration
while (BU_end < BU_guess)
{
    X = X - 0.001;
    BU_guess = burnupcalc(DataReader(test2,X), N, 0.01 );
}


while (BU_end > BU_guess)
{
    X = X + 0.001;
    BU_guess = burnupcalc(DataReader(test2,X), N, 0.01 );
}
    return X;




}





























