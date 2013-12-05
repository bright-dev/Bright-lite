#include<iostream>
#include<vector>
#include<regex>
#include<iterator>
#include<fstream>
#include<algorithm>
#include "structures.h"
#include "origenBuilder.h"

using namespace std;

double intpol(double y0, double y1, double x0, double x1, double x);

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

double burnupcalc(isoInformation tempone, int N, double tolerance)
{

double mass = 1; //in kilograms
double BU_f, time_f; // burnup and time when k reaches one, time in days
double BU_total = 0;
double BU_finder = 0;
double BU_n = 0;    // estimated burnup of the n th batch
double time_finder = 0;
double k_total = 0;
double leakage = 0.99; // neutron non-leakage probability
int i = 0;
double k_batch [N]; // the k-value of each batch
double t_batch;
double x0 = 0;
double x1 = 0;
int m = 0;
while (i< tempone.neutron_prod.size()) // builds the k-values from neutron production and destruction
{
    tempone.k_inf.push_back(tempone.neutron_prod[i]*leakage/tempone.neutron_dest[i]);
    i++;
}

i=0;




while (tempone.k_inf[i] > 1)
    i++;            // finds the number of entry when k drops under 1


BU_f = intpol(tempone.BUd[i-1],tempone.BUd[i],tempone.k_inf[i-1],tempone.k_inf[i],1.0);

time_f = intpol(tempone.time[i-1],tempone.time[i],tempone.k_inf[i-1],tempone.k_inf[i],1);


int j=0;
while (j < i)
{
    BU_total = tempone.BUd[j] + BU_total;
    j++;
}

BU_total = BU_total + BU_f; // adds the last value when k reaches one


if (N == 1)
    return BU_total; // returns the single batch BU in the case of single batch

BU_total = 2*N*BU_total/(N+1); //linear approximation of mutli batch burnup, used to find a good initial guess of the max burnup

while (1)
{

    j = 0;
    while (j<N)
    {
        BU_finder = 0;
        BU_n = BU_total*(j+1)/N;

            i=0;
                    while (BU_finder + tempone.BUd[i] < BU_n)  // finds the discrete point i corresponding to the burnup
                {
                    BU_finder = BU_finder + tempone.BUd[i];
                    i++;
                }

        while (m < i)
        {
            x0 = x0 + tempone.BUd[m];  // sums the burnup for total burnup
            m++;
        }
        x1 = x0 + tempone.BUd[m];     // adds on more discrete point for linear interpolation

        k_batch[j] = intpol(tempone.k_inf[i], tempone.k_inf[i+1], x0, x1, BU_n);  //finds the k of the batch


        j++;
        x0 = 0;
        x1= 0;
        m = 0;
    }
    j = 0;


    while (j < N) //sums the k values of every batch
    {
        k_total = k_total + k_batch[j];
        j++;
    }
    k_total = k_total/N;

    //cout << k_total << endl;

    if (abs(1 - k_total) < 0.0001 ) //breaks out of loop if k is close enough, tolerance value passed to the function can be used here
    break;


    BU_total = intpol(0,BU_total,tempone.k_inf[0],k_total,1); // updates the guess using (k(0),0) and (k_total, BU_total)



    k_total = 0;

}

return BU_total;

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



double intpol(double y0, double y1, double x0, double x1, double x)
{
    // linear interpolation function
    double y = y0 + (y1 - y0)*(x - x0)/(x1 - x0);
    return y;
}

























