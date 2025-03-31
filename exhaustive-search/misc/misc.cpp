
#include "includes.h"

map<long long, vector<Eiseinstein>> fp_enumeration(long long N){
    // return a dictionary of Eiseinstein integer up to norm N
    // The dictionary has the {norm:vector of Eiseinstein integers}
    map<long long, vector<Eiseinstein>> results;
    double qmin, qmax;

    double half = 0.5;
    double p_bound = 2*sqrt(static_cast<double>(N));

    for (double p = ceil(-p_bound)/2; p <= floor(p_bound)/2; p+= 0.5 ){

        double q_bound = 2*sqrt((N - p*p)/3);
        
        qmin = ceil(-q_bound)/2;
        qmax = floor(q_bound)/2;

        if (abs(static_cast<long long>(2*qmin) % 2) != abs(static_cast<long long>(2*p) % 2)){
            qmin += half;
        }
        if (abs(static_cast<long long>(2*qmax) % 2) != abs(static_cast<long long>(2*p) % 2)){
            qmax -= half;
        }

        for (double q = qmin; q <= qmax; q++){

            long long norm = static_cast<long long>(p*p + 3*q*q);
            results[norm].push_back(Eiseinstein(static_cast<long long>(p + q), static_cast<long long>(2*q)));

        }

    }

    return results;

}

// compute Re(x e^(i theta))
double product_realpart(Eiseinstein x, double theta){
    double p = x.a - 0.5*x.b;
    double q = 0.5 * x.b;
    double result = p * cos(theta) - q * sqrt(3.0) * sin(theta);
    return result;
}