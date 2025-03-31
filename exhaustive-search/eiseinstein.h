#ifndef EISEINSTEIN_H
#define EISEINSTEIN_H

#include "includes.h" 
using namespace std;

class Eiseinstein {
public:
    long long a;
    long long b;

    Eiseinstein(long long aval, long long bval);
    long long norm() const;
    Eiseinstein conjugate() const;
    Eiseinstein operator +(const Eiseinstein& other);
    Eiseinstein operator -(const Eiseinstein& other);
    Eiseinstein operator *(const Eiseinstein& other);
    Eiseinstein operator-() const;
    Eiseinstein pow(int k);
    pair<Eiseinstein, Eiseinstein> operator /(const Eiseinstein& other);
    pair<Eiseinstein, Eiseinstein> operator /(const long long& other);
    //pair<double, double> value();
    //double LocalizedDotProduct(double alpha, int k);
    int get_sde(int exponent);
    void print() const;
};

#endif // EISEINSTEIN_H
