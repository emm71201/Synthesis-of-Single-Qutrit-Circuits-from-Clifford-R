#ifndef EISEINSTEIN_H
#define EISEINSTEIN_H



//#include <gmp.h>
//#include <gmpxx.h>
#include "mpreal/mpreal.h"
#include <iostream>
#include <complex>

using namespace mpfr;
using namespace std;

class Eiseinstein {
public:
    mpreal a;
    mpreal b;

    Eiseinstein();
    Eiseinstein(mpreal aval, mpreal bval);
    mpreal norm() const;
    Eiseinstein conjugate() const;
    Eiseinstein operator +(const Eiseinstein& other);
    Eiseinstein operator -(const Eiseinstein& other);
    Eiseinstein operator-() const;
    Eiseinstein operator *(const Eiseinstein& other);
    Eiseinstein pow(int k);
    pair<Eiseinstein, Eiseinstein> operator /(const Eiseinstein& other);
    pair<Eiseinstein, Eiseinstein> operator /(const mpreal& other);
    complex<mpreal> value();

    void print() const;
};

#endif // EISEINSTEIN_H
