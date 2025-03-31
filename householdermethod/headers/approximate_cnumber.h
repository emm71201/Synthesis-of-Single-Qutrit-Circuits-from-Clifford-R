// approximate_Cnumber.h

#ifndef APPROXIMATE_CNUMBER_H
#define APPROXIMATE_CNUMBER_H

//#include <gmp.h>
//#include <gmpxx.h>
#include "mpreal/mpreal.h"
#include <cassert>
#include "eiseinstein.h"
#include <iostream>
#include <vector>
#include <tuple>
using namespace mpfr;

// Function declaration for approx_cnumber
vector<Eiseinstein> approx_cnumber(mpreal alpha, mpreal lambda1, mpreal lambda2);

#endif // APPROXIMATE_CNUMBER_H
