#ifndef NORM_EQ_SOLVER_H
#define NORM_EQ_SOLVER_H


//#include <gmp.h>
//#include <gmpxx.h>
#include "mpreal/mpreal.h"
#include "eiseinstein.h"

using namespace mpfr;
using namespace std;

void simple_solver(mpreal norm, Eiseinstein& x3, int& state_found);

#endif // NORM_EQ_SOLVER_H
