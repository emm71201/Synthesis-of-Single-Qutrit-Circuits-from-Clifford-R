// state_approximation.h

#ifndef STATE_APPROXIMATION_H
#define STATE_APPROXIMATION_H

#include <complex>
#include <iostream>
#include <iomanip>
#include "mpreal/mpreal.h"
#include "eiseinstein.h"
#include "enumeration.h"
#include "state.h"
#include "norm_eq_solver.h"

using namespace mpfr;
using namespace std;

mpreal error(mpreal alpha, State state);
// Function declaration for approximate_state
void approximate_state(mpreal alpha, mpreal state_tol, mpreal matrix_tol, State& state);

#endif // STATE_APPROXIMATION_H
