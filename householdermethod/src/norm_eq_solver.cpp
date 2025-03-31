

#include "norm_eq_solver.h"

void simple_solver(mpreal norm, Eiseinstein& x3, int& state_found){

    if (norm >= 0){
    
        Eiseinstein result(0,0);
        mpreal sqrt_three = mpfr::sqrt(3.0);
        mpreal sqrt_norm = mpfr::sqrt(norm);

        mpreal qmax_tmp = sqrt_norm / sqrt_three;
        mpreal qmin = ceil(-2.0 * qmax_tmp) / 2.0;
        mpreal qmax = floor(2.0 * qmax_tmp) / 2.0;

        for (mpreal q = qmin; q <= qmax; q += 0.5){

            mpreal qsquared = q*q;

            //mpreal delta = 4 * (norm - 3 * qsquared);
            mpreal sqrt_delta = sqrt(4 * (norm - 3 * qsquared));

            mpreal p = 0.5 * sqrt_delta;
            
            if (floor(p + q) == p + q){
                mpreal a = static_cast<mpreal>(p + q);
                mpreal b = static_cast<mpreal>(2*q);
                x3 = Eiseinstein(a,b);
                state_found = 1;
                break;
            }

            p = -0.5 * sqrt_delta;
            if (floor(p + q) == p + q){
                mpreal a = static_cast<mpreal>(p + q);
                mpreal b = static_cast<mpreal>(2*q);
                x3 = Eiseinstein(a,b);
                state_found = 1;
                break;
            }

            if (state_found == 1) break;

        }
    }

}

