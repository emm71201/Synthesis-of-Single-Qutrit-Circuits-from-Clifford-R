
// this is state_approximation.cpp
// generate the header file 

#include "state_approximation.h"
#include "state.h"

mpreal error(mpreal alpha, State state){

    int f = state.get_sde();
    mpreal maxnorm = mpfr::pow(3.0, f);
    mpreal cosalpha = cos(alpha);
    mpreal sinalpha = sin(alpha);
    mpreal rex1 = state.get(0).value().real();
    mpreal imx1 = state.get(0).value().imag();
    mpreal rex2 = state.get(1).value().real();
    mpreal imx2 = state.get(1).value().imag();

    mpreal re_temp = cosalpha * rex1 + sinalpha * imx1 - rex2;
    mpreal im_temp = cosalpha * imx1 - sinalpha * rex1 - imx2;

    mpreal temp_sq = re_temp * re_temp + im_temp * im_temp;
    temp_sq = 0.5 * temp_sq / maxnorm;
    return sqrt( 8 * (1 - temp_sq) );
}


void approximate_state(mpreal alpha, mpreal state_tol, mpreal matrix_tol, State& state){

    int statefound = 0;
    mpreal maxnorm = 1;
    Eiseinstein chi(1,2);
    // setting up the denominator exponent
    int f = 0;
    int ceil_fhalf = 0;
    //int floor_fhalf = 0;
    int parity = 0;
    Eiseinstein factor(1,0);
    Eiseinstein x3(0,0);


    mpreal x3_norm;

    mpreal r1_initial = 1 - 0.5 * state_tol * state_tol;
    mpreal r2 = 1;
    mpreal r1 = r1_initial;


    while (statefound == 0){

        
        vector<pair<Eiseinstein, Eiseinstein>> candidates = enumerate(alpha, r1, r2, factor);
        int cand_size = candidates.size();

        for (int ii = 0; ii < cand_size; ii++){
            if (statefound) break;
            state.set_sde(f);
            state.set_eis(0, candidates[ii].first);
            state.set_eis(1, candidates[ii].second);
            // crucially, we check the final error 
            mpreal final_error = error(alpha, state);
            if (final_error > matrix_tol) continue;


            mpreal x3_norm = maxnorm - state.eis_norm(0) - state.eis_norm(1);
            simple_solver(x3_norm, x3, statefound);
        
            if (statefound == 1){
                state.set_sde(f);
                state.set_eis(2, x3);
            }

        }

        // updating f and all the settings
        f += 1;
        ceil_fhalf = ceil(f/2.0);
        //floor_fhalf = floor(f/2.0);
        parity = (parity + 1) % 2;
        if (ceil_fhalf % 2 == 0) {
            factor = chi.pow(parity);
        } else {
            factor = -chi.pow(parity);
        }
        maxnorm = 3.0 * maxnorm;
        r2 = pow(3.0, ceil_fhalf);
        r1 = r2 * r1_initial;

        // r2 = pow(3.0, floor_fhalf);
        // r1 = pow(3.0, ceil_fhalf) * r1_initial;
        
    }

    
}