
#include "approximate_cnumber.h"

// alpha in first or 4th quadrant --> cosalpha >= 0
vector<Eiseinstein> approx_cnumber(mpreal alpha, mpreal lambda1, mpreal lambda2){

    assert((lambda2 >= lambda1) && "Invalid values in approximate cnumber ");

    vector<Eiseinstein> results;

    mpreal lambda2_squared = lambda2 * lambda2;
    mpreal lambda1_squared = lambda1 * lambda1;
    mpreal cosalpha = cos(alpha);
    mpreal sinalpha = sin(alpha);

    mpreal sqrt_three = mpfr::sqrt(3.0);
    mpreal delta = acos(lambda1 / lambda2);
    mpreal phi1 = alpha - delta;
    mpreal phi2 = alpha + delta;
    mpreal xmin, xmax, ymin, ymax, qmin, qmax;
    mpreal x_intersect1 = lambda2 * cos(phi1);
    mpreal x_intersect2 = lambda2 * cos(phi2);
    mpreal x_intersect_max = max(x_intersect1, x_intersect2);

    if (alpha == 0){
        xmin = lambda1;
        xmax = lambda2;
    } else 
    {
        tie(xmin, xmax) = minmax(x_intersect1, x_intersect2);
        // check if (r2,0) works and adjust xmax accordingly
        if (lambda2*cosalpha >= lambda1) xmax = lambda2;
    }

    mpreal half = 0.5;

    for (mpreal p = ceil(2.0 * xmin) / 2.0; p <= floor(2.0 * xmax) / 2.0; p += half) {
        if (alpha == 0){
            ymin = -sqrt(lambda2_squared - p * p);
            ymax = sqrt(lambda2_squared - p * p);
        }
        else {

            if (p <= x_intersect_max){

                if (alpha > 0){
                    ymin = (lambda1 - p * cosalpha) / sinalpha;
                    ymax = sqrt(lambda2_squared - p * p);
                }
                else {
                    ymin = -sqrt(lambda2_squared - p * p);
                    ymax = (lambda1 - p * cosalpha) / sinalpha;
                }
            } else {

                ymin = -sqrt(lambda2_squared - p * p);
                ymax = sqrt(lambda2_squared - p * p);
            }

        }

        qmin = ceil(2 * ymin / sqrt_three) / 2;
        qmax = floor(2 * ymax / sqrt_three) / 2;

        if (floor(p + qmin) != p + qmin) qmin += half;


        //if (abs(static_cast<long long>(2 * qmin) % 2) != abs(static_cast<long long>(2 * p) % 2)) qmin += half ;
        //if (abs(static_cast<long long>(2 * qmax) % 2) != abs(static_cast<long long>(2 * p) % 2)) qmax -= half;

        for (mpreal q = qmin; q <= qmax; q++) {

            #pragma omp critical
            {
                results.push_back(Eiseinstein(p+q, 2*q));
            }
            
        }

    }
    

    return results;
}