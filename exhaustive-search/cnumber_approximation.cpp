
#include "cnumber_approximation.h"

double compute_error_sq(double alpha, int k, Eiseinstein x){
    double threetok = pow(3.0, k);
    Eiseinstein factor = x * Eiseinstein(-1,-2).pow(k);
    double q = 0.5*factor.b;
    double p = factor.a - q;
    //double Re = p * cos(alpha) + q * sqrt(3.0) * sin(alpha);
    double Re = p * cos(alpha) + q * 1.73205080757 * sin(alpha);

    return 2 * (1 - Re/threetok);
}

double compute_error_sq_precompute(double alpha, int k, Eiseinstein factor,double threetok){
    double q = 0.5*factor.b;
    double p = factor.a - q;
    //double Re = p * cos(alpha) + q * sqrt(3.0) * sin(alpha);
    double Re = p * cos(alpha) + q * 1.73205080757 * sin(alpha);
    return 2 * (1 - Re/threetok);
}

double compute_error_sq_zero_precompute(int k, Eiseinstein factor,double threetok){
    double p = factor.a - 0.5*factor.b;

    return 2 * (1 - p/threetok);
}



// The rotation angle: -pi/2 <= theta <= pi/2
// alpha = theta/2   ---->   -pi/4 <= alpha <= pi/4. 
// We will also assume  tol < sqrt(2)  ------> otherwise, there is a trivial solution (I need to confirm it)
vector<Eiseinstein> approximate_cnumber(double alpha, int k, double tol, cpp_dec_float_50 three, cpp_dec_float_50 three_sqrt, cpp_dec_float_50 radius_two, cpp_dec_float_50 radius_one){

    int f = ceil(k/2.0);
    int parity = k % 2;
    Eiseinstein factor = Eiseinstein(-1, 0).pow(f) * Eiseinstein(1,2).pow(parity);
    vector<Eiseinstein> results;
    
    cpp_dec_float_50 delta = acos(1 - pow(tol,2.0)/2.0);
    cpp_dec_float_50 phi1 = alpha - delta;
    cpp_dec_float_50 phi2 = alpha + delta;
    cpp_dec_float_50 xmin, xmax, ymin, ymax, qmin, qmax;
    cpp_dec_float_50 x_intersect1 = radius_two * cos(phi1);
    cpp_dec_float_50 x_intersect2 = radius_two * cos(phi2);
    cpp_dec_float_50 x_intersect_max = max(x_intersect1, x_intersect2);

    if (alpha == 0){
        xmin = radius_one;
        xmax = radius_two;
    } else 
    {
        tie(xmin, xmax) = minmax(x_intersect1, x_intersect2);
        // check if (r2,0) works and adjust xmax accordingly
        if (radius_two*cos(alpha) >= radius_one) xmax = radius_two;
    }

    cpp_dec_float_50 half = 0.5;
    for (cpp_dec_float_50 p = ceil(2.0 * xmin) / 2.0; p <= floor(2.0 * xmax) / 2.0; p += half) {
        if (alpha == 0){
            ymin = -sqrt(radius_two * radius_two - p * p);
            ymax = sqrt(radius_two * radius_two - p * p);
        }
        else {

            if (p <= x_intersect_max){

                if (alpha > 0){
                    ymin = (radius_one - p * cos(alpha)) / sin(alpha);
                    ymax = sqrt(radius_two * radius_two - p * p);
                }
                else {
                    ymin = -sqrt(radius_two * radius_two - p * p);
                    ymax = (radius_one - p * cos(alpha)) / sin(alpha);
                }
            } else {

                ymin = -sqrt(radius_two * radius_two - p * p);
                ymax = sqrt(radius_two * radius_two - p * p);
            }

        }

        qmin = ceil(2 * ymin / three_sqrt) / 2;
        qmax = floor(2 * ymax / three_sqrt) / 2;


        if (abs(static_cast<long long>(2 * qmin) % 2) != abs(static_cast<long long>(2 * p) % 2)) qmin += half ;
        if (abs(static_cast<long long>(2 * qmax) % 2) != abs(static_cast<long long>(2 * p) % 2)) qmax -= half;

        for (cpp_dec_float_50 q = qmin; q <= qmax; q++) {

            Eiseinstein x_tmp = Eiseinstein(static_cast<long long>(p + q), static_cast<long long>(2*q));
            auto x_data = x_tmp/factor;
            if (x_data.second.norm()==0) {
                Eiseinstein x = x_data.first;
                //results[x.norm()].push_back(x);
                results.push_back(x);
            }
        }

    }

    // Sorting according to the error
    sort(results.begin(), results.end(), [alpha, k](const Eiseinstein& a, const Eiseinstein& b) {
        return compute_error_sq(alpha, k, a) < compute_error_sq(alpha, k, b);
    });

    return results;
}


map<long long, vector<Eiseinstein>> fincke_pohst(int k, double tol, double threetok){

    double maxnorm = floor(threetok * pow(tol, 2) * (1 - pow(tol/2.0, 2)));

    double sqrt_maxnorm = sqrt(maxnorm);

    double maxp = floor(2.0 * sqrt_maxnorm)/2.0;

    map<long long, vector<Eiseinstein>> results;


    for (double p = -maxp; p <= maxp; p += 0.5){

        double psqure = p*p;

        double maxq = floor(2.0 * sqrt((maxnorm - psqure)/3.0))/2.0;

        if (static_cast<long long>(2 * maxq) % 2 != abs(static_cast<long long>(2 * p) % 2)) maxq -= 0.5;

        cout << p << " " << -maxq << endl;


        for (double q = -maxq; q <= maxq; q++){

            double qsquare = q*q;
            long long norm = static_cast<long long>(psqure + 3*qsquare);

            //cout << p << " " << q << " " << norm << endl;

            results[norm].push_back(Eiseinstein(static_cast<long long>(p+q), static_cast<long long>(2*q)));

        }


    }

    return results;

}


