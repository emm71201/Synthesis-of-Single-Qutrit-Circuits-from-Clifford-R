#ifndef CNUMBER_APPROXIMATION_H
#define CNUMBER_APPROXIMATION_H

#include "includes.h" // All necessary includes are handled here

double compute_error_sq(double alpha, int k, Eiseinstein x);
double compute_error_sq_precompute(double alpha, int k, Eiseinstein x, double threetok);
double compute_error_sq_zero_precompute(int k, Eiseinstein x, double threetok);
vector<Eiseinstein> approximate_cnumber(double alpha, int k, double tol, cpp_dec_float_50 three, cpp_dec_float_50 three_sqrt, cpp_dec_float_50 radius_two, cpp_dec_float_50 radius_one);
map<long long, vector<Eiseinstein>> fincke_pohst(int k, double tol, double threetok);

#endif // CNUMBER_APPROXIMATION_H
