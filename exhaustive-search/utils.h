#ifndef UTILS_H
#define UTILS_H

#include "includes.h" 

int minimum_sde(double theta, double tol);
vector<Eiseinstein> solve_norm_eq(long long N);
void solution(double theta, int k, double tol, vector<Eiseinstein>& matrix);

#endif // UTILS_H
