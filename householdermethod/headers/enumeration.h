#ifndef ENUMERATION_H
#define ENUMERATION_H

#include <iostream>
#include <vector>
#include <utility>
#include "mpreal/mpreal.h"
#include "eiseinstein.h"
#include "approximate_cnumber.h"
using namespace std;
using namespace mpfr;

vector<pair<Eiseinstein, Eiseinstein>> enumerate(mpreal alpha, mpreal r1, mpreal r2, Eiseinstein factor);

#endif // ENUMERATION_H
