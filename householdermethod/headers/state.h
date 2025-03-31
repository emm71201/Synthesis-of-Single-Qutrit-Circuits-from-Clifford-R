// state.h

#ifndef STATE_H
#define STATE_H

#include <vector>
#include "eiseinstein.h"
#include "mpreal/mpreal.h"
//#include "matrix.h"

using namespace mpfr;
using namespace std;

class State {
private:
    vector<Eiseinstein> data; // Vector of length 3
    int f;

public:
    // Constructor
    State(const vector<Eiseinstein>& initialData, int initialF);
    
    // get data[i]
    Eiseinstein get(int i);
    int get_sde();
    mpreal eis_norm(int i);
    
    // set data[i]
    void set_eis(int i, Eiseinstein x);
    void set_sde(int sde);

    mpreal norm() const;

    // outer product of the state with itself
    //Matrix outer() const;

    void print();



};

#endif // STATE_H
