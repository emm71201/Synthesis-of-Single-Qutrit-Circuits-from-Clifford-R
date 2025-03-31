// state.cpp

#include "state.h"

// Constructor
State::State(const vector<Eiseinstein>& initialData, int initialF) : data(initialData), f(initialF) {
    if (data.size() != 3) {
        throw invalid_argument("Data vector must have exactly 3 elements.");
    }
}

Eiseinstein State::get(int i){
    if (i >= 3) {
        throw invalid_argument("argument must be < 3");
    }
    return data[i];
}

int State::get_sde(){
    return f;
}
mpreal State::eis_norm(int i){
    return data[i].norm();
}

void State::set_eis(int i, Eiseinstein x){

    if (i >= 3) {
        throw invalid_argument("argument must be < 3");
    }

    data[i] = x;

}

void State::set_sde(int sde){
    f = sde;
}



mpreal State::norm() const {

    mpreal res = 0;

    for (int i = 0; i < 3; i++){

        res += static_cast<mpreal>(data[i].norm());

    }

    return res / pow(3.0, f);

}

// Matrix State::outer() const {

//     Matrix result;
//     mpreal denominator = pow(3.0, f);

//     for (int i = 0; i < 3; i++){

//         for (int j = 0; j < 3; j++){

//             Eiseinstein val2 = data[i].conjugate() * data[j];

//             result.setElement(i, j, val2.value() / denominator);

//         }
//     }

//     return result;

// }

void State::print(){
    for (int i = 0; i < 3; i++){
        data[i].print();
    }
}
