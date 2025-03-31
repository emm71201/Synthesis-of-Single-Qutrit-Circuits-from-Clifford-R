
#include "includes.h"
#include "cnumber_approximation.h"
#include "utils.h"
using namespace std;

void process(double theta, double tol){
    vector<Eiseinstein> matrix;
    int k = minimum_sde(theta, tol);
    //cout << "CURRENT DENOMINATOR: " << k <<endl;
    solution(theta, k, tol, matrix);
    while (matrix.empty()){
        k += 1;
        //cout << "CURRENT DENOMINATOR: " << k <<endl;
        solution(theta, k, tol,matrix);
    }
}

int main(int argc, char* argv[]) {
    //Check if the correct number of command-line arguments is given
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <theta> <tol>" << endl;
        return 1;  // Return error code
    }

    // Convert command-line arguments to double
    double theta = atof(argv[1]);
    double tol = atof(argv[2]);
    //int k = atof(argv[3]);
    //double threetok = pow(3.0, k);

    process(theta, tol);

    // map<long long, vector<Eiseinstein>> fp = fincke_pohst(k, tol, threetok);

    // for (auto it : fp){
    //     cout << "Norm: " << it.first << endl;
    //     for (size_t i = 0; i < it.second.size(); i++){
    //         it.second[i].print();
    //     }
    // }


    



	return 0;
}
