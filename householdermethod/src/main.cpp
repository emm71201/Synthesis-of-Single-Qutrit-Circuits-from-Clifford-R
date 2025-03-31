#include <cmath>
#include "mpreal/mpreal.h"
#include <vector>
#include <utility>
#include <iomanip>
#include "eiseinstein.h"
#include "approximate_cnumber.h"
#include "enumeration.h"
#include "state_approximation.h"
#include "state.h"
//#include "matrix.h"
#include "norm_eq_solver.h"
using namespace mpfr;
using namespace std;

// the state we are approximating is: 1/sqrt(2) (cos theta/2, sin theta/2, -1, 0)
// set alpha = theta/2
void solution(mpreal theta, mpreal tol, mpreal contraction){

    mpreal sqrt_two = mpfr::sqrt(2.0);
    mpreal state_tol = tol / (2.0 * sqrt_two * contraction);
    mpreal alpha = 0.5 * theta;

    vector<Eiseinstein> statetmp = {Eiseinstein(0,0), Eiseinstein(0,0), Eiseinstein(0,0)};
    State state(statetmp, 0);

    approximate_state(alpha, state_tol, tol, state);

    mpreal error_analytic = error(alpha, state);

    // get the unit vector to for the exact synthesis
    int f = state.get_sde();
    Eiseinstein sign(1,0);
    if (f % 2 != 0) sign = Eiseinstein(-1, 0);
    mpreal maxnorm = mpfr::pow(3.0, f);
    mpreal tmp_x1 = maxnorm - 2 * state.eis_norm(0);
    Eiseinstein x1 =  sign * Eiseinstein(tmp_x1, 0);
    Eiseinstein x2 = sign * Eiseinstein(-2.0, 0) * state.get(0).conjugate() * state.get(1);
    Eiseinstein x3 = sign * Eiseinstein(-2.0, 0) * state.get(0).conjugate() * state.get(2);


    //display the entries of the first column and the sde: x1.a, x1.b, x2.a, x2.b, x3.a, x3.b and 2f
    //cout << "answer: " << x1.a << "," << x1.b << ","<< x2.a << "," << x2.b << "," << x3.a << "," << x3.b << "," << 2 * f << endl;
    //display error results. The sde shown is multiplied by 2
    cout << std::fixed << std::setprecision(10) << theta << "," << std::defaultfloat << 2 * f << "," << tol << "," << error_analytic <<  endl;

}


int main(int argc, char* argv[]) {
    // Check if the user provided the required arguments
    if (argc != 4) {
        std::cerr << "Usage: ./main <theta> <tol>" << std::endl;
        return 1;
    }

    mpfr_prec_t MAX_NBITS = 512; // we start with a very high precision : about 155 digits
    mpreal::set_default_prec(MAX_NBITS);
    // get tol
    mpreal tol =  mpreal(argv[2]); 
    // contraction is a factor between 0 and 1 to relatx the 2 * sqrt(2) in the state error
    mpreal contraction = mpreal(argv[3]); 

    // compute the number of bits we actually need: the numbers we will need won't exceed about O(1/epsilon^2) but I will use O(1/epsilon^3) because 
    // classical memory is not an issue
    //mpreal NBITS_temp = abs(-3.0* log(tol)/log(2.0));
    //mpfr_prec_t NBITS = static_cast<mpfr_prec_t>( min(mpreal::get_default_prec(), NBITS_temp) );
    //mpreal::set_default_prec(NBITS);


    //tol.setPrecision(NBITS); // adjust down the precision for tol to what we need
    mpreal theta = mpreal(argv[1]);
    
    // Call the solution function with the parsed values
    solution(theta, tol, contraction);    

    return 0;
}

