
#include "utils.h"
#include "cnumber_approximation.h"


int minimum_sde(double theta, double tol){
    cpp_dec_float_50 rhs = cpp_dec_float_50(1)/(2 * pow(tol, 2) * (1 - pow(tol,2)*0.25) * acos(1 - pow(tol,2)/2));
    cpp_dec_float_50 tmp = ceil(0.5 * (log(rhs)/log(sqrt(3.0)) + 1));

    return static_cast<int>(tmp);
}

vector<Eiseinstein> solve_norm_eq(long long N) {
    vector<Eiseinstein> results;
    double fourN = 4 * static_cast<double>(N);
    double sqrt4N3 = sqrt(fourN / 3);
    double lower_bound = ceil(-sqrt4N3);
    double upper_bound = floor(sqrt4N3);

    for (double y = lower_bound; y <= upper_bound; y++) {
        double delta = fourN - 3 * y * y;
        double x1 = (y + sqrt(delta)) / 2;
        double x2 = (y - sqrt(delta)) / 2;

        if (fmod(x1, 1.0) == 0.0) {
            Eiseinstein tmp(static_cast<long long>(x1), static_cast<long long>(y));
            results.push_back(tmp);
        }
        if (fmod(x2, 1.0) == 0.0) {
            Eiseinstein tmp(static_cast<long long>(x2), static_cast<long long>(y));
            results.push_back(tmp);
        }
    }
    return results;
}

void solution(double theta, int k, double tol, vector<Eiseinstein>& matrix){
    bool matrix_found = false;

    // setup
    cpp_dec_float_50 three = 3.0;
    cpp_dec_float_50 three_sqrt = sqrt(three);
    cpp_dec_float_50 radius_two = pow(three, ceil(k/2));
    cpp_dec_float_50 radius_one = radius_two * (1 - pow(tol,2)/2);
    double radius_two_double = pow(sqrt(3.0), k);
    double tol_sq = tol*tol;
    long long threeToK = static_cast<long long>(pow(3,k));
    Eiseinstein factor = Eiseinstein(2,1).pow(k);
    vector<Eiseinstein> units = {Eiseinstein(1,0),Eiseinstein(-1,0),Eiseinstein(0,1),Eiseinstein(0,-1),Eiseinstein(1,1),Eiseinstein(-1,-1)};


    vector<Eiseinstein> x1_vector = approximate_cnumber(-theta/2, k, tol, three, three_sqrt, radius_two, radius_one);
    long long x1_vector_size = x1_vector.size();

    // We will try to save some computation of the norm equation
    map<long long, vector<Eiseinstein>> solution_normeq;

    #pragma omp parallel for
    for (int i = 0; i < x1_vector_size; i++){
        if (matrix_found) continue;
        Eiseinstein x1 = x1_vector[i];
        Eiseinstein x1_conj = x1.conjugate();
        long long x1_norm = x1.norm();
        long long eta1 = threeToK - x1_norm; // Used in the quadratic equation  eta |x2|^4 + (|beta|^2 - |alpha|^2 - eta^2)|x2|^2 + eta * |alpha|^2 = 0
        double x1_norm_sqrt = sqrt(static_cast<double>(x1_norm));


        double err1_sq = compute_error_sq(-theta/2, k, x1);
        double tol2 = sqrt(tol_sq - err1_sq);
        
        vector<Eiseinstein> y2_vector = approximate_cnumber(theta/2, k, tol2, three, three_sqrt, radius_two, radius_one);
        long long y2_vector_size = y2_vector.size();


        for (int j = 0; j < y2_vector_size; j++){
            if (matrix_found) continue;
            Eiseinstein y2 = y2_vector[j];
            Eiseinstein y2_conj = y2.conjugate();
            long long y2_norm = y2.norm();
            long long eta2 = threeToK - y2_norm;
            double y2_norm_sqrt = sqrt(static_cast<double>(y2_norm));


            double err2_sq = compute_error_sq(theta/2, k, y2);
            if (sqrt(err2_sq) > tol2) continue;

            double tol3 = sqrt(tol_sq - err1_sq - err2_sq);

            
            vector<Eiseinstein> z3_vector = approximate_cnumber(0, k, tol3, three, three_sqrt, radius_two, radius_one);
            long long z3_vector_size = z3_vector.size();


            for (int kz = 0; kz < z3_vector_size; kz++){
                if (matrix_found) continue;
                Eiseinstein z3 = z3_vector[kz];
                Eiseinstein z3_conj = z3.conjugate();
                long long z3_norm = z3.norm();
                long long eta3 = threeToK - z3_norm;
                double z3_norm_sqrt = sqrt(static_cast<double>(z3_norm));

                double err3_sq = compute_error_sq(0, k, z3);                

                // Horn conditions
                if (x1_norm_sqrt + y2_norm_sqrt - z3_norm_sqrt > radius_two_double) continue;
                if (x1_norm_sqrt + z3_norm_sqrt - y2_norm_sqrt > radius_two_double) continue;
                if (y2_norm_sqrt + z3_norm_sqrt - x1_norm_sqrt > radius_two_double) continue;

                for (int ia=0; ia < 6; ia++){
                    if (matrix_found) continue;
                    Eiseinstein alpha = x1_conj * y2_conj - units[ia] * factor * z3;
                    long long alpha_norm = alpha.norm();
                    if (alpha_norm > eta1 * eta2 || alpha_norm > eta1*eta1 || alpha_norm > eta2*eta2) continue;

                    //for (int ib=0; ib < 6; ib++){
                    if (matrix_found) continue;
                    Eiseinstein beta = x1_conj * z3_conj - units[ia] * factor * y2;
                    long long beta_norm = beta.norm();
                    if (beta_norm > eta1*eta3 || beta_norm > eta1*eta1 || beta_norm > eta3*eta3) continue;

                    // Get also gamma
                    Eiseinstein gamma = y2_conj * z3_conj - units[ia] * factor * x1;
                    long long gamma_norm = gamma.norm();
                    if (gamma_norm > eta2 * eta3 || gamma_norm > eta2*eta2 || gamma_norm > eta3*eta3) continue;

                    // Build the quadratic equation
                    long long a1 = eta1;
                    long long a2 = beta_norm - alpha_norm - eta1*eta1;
                    long long a3 = eta1 * alpha_norm;
                    double discr = a2*a2 - 4 * a1 * a3;

                    
                    if (discr < 0) continue;
                    double val1 = (-a2 + sqrt(discr))/(2*a1);
                    double val2 = (-a2 - sqrt(discr))/(2*a1);
                    vector<long long> x2_norms;
                    if (floor(val1) == val1 && val1 > 0 && val1 <= eta1 && val1 <= eta2) x2_norms.push_back(val1);
                    if (floor(val2) == val2 && val2 > 0 && val2 <= eta1 && val2 <= eta2) x2_norms.push_back(val2);
                    for (int ix2_n=0; ix2_n<x2_norms.size();ix2_n++){
                        if (matrix_found) continue;
                        long long x2_norm = x2_norms[ix2_n];
                        long long x3_norm = eta1 - x2_norm;
                        long long z2_norm = eta2 - x2_norm;
                        long long z1_norm = eta3 - z2_norm;
                        long long y1_norm = eta1 - z1_norm;
                        long long y3_norm = eta2 - y1_norm;


                        auto x2_iter = solution_normeq.find(x2_norm);
                        if (x2_iter != solution_normeq.end()){
                            vector<Eiseinstein> x2_vector = x2_iter->second;
                        } else{
                            vector<Eiseinstein> x2_vector = solve_norm_eq(x2_norm);
                            #pragma omp critical
                            {
                                solution_normeq[x2_norm] = x2_vector;
                            }
                            
                        }
                        auto x3_iter = solution_normeq.find(x3_norm);
                        if (x3_iter != solution_normeq.end()){
                            vector<Eiseinstein> x3_vector = x3_iter->second;
                        } else{
                            vector<Eiseinstein> x3_vector = solve_norm_eq(x3_norm);
                            #pragma omp critical
                            {
                                solution_normeq[x3_norm] = x3_vector;
                            }
                            
                        }


                        vector<Eiseinstein> x2_vector = solve_norm_eq(x2_norm);
                        vector<Eiseinstein> x3_vector = solve_norm_eq(x3_norm);

                        for (int ix2=0; ix2 < x2_vector.size(); ix2++){
                            if (matrix_found) continue;
                            Eiseinstein x2 = x2_vector[ix2];
                            Eiseinstein x2_conj = x2.conjugate();
                            auto y1_conj_data = alpha/x2_conj;
                            if (y1_conj_data.second.norm() != 0) continue;
                            //Eiseinstein y1 = y1_conj_data.first.conjugate();
                            Eiseinstein y1_conj = y1_conj_data.first;
                            Eiseinstein y1 = y1_conj.conjugate();

                            Eiseinstein tmpy3 = -(x1_conj*y1 + x2_conj*y2);
                                

                            for (int ix3=0; ix3 < x3_vector.size(); ix3++){
                                if (matrix_found) continue;
                                Eiseinstein x3 = x3_vector[ix3];
                                Eiseinstein x3_conj = x3.conjugate();
                                auto z1_conj_data = beta/x3_conj;
                                if (z1_conj_data.second.norm() != 0) continue;
                                Eiseinstein z1 = z1_conj_data.first.conjugate();
                                auto y3_data = tmpy3/x3_conj;
                                if (y3_data.second.norm() != 0) continue;
                                Eiseinstein y3 = y3_data.first;
                                Eiseinstein y3_conj = y3.conjugate();

                                //for (int ig = 0; ig < 6; ig++){
                                //if (matrix_found) continue;
                                //Eiseinstein gamma = y2_conj * z3_conj - units[ia] * factor * x1;
                                //if (gamma.norm() != y3_norm * z2_norm) continue;


                                auto z2_conj_data = gamma/y3_conj;
                                if (z2_conj_data.second.norm() != 0) continue;
                                Eiseinstein z2 = z2_conj_data.first.conjugate();

                                // Final check
                                Eiseinstein cross1 = x1_conj*y1 + x2_conj*y2 + x3_conj*y3;
                                Eiseinstein cross2 = x1_conj*z1 + x2_conj*z2 + x3_conj*z3;
                                Eiseinstein cross3 = y1_conj*z1 + y2_conj*z2 + y3_conj*z3;

                                if (cross1.norm()!=0 || cross2.norm() != 0 || cross3.norm() != 0) continue; 
                                matrix_found = true;

                                // Now compute the final error
                                //double err3_sq = compute_error_sq(0, k, z3);
                                double final_error = sqrt(err1_sq + err2_sq + err3_sq);
                                #pragma omp critical
                                {

                                    cout << tol << " " << tol2 << " " << tol3 << endl;

                                    cout << sqrt(err1_sq) << " " << sqrt(err2_sq) << " " << sqrt(err3_sq) << endl;

                                    matrix = {x1,x2,x3,y1,y2,y3,z1,z2,z3};
                                    cout << theta << "," << k << "," << tol << "," << final_error << endl;
                                }
                                
                            }

                        }

                    }
                }
                
            }
        }
    }
}





