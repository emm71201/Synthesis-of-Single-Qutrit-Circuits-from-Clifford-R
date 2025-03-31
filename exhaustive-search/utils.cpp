
#include "utils.h"
#include "cnumber_approximation.h"


int minimum_sde(double theta, double tol){
    cpp_dec_float_50 tol_sqrt = pow(tol,2);
    cpp_dec_float_50 rhs = cpp_dec_float_50(1)/(2 * tol_sqrt * (1 - tol_sqrt*0.25) * acos(1 - tol_sqrt/2));
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
        double delta = sqrt(fourN - 3 * y * y);
	if (fmod(delta,1.0) == 0.0){
        double x1 = (y + delta) / 2;
        double x2 = (y - delta) / 2;
//	cout << "norm_eq: " << x1 << " " << x2 << " " << fmod(x1,1.0) << " " << fmod(x2,1.0) << endl;
        if (fmod(x1, 1.0) == 0.0) {
//	    cout << "success" << endl;
            Eiseinstein tmp(static_cast<long long>(x1), static_cast<long long>(y));
            results.push_back(tmp);
            Eiseinstein tmp2(static_cast<long long>(x2), static_cast<long long>(y));
            results.push_back(tmp2);
        }
	}
    }
    return results;
}

void solution(double theta, int k, double tol, vector<Eiseinstein>& matrix){
    bool matrix_found = false;

    cpp_dec_float_50 three = 3.0;
    cpp_dec_float_50 three_sqrt = sqrt(three);
    cpp_dec_float_50 radius_two = pow(three, ceil(k/2.0));
    cpp_dec_float_50 radius_one = radius_two * (1 - pow(tol,2.0)/2.0);

    double tol_sq = tol*tol;
    int f = ceil(k/2);
    int parity = k % 2;
    double three_to_f = pow(3.0, f);
    double three_to_k = pow(3.0, k);
    
    double three_sqrt_parity = pow(sqrt(3.0), -parity);
    vector<Eiseinstein> x1_vector = approximate_cnumber(-theta/2.0, k, tol, three, three_sqrt, radius_two, radius_one);
    vector<Eiseinstein> y2_vector = approximate_cnumber(theta/2.0, k, tol, three, three_sqrt, radius_two, radius_one);
    vector<Eiseinstein> z3_vector = approximate_cnumber(0, k, tol, three, three_sqrt, radius_two, radius_one);

    int x1_vector_size = x1_vector.size();
    int y2_vector_size = y2_vector.size();
    int z3_vector_size = z3_vector.size();
   
    Eiseinstein factor = Eiseinstein(2,1).pow(k);
    Eiseinstein prefactor = Eiseinstein(-1,-2).pow(k);
   
//    cout << x1_vector_size << " " << y2_vector_size<< " " << z3_vector_size << " " << endl;
/*    for (int i = 0; i < x1_vector_size; i++){Eiseinstein x1 = x1_vector[i]; cout << "x1 " << k << " " << x1.a << " " << x1.b << endl;}
    for (int i = 0; i < y2_vector_size; i++){Eiseinstein y2 = y2_vector[i]; cout << "y2 " << k << " " << y2.a << " " << y2.b << endl;}
    for (int i = 0; i < z3_vector_size; i++){Eiseinstein z3 = z3_vector[i]; cout << "z3 " << k << " " << z3.a << " " << z3.b << endl;}
*/    /*    for (int i = 0; i < x1_vector_size; i++){Eiseinstein x1 = x1_vector[i]; cout << "x1 " << k << " " << compute_error_sq_precompute(-theta/2, k, x1*prefactor,three_to_k) << endl;}
    for (int i = 0; i < y2_vector_size; i++){Eiseinstein y2 = y2_vector[i]; cout << "y2 " << k << " " << compute_error_sq_precompute(theta/2, k, y2*prefactor,three_to_k) << endl;}
    for (int i = 0; i < z3_vector_size; i++){Eiseinstein z3 = z3_vector[i]; cout << "z3 " << k << " " << compute_error_sq_zero_precompute(k, z3*prefactor,three_to_k) << endl;}
*/
    double radius_two_double = pow(sqrt(3.0), k);
    long long threeToK = pow(3.0, k);

    // We will try to save some computation of the norm equation
    map<long long, vector<Eiseinstein>> solution_normeq;

    vector<Eiseinstein> units = {Eiseinstein(1,0),Eiseinstein(-1,0),Eiseinstein(0,1),Eiseinstein(0,-1),Eiseinstein(1,1),Eiseinstein(-1,-1)};

    int sizex_sizey = x1_vector_size * y2_vector_size;

    //for (int i = 0; i < x1_vector_size; i++){
    #pragma omp parallel for
    for (int ij = 0; ij < sizex_sizey; ij ++){
        if (matrix_found) continue;



        // get i and j index from the overal ij index
        int i = ij / y2_vector_size;
        int j = ij - i * y2_vector_size;


        // deal with x1

        Eiseinstein x1 = x1_vector[i];
        Eiseinstein x1_conj = x1.conjugate();
        int x1_sde = x1.get_sde(k);
        // compute error and new tolerance
        double err1_sq = compute_error_sq_precompute(-theta/2, k, x1*prefactor,three_to_k);
        double tol2_sq = tol_sq - err1_sq;



        //compute additional quantities
        long long x1_norm = x1.norm();
        double x1_norm_sqrt = sqrt(static_cast<double>(x1_norm));
        long long eta1 = threeToK - x1_norm; // Used in the quadratic equation  eta |x2|^4 + (|beta|^2 - |alpha|^2 - eta^2)|x2|^2 + eta * |alpha|^2 = 0


        // (2) deal with y2


        //for (int j = 0; j < y2_vector_size; j++){


        if (matrix_found) continue;

        Eiseinstein y2 = y2_vector[j];
        Eiseinstein y2_conj = y2.conjugate();

        // compute error and new tolerance
        double err2_sq = compute_error_sq_precompute(theta/2, k, y2*prefactor,three_to_k);
        // break if we have already surpassed tol2


        // we need to handle this carefully. // if err2_sq surpasses tol2_sq ... move to the next x1 
        if (err2_sq > tol2_sq) {

            ij = (i + 1) * y2_vector_size - 1;
            continue;

        }

        if (err2_sq > tol2_sq) cout << "we should not see this" << endl;



        //Eiseinstein y2_conj = y2.conjugate();
        int y2_sde = y2.get_sde(k);

        if (x1_sde != y2_sde) continue;

        double tol3_sq = tol_sq - err1_sq - err2_sq;



        // helpful quantities
        long long y2_norm = y2.norm();
        double y2_norm_sqrt = sqrt(static_cast<double>(y2_norm));
        long long eta2 = threeToK - y2_norm;

        for (int kz = 0; kz < z3_vector_size; kz++){
            if (matrix_found) continue;

            Eiseinstein z3 = z3_vector[kz];
            // compute error and check against the tolerance
            //double err3_sq = compute_error_sq_precompute(0, k, z3*prefactor,three_to_k);
            double err3_sq = compute_error_sq_zero_precompute(k, z3*prefactor,three_to_k);
            if (err3_sq > tol3_sq) break;
            Eiseinstein z3_conj = z3.conjugate();
            int z3_sde = z3.get_sde(k);
            if (x1_sde != z3_sde) continue; 

            // helpful quantities
            long long z3_norm = z3.norm();
            long long eta3 = threeToK - z3_norm;
            double z3_norm_sqrt = sqrt(static_cast<double>(z3_norm));

            // Horn conditions: Necessary and sufficient conditions for x1, y2, z3 to
            // be diagonal entries of a unitary matrix
            if (x1_norm_sqrt + y2_norm_sqrt - z3_norm_sqrt > radius_two_double) continue;
            if (x1_norm_sqrt + z3_norm_sqrt - y2_norm_sqrt > radius_two_double) continue;
            if (y2_norm_sqrt + z3_norm_sqrt - x1_norm_sqrt > radius_two_double) continue;

	        Eiseinstein x1y2_conj = x1_conj * y2_conj;
	        Eiseinstein z3_factor = factor * z3;
            for (int ia=0; ia < 6; ia++){
                if (matrix_found) continue;
                Eiseinstein alpha = x1y2_conj - units[ia] * z3_factor;
                long long alpha_norm = alpha.norm();
                //if (alpha_norm > eta1 * eta2 || alpha_norm > eta1*eta1 || alpha_norm > eta2*eta2) continue;
                if (alpha_norm > eta1*eta1 || alpha_norm > eta2*eta2) continue;

                Eiseinstein beta = x1_conj * z3_conj - units[ia] * factor * y2;
                long long beta_norm = beta.norm();
                //if (beta_norm > eta1*eta3 || beta_norm > eta1*eta1 || beta_norm > eta3*eta3) continue;
                if (beta_norm > eta1*eta1 || beta_norm > eta3*eta3) continue;

                // Get also gamma
                Eiseinstein gamma = y2_conj * z3_conj - units[ia] * factor * x1;
                long long gamma_norm = gamma.norm();
                //if (gamma_norm > eta2 * eta3 || gamma_norm > eta2*eta2 || gamma_norm > eta3*eta3) continue;
                if (gamma_norm > eta2*eta2 || gamma_norm > eta3*eta3) continue;

                // Build the quadratic equation
                long long a1 = eta1;
                long long a2 = beta_norm - alpha_norm - eta1*eta1;
                long long a3 = eta1 * alpha_norm;
                double discr = a2*a2 - 4 * a1 * a3;

                if (discr < 0) continue;
                double val1 = (-a2 + sqrt(discr))/(2*a1);
                double val2 = (-a2 - sqrt(discr))/(2*a1);
                vector<long long> x2_norms;
	            double eta_min = min(eta1,eta2);
                if (floor(val1) == val1 && val1 > 0 && val1 <= eta_min) x2_norms.push_back(val1);
                if (floor(val2) == val2 && val2 > 0 && val2 <= eta_min) x2_norms.push_back(val2);
                for (int ix2_n=0; ix2_n<x2_norms.size();ix2_n++){
                    if (matrix_found) continue;
                    long long x2_norm = x2_norms[ix2_n];
                    long long x3_norm = eta1 - x2_norm;
                    long long z2_norm = eta2 - x2_norm;
                    long long z1_norm = eta3 - z2_norm;
                    long long y1_norm = eta1 - z1_norm;
                    long long y3_norm = eta2 - y1_norm;

                    
                    // Trying to save results of solving the norm equations
                    vector<Eiseinstein> x2_vector;
                    vector<Eiseinstein> x3_vector;
                    auto x2_iter = solution_normeq.find(x2_norm);
                    if (x2_iter != solution_normeq.end()){
                        x2_vector = x2_iter->second;

                    } 
                    else
                    {
                        x2_vector = solve_norm_eq(x2_norm);
                        #pragma omp critical
                        {
                            solution_normeq[x2_norm] = x2_vector;
                        }
                        
                    }

                    auto x3_iter = solution_normeq.find(x3_norm);

                    if (x3_iter != solution_normeq.end()){
                        x3_vector = x3_iter->second;

                    } 
                    else
                    {

                        x3_vector = solve_norm_eq(x3_norm);
                        #pragma omp critical
                        {
                            solution_normeq[x3_norm] = x3_vector;
                        }
                        
                    }
                    // End of method of saving results

                    // Start of the method of solving the norm equation without saving results
                    //vector<Eiseinstein> x2_vector = solve_norm_eq(x2_norm);
                    //vector<Eiseinstein> x3_vector = solve_norm_eq(x3_norm);



//			cout << "x2 x3: " << x2_vector.size() << " " << x3_vector.size() << endl;
//    for (int i = 0; i < x2_vector.size(); i++){Eiseinstein x2 = x2_vector[i]; cout << "x2 " << k << " " << x2.a << " " << x2.b << endl;}
//    for (int i = 0; i < x3_vector.size(); i++){Eiseinstein x3 = x3_vector[i]; cout << "x3 " << k << " " << x3.a << " " << x3.b << endl;}
                    
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

                            double final_error = sqrt(err1_sq + err2_sq + err3_sq);
                            #pragma omp critical
                            {

                                x1.print();
                                x2.print();
                                x3.print();
                                y1.print();
                                y2.print();
                                y3.print();
                                z1.print();
                                z2.print();
                                z3.print();

                                matrix = {x1,x2,x3,y1,y2,y3,z1,z2,z3};
			                     cout << "errors: " << err1_sq << " " << err2_sq << " " << err3_sq << endl;
                                cout << "answer: " << x1.a << " " << x1.b << " " << y2.a << " " << y2.b << " " << z3.a << " " << z3.b << endl;
                                cout << theta << "," << k << "," << tol << "," << final_error << endl;
                            }
                        }
                    }
                }
            }

        }

    }

}





