
#include "enumeration.h"


vector<pair<Eiseinstein, Eiseinstein>> enumerate(mpreal alpha, mpreal r1, mpreal r2, Eiseinstein factor){
    
    // return vector of (x1, x2) such that |x1|^2 + |x2|^2 <= r2^2 and such that Re[u1_conjugate x1 + u2_conj * x2] >= r1

    vector<pair<Eiseinstein, Eiseinstein>> results;
    mpreal r1_squared = r1 * r1;
    mpreal r2_squared = r2 * r2;
    mpreal cosalpha = cos(alpha);
    mpreal sinalpha = sin(alpha);
    mpreal three_sqrt = mpfr::sqrt(3.0);
    mpreal two_sqrt = mpfr::sqrt(2.0);
    mpreal half = 0.5;
    mpreal q2_min, q2_max, q2, p2_min, p2_max, p2, y1_min, y1_max, p1_min, p1_max, p1, y2_min, y2_max, q1_min, q1_max, q1;


    // sample q2
    mpreal T4_squared = r2_squared - r1_squared;
    mpreal T4 = sqrt(T4_squared);
    q2_min = ceil(-2.0 * T4/three_sqrt) / 2.0;
    q2_max = floor(2.0 * T4/three_sqrt) / 2.0;

    for (q2 = q2_min; q2 <= q2_max; q2 += half){
        mpreal y4 = three_sqrt * q2;
        mpreal y4_squared = y4*y4;

        mpreal T3_squared = r2_squared - r1_squared - y4_squared;
        mpreal T3 = sqrt(T3_squared);

        p2_min = ceil(-two_sqrt * (T3 + r1)) / 2.0;
        p2_max = floor(two_sqrt * (T3 - r1)) / 2.0;

        if (floor(p2_min + q2) != p2_min + q2) p2_min += half;

        for (p2 = p2_min; p2 <= p2_max; p2++){
            mpreal y3 = p2; // not needed but it helps me think
            mpreal y3_squared = y3*y3;

            // we have x2
            Eiseinstein x2_prime(p2 + q2, 2*q2);
            auto x2_data = x2_prime / factor;
            if (x2_data.second.norm() != 0) continue;
            Eiseinstein x2 = x2_data.first;


            // sample the p1 and q1
            mpreal lambda1 = two_sqrt * r1 + y3;
            mpreal lambda1_squared = lambda1 * lambda1;
            mpreal lambda2_squared = r2_squared - y4_squared - y3_squared;
            mpreal lambda2 = sqrt(lambda2_squared);

            mpreal delta = acos(lambda1/lambda2);
            mpreal phi1 = alpha - delta;
            mpreal phi2 = alpha + delta;
            // compute lambda2 * cos(phi1) and lambda2 * cos(phi2)
            // the smallest value should be assigned to y1_min and the largest to y1_max
            mpreal tmp_value1 = lambda2 * cos(phi1);
            mpreal tmp_value2 = lambda2 * cos(phi2);
            tie(y1_min, y1_max) = minmax(tmp_value1, tmp_value2);
            if (phi1 * phi2 < 0) y1_max = lambda2;
            p1_min = ceil(2 * y1_min) / 2.0;
            p1_max = floor(2 * y1_max) / 2.0;

            for (p1 = p1_min; p1 <= p1_max; p1 += half){
                mpreal y1 = p1;
                mpreal y1_squared = y1*y1;

                
                mpreal T2 = sqrt(lambda2_squared - y1_squared);
                mpreal y2_min = -T2;
                mpreal y2_max = T2;
                if (sinalpha > 0) y2_min = (lambda1 - y1 * cosalpha) / sinalpha;
                if (sinalpha < 0) y2_max = (lambda1 - y1 * cosalpha) / sinalpha;

                q1_min = ceil(2 * y2_min/three_sqrt) / 2.0;
                q1_max = floor(2 * y2_max/three_sqrt) / 2.0;
                if (floor(q1_min + p1) != q1_min + p1) q1_min += half;

                for (q1 = q1_min; q1 <= q1_max; q1++){

                    Eiseinstein x1_prime(p1 + q1, 2*q1);
                    if (x1_prime.norm() > lambda2_squared) continue;
                    auto x1_data = x1_prime / factor;
                    if (x1_data.second.norm() != 0) continue;
                    Eiseinstein x1 = x1_data.first;

                    results.push_back(make_pair(x1, x2));

                }

            }

        }
    }
    return results;
}


// vector<Eiseinstein> x1_candidates = approx_cnumber(alpha, lambda1, lambda2);
// int x1_cand_size = x1_candidates.size();
// for (int ii = 0; ii < x1_cand_size; ii++){

//     // now we have x1_prime
//     Eiseinstein x1_prime = x1_candidates[ii];

//     auto x1_data = x1_prime / factor;
//     if (x1_data.second.norm() != 0) continue;
//     Eiseinstein x1 = x1_data.first;

//     results.push_back(make_pair(x1, x2));

// }  
    

