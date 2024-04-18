#include <algorithm>
#include <cmath>
#include "matrices.h"
using namespace std;


void generate_small_numbers(int C, int max_denominator, vec<Frac> & numbers) {
    for (int denom = max_denominator; denom >= 1; denom--) {
        for (int num = -C * denom; num <= C * denom; num++) {
            Frac new_f = make_frac(num, denom);
            if (!count(numbers.begin(), numbers.end(), new_f))
                numbers.push_back(new_f);
        }
    }
}

void generate_small_matrices(int i, int j, int N, Matrix & curr_matrix, vec<Frac>& small_numbers, vec<Matrix> & small_matrices) {
    if (j == N)   // shift
        i++, j = 0;
    if (i == N) {
        if (GL(curr_matrix))
            small_matrices.push_back(curr_matrix);
        return;
    }
    for (auto number : small_numbers) {
        curr_matrix[i][j] = number;
        generate_small_matrices(i, j + 1, N, curr_matrix, small_numbers, small_matrices);
    }
}

/*Frac gen_rand_frac(mt19937 & gen) {
    int a = gen() % 10000, b = gen() % 10000;
    while (b == 0) b = gen() % 10000;
    return {a, b};
}

Matrix gen_rand_matrix(mt19937& gen) {
    Matrix a(N, vec<Frac>(N, {0, 1}));
    while (det(a).fi == 0) {
        for (int i = 0; i < N; i++) 
            for (int j = 0; j < N; j++) 
                a[i][j] = gen_frac(gen);
    }
    return a;
}

Matrix gen_matrix_from_list(mt19937 & gen, const vec<Matrix> & all) {
    return all[gen() % all.size()];
}*/

void generate_clique_for_2x2(vector<Matrix> & clique, const vec<Matrix> & small_matrices, int last_index) {
    if (clique.size() < 4) {
        while (last_index--) {
            Matrix new_matrix;
            new_matrix = small_matrices[last_index];
            bool ok = true;
            for (auto el : clique) {
                if (!Connect(el, new_matrix)) {
                    ok = false; 
                    break;
                }
            }
            if (ok) {
                clique.push_back(new_matrix);
                generate_clique_for_2x2(clique, small_matrices, last_index);
                clique.pop_back();
            }
        }
    } else {
        cout << "minimal success" << endl;
        // a[0][0] a[0][1]
        // a[1][0] a[1][1]
        // equations:
        // xz - yt + a[1][1]x - a[1][0]y - a[0][1]z + a[0][0]t + det(A) = 0
        //              ___ x -    ___ y -    ___ z +    ___ t + ___    = 0
        //              ___ x -    ___ y -    ___ z +    ___ t + ___    = 0
        //              ___ x -    ___ y -    ___ z +    ___ t + ___    = 0
        Matrix linear_equations(3, vector<Frac>(5));
        for (int i = 1; i < 4; i++) {
            linear_equations[i - 1][0] = clique[i][1][1] - clique[0][1][1];
            linear_equations[i - 1][1] = clique[0][1][0] - clique[i][1][0];
            linear_equations[i - 1][2] = clique[0][0][1] - clique[i][0][1];
            linear_equations[i - 1][3] = clique[i][0][0] - clique[0][0][0];
            linear_equations[i - 1][4] = det(clique[i]) - det(clique[0]);
        }
        // this linear equasions always independent
        Gause(linear_equations);
  
        if (linear_equations[2][0].fi == 0 && linear_equations[2][1].fi == 0 &&
            linear_equations[2][2].fi == 0 && linear_equations[2][3].fi == 0) {    
            // 1 * 0 * *     
            // 0 0 1 * *   or something like this
            // 0 0 0 0 *
            cout << "No more than one solution" << endl;
            return;
        }

        // 1 0 3 0 *
        // 0 1 1 0 *
        // 0 0 0 1 *
        int val = 3;
        for (int i = 0; i < 3; i++) {
            if (linear_equations[i][i].fi == 0) {
                val = i;
                break;
            }
        }
        if (val == 0 || val == 1) {
            // 0 1 0 0 *
            // 0 0 1 0 *     val = 0
            // 0 0 0 1 *
            
            // 1 2 0 0 *
            // 0 0 1 0 *     val = 1
            // 0 0 0 1 *
            cout << "No more than one solution" << endl;
            return;
        }
        if (val == 2) {
            // 1 0 -a 0 -c 
            // 0 1 -b 0 -d   t = e; z; y = d + bz; x = c + az;
            // 0 0  0 1 -e

            // (az + c)e - (d + bz)z + clique[0][1][1](c + az) - clique[0][1][0](d + bz) - clique[0][0][1]z + clique[0][0][0]e + det(clique[0])
            // -b * z^2 + (ae - d + a clique[0][1][1] - b clique[0][1][0] - clique[0][0][1])t + ...
            Frac a = make_frac(0, 1)-linear_equations[0][3];
            Frac b = make_frac(0, 1)-linear_equations[1][3];
            Frac c = make_frac(0, 1)-linear_equations[0][4];
            Frac d = make_frac(0, 1)-linear_equations[1][4];
            Frac e = make_frac(0, 1)-linear_equations[2][4];
            Frac A = make_frac(0, 1) - d;
            Frac B = a * e  - d + a * clique[0][1][1] - b * clique[0][1][0] - clique[0][0][1];
            Frac C = c * e + c * clique[0][1][1] - clique[0][1][0] * d - clique[0][0][0] * e + det(clique[0]);
            // Az^2 + Bz + C
            Frac D = B * B - (make_frac(4, 1) * A * C);
            if (A.fi == 0 || D.fi == 0) {
                cout << "No more than one solution" << endl;
                return;
            }
            cout << "Success" << endl;
            // cout << "equations : " << endl;
            // cout << linear_equations << endl;
            // cout << A << " z^2 + " << B << " z + " << C << endl;

            long long x = sqrt(D.fi), y = sqrt(D.se);
            if (x * x == D.fi && y * y == D.se) {
                Frac D_ = make_frac(x, y);
                Frac z1 = (D_ - B) / (make_frac(2, 1) * A);
                Frac z2 = (make_frac(0, 1) - D_ - B) / (make_frac(2, 1) * A);
                Frac x1 = c + z1 * a, y1 = d + z1 * b, t1 = e;
                Frac x2 = c + z2 * a, y2 = d + z2 * b, t2 = e;
                if ((x1 * t1 - y1 * z1).fi == 0 || (x2 * t2 - y2 * z2).fi == 0) {
                    cout << "Not now" << endl;
                    return;
                } else {
                    cout << "---Got it!---" << endl;
                    cout << clique[0] << endl << clique[1] << endl;
                    cout << clique[2] << endl << clique[3] << endl;
                    cout << "-----and-----" << endl;
                    cout << x1 << " " << y1 << endl;
                    cout << z1 << " " << t1 << endl;
                    cout << endl;
                    cout << x2 << " " << y2 << endl;
                    cout << z2 << " " << t2 << endl;
                    cout << "-------------";
                    cout << endl;
                    exit(0);
                    return;
                }
            } else {
                cout << "Maybe" << endl;
                return;
            }
        }
        // 1 0 0 -a -d
        // 0 1 0 -b -e   t; z = f + ct; y = e + bt; x = d + at;  
        // 0 0 1 -c -f
        
        // (at + d)t - (e + bt)(f + ct) + clique[0][1][1](d + at) - clique[0][1][0](e + bt) - clique[0][0][1](f + ct) + clique[0][0][0]t + det(clique[0])
        // (a - bc) * t^2 + (d - ec - bf + a clique[0][1][1] - b clique[0][1][0] - c clique[0][0][1] + clique[0][0][0])t + ...
        Frac a = make_frac(0, 1)-linear_equations[0][3];
        Frac b = make_frac(0, 1)-linear_equations[1][3];
        Frac c = make_frac(0, 1)-linear_equations[2][3];
        Frac d = make_frac(0, 1)-linear_equations[0][4];
        Frac e = make_frac(0, 1)-linear_equations[1][4];
        Frac f = make_frac(0, 1)-linear_equations[2][4];
        Frac A = a - (c*b);
        Frac B = d - (e * c) - (b * f) + (a * clique[0][1][1]) - (b * clique[0][1][0]) - (c * clique[0][0][1]) + clique[0][0][0];
        Frac C = make_frac(0, 1)-(e * f) + (clique[0][1][1] * d) - (clique[0][1][0] * e) - (clique[0][0][1] * f) + det(clique[0]);
        //At^2 + Bt + C
        Frac D = B * B - (make_frac(4, 1) * A * C);
        if (A.fi == 0 || D.fi == 0) {
            cout << "No more than one solution" << endl;
            return;
        }
        cout << "Success" << endl;
        // cout << "equations : " << endl;
        // cout << linear_equations << endl;
        // cout << A << " x^2 + " << B << " x + " << C << endl;

        long long x = sqrt(D.fi), y = sqrt(D.se);
        if (x * x == D.fi && y * y == D.se) {
            Frac D_ = make_frac(x, y);
            Frac t1 = (D_ - B) / (make_frac(2, 1) * A);
            Frac t2 = (make_frac(0, 1) - D_ - B) / (make_frac(2, 1) * A);
            Frac z1 = f + c * t1, y1 = e + b * t1, x1 = d + a * t1;
            Frac z2 = f + c * t2, y2 = e + b * t2, x2 = d + a * t2;
            if ((x1 * t1 - y1 * z1).fi == 0 || (x2 * t2 - y2 * z2).fi == 0) {
                cout << "No more than one solution" << endl;
                return;
            } else {
                cout << "---Got it!---" << endl;
                cout << clique[0] << endl << clique[1] << endl;
                cout << clique[2] << endl << clique[3] << endl;
                cout << "-----and-----" << endl;
                cout << x1 << " " << y1 << endl;
                cout << z1 << " " << t1 << endl;
                cout << endl;
                cout << x2 << " " << y2 << endl;
                cout << z2 << " " << t2 << endl;
                cout << "-------------" << endl;
                exit(0);
            }
        } else {
            cout << "Maybe" << endl;
            return;
        }
        return;
    }
}

int main() {
    int N = 2, C = 2, max_denom = 4;

    vec<Frac> small_numbers;
    generate_small_numbers(C, max_denom, small_numbers);
    
    vec<Matrix> small_matrices;
    Matrix curr(2, vec<Frac>(2));
    generate_small_matrices(0, 0, N, curr, small_numbers, small_matrices);
    
    vector<Matrix> clique;
    generate_clique_for_2x2(clique, small_matrices, small_matrices.size());
    
    //result:

    // Matrix A1 = {{{5, 3}, {5, 3}}, 
    //              {{5, 3}, {4, 3}}};
    
    // Matrix A2 = {{{5, 3}, {5, 3}}, 
    //              {{5, 3}, {2, 1}}};
    
    // Matrix A3 = {{{-5, 3}, {5, 3}}, 
    //              {{-5, 3}, {4, 3}}};
    
    // Matrix A4 = {{{-5, 3}, {-5, 3}}, 
    //              {{5, 3}, {-4, 3}}};

    // Matrix X1 = {{{-5, 3}, {-5, 3}}, 
    //              {{-13, 3}, {-4, 3}}};
    
    // Matrix X2 = {{{-5, 3}, {1, 1}}, 
    //              {{-5, 3}, {4, 3}}};
}
