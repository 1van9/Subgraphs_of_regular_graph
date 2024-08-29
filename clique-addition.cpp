#include <algorithm>
#include <random>
#include "any-matrix.h"
#include <bitset>

using namespace std;

// void generate_small_numbers(int C, int max_denominator, vec<Frac> & numbers) {
//     for (int denom = max_denominator; denom >= 1; denom--) {
//         for (int num = -C * denom; num <= C * denom; num++) {
//             Frac new_f = make_frac(num, denom);
//             if (!count(numbers.begin(), numbers.end(), new_f))
//                 numbers.push_back(new_f);
//         }
//     }
// }

void generate_small_matrices(int i, int j, int N, Matrix<Poly<int>> & curr_matrix, vec<Poly<int>>& small_numbers, vec<Matrix<Poly<int>>> & small_matrices) {
    if (j == N)   //shift
        i++, j = 0;
    if (i == N) {
        if (curr_matrix.GL())
            small_matrices.push_back(curr_matrix);
        return;
    }
    for (auto number : small_numbers) {
        curr_matrix[i][j] = number;
        generate_small_matrices(i, j + 1, N, curr_matrix, small_numbers, small_matrices);
    }
}

// inline Frac gen_rand_frac(mt19937 & gen) {
//     int a = gen() % 10000, b = gen() % 10000;
//     while (b == 0) b = gen() % 10000;
//     return {a, b};
// }

// inline Matrix gen_rand_matrix(int N, mt19937& gen) {
//     Matrix a(N, vec<Frac>(N, {0, 1}));
//     while (det(a).fi == 0) {
//         for (int i = 0; i < N; i++) 
//             for (int j = 0; j < N; j++) 
//                 a[i][j] = gen_rand_frac(gen);
//     }
//     return a;
// }

inline MatRf gen_matrix_from_list(mt19937 & gen, const vec<MatRf> & all) {
    return all[gen() % all.size()];
}

vector<bool> was_add(6);
mt19937 gen(4385);
int cnt = 3;
const int sm_mat_sz = 2112;

Rf to_rf(const Poly<int> & x) {
    vec<Q> coefs;
    for (auto a : x.coefficients)
        coefs.emplace_back(a);
    Rf y(Poly<Q>(coefs), Poly<Q>(Q(1)));
    return y;
}

MatRf to_rf(const Matrix<Poly<int>> & x) {
    MatRf y(x.n, x.m);
    for (int i = 0; i < x.n; i++)
        for (int j = 0; j < x.m; j++)
            y[i][j] = to_rf(x[i][j]);
    return y;
}

void generate_clique_for_2x2(vec<int> & clique, const vec<Matrix<Poly<int>>> & small_matrices, int last_index, 
                             const vec<bitset<sm_mat_sz>> & g, const bitset<sm_mat_sz> & neighbours) {
    if (!cnt)
        return;
    if (clique.size() < 4) {
        while (last_index--) {
            if (neighbours[last_index]) {
                clique.push_back(last_index);
                generate_clique_for_2x2(clique, small_matrices, last_index, g, neighbours & g[last_index]);
                clique.pop_back();
            }
        }
    } else if (clique.size() == 4) {
        vec<MatRf> cl;
        for (auto i : clique)
            cl.push_back(to_rf(small_matrices[i]));
        
        vec<MatRf> add = is_additionable(cl);
        for (auto el : add) {
            cl.push_back(el);
            int add = addition(cl);
            cout << "find! addition = " << add << endl;
            if (add == 1 && cnt--) {
                for (auto el : cl) {
                    cout << el << endl;
                }
                cout << "---------------" << endl;
            }
            if (add == 3 && cnt--) {
                for (auto el : cl) {
                    cout << el << endl;
                }
                cout << "---------------" << endl;
            }
            if (add == 5 && cnt--) {
                for (auto el : cl) {
                    cout << el << endl;
                }
                cout << "---------------" << endl;
            }
            cl.pop_back();
        }
    }
    return;
}



int main() {
    int N = 2;
    Poly<int> y({0, 1});
    vec<Poly<int>> small_numbers = {Poly<int>(0), Poly<int>(1), -Poly<int>(1), y, -y, y + y - Poly<int>(1), -y-y + Poly<int>(1)};
    // generate_small_numbers(3, 2, small_numbers);
    Matrix<Poly<int>> curr(2);
    vec<Matrix<Poly<int>>> small_matrices;
    generate_small_matrices(0, 0, N, curr, small_numbers, small_matrices);

    vec<bitset<sm_mat_sz>> g(sm_mat_sz);
    int edges = 0;
    for (int i = 0; i < sm_mat_sz; i++) {
        for (int j = 0; j < sm_mat_sz; j++)
            if (small_matrices[i].Connect(small_matrices[j]))
                g[i][j] = 1, edges++;
        if (i % 100 == 0) {
            cout << "o";
        }
    }
    cout << "k" << endl;
    vec<int> clique;
    bitset<sm_mat_sz> nei;
    for (int i = 0; i < sm_mat_sz; i++)
        nei[i] = 1;
    generate_clique_for_2x2(clique, small_matrices, small_matrices.size(), g, nei);
}
