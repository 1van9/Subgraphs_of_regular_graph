#include "matrix.h"
#include "Zp.h"
#include "rational.h"
#include "poly.h"
#include "isomorphism.h"
#include "clique-addition.h"

using Q = Frac<long long>;
using Z3 = Zp<3>;
using Rf = Frac<Poly<Q>>;
using MatRf = Matrix<Rf>;
using MatQ = Matrix<Q>;
using MatN = Matrix<long long>;
using Mat3 = Matrix<Z3>;
using CliqueRf = vec<MatRf>;
using CliqueQ = vec<MatQ>;
using Clique3 = vec<Mat3>;


Rf x = Rf(vec<Q>({0, 1}));

MatQ get(const MatRf & a, const Q & v) {
    MatQ b(a.n, a.m);
    for (int i = 0; i < a.n; i++) {
        for (int j = 0; j < a.m; j++) {
            b[i][j] = a[i][j].fi(v) / a[i][j].se(v);
        }
    }
    return b;
}

CliqueQ get(const CliqueRf & a, const Q & v) {
    CliqueQ b(a.size());
    for (int i = 0; i < a.size(); i++) {
        b[i] = get(a[i], v);
    }
    return b;
}

template<typename T>
void print(vec<Matrix<T>> & a) {
    cout << "print matrix : " << endl; 
    for (auto e : a)
        cout << e << endl;
}

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

template<typename T>
void generate_small_matrices(int i, int j, int N, Matrix<T> & curr_matrix, const vec<T>& numbers, vec<Matrix<T>> & matrices) {
    if (j == N)   //shift
        i++, j = 0;
    if (i == N) {
        if (curr_matrix.GL())
            matrices.push_back(curr_matrix);
        return;
    }
    for (auto number : numbers) {
        curr_matrix[i][j] = number;
        generate_small_matrices(i, j + 1, N, curr_matrix, numbers, matrices);
    }
}


template<typename T>
vec<Matrix<T>> get_matrixes(const vec<T> & nums, int N) {
    Matrix<T> curr(N, N);
    vec<Matrix<T>> allm;
    generate_small_matrices(0, 0, N, curr, nums, allm);
    return allm;
}
