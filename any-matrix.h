#include "matrix.h"
#include "rational.h"
#include "poly.h"
#include "isomorphism.h"
#include "clique-addition.h"

using Q = Frac<long long>;
using Rf = Frac<Poly<Q>>;
using MatRf = Matrix<Rf>;
using MatQ = Matrix<Q>;
using MatN = Matrix<long long>;
using CliqueRf = vec<MatRf>;
using CliqueQ = vec<MatQ>;

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
