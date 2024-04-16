#include <iostream>
#include <random>
#include <vector>
#include <cassert>
#include <algorithm>
#include <set>
using namespace std;

namespace {
    
#define vec vector
#define fi first
#define se second
#define mfrac make_pair
    
using Frac = pair<long long, long long>;
using Matrix = vec<vec<Frac>>;

long long gcd(long long a, long long b) {
    if (b == 0) return a;
    return gcd(b, a % b);
}

Frac upd(const Frac & a) {
    if (a.fi == 0 && a.se == 0) return a;
    if (a.se < 0)
        return {-a.fi / gcd(abs(a.fi), -a.se), -a.se / gcd(abs(a.fi), -a.se)};
    else
        return {a.fi / gcd(abs(a.fi), a.se), a.se / gcd(abs(a.fi), a.se)};
}
Frac operator+ (const Frac & a, const Frac & b) {
    return upd({a.fi * b.se + b.fi * a.se, a.se * b.se});
}
Frac operator+= (Frac & a, const Frac & b) {
    return a = upd({a.fi * b.se + b.fi * a.se, a.se * b.se});
}
Frac operator- (const Frac & a, const Frac & b) {
    return upd({a.fi * b.se - b.fi * a.se, a.se * b.se});
}
Frac operator-= (Frac & a, const Frac & b) {
    return a = upd({a.fi * b.se - b.fi * a.se, a.se * b.se});
}
Frac operator* (const Frac & a, const Frac & b) {
    return upd({a.fi * b.fi, a.se * b.se});
}
Frac operator*= (Frac & a, const Frac & b) {
    return a = upd({a.fi * b.fi, a.se * b.se});
}
Frac operator/ (const Frac & a, const Frac & b) {
    return upd({a.fi * b.se, a.se * b.fi});
}
Frac operator/= (Frac & a, const Frac & b) {
    return a = upd({a.fi * b.se, a.se * b.fi});
}
bool operator!=(const Frac & a, const Frac & b) {
    return a.fi * b.se != a.se * b.fi;
}
bool operator==(const Frac & a, const Frac & b) {
    return a.fi * b.se == a.se * b.fi;
}
Matrix operator+ (const Matrix & a, const Matrix & b) {
    return {{a[0][0] + b[0][0], a[0][1] + b[0][1]}, {a[1][0] + b[1][0], a[1][1] + b[1][1]}};
}
Matrix operator- (const Matrix & a, const Matrix & b) {
    return {{a[0][0] - b[0][0], a[0][1] - b[0][1]}, {a[1][0] - b[1][0], a[1][1] - b[1][1]}};
}
Matrix operator* (const Matrix & a, const Matrix & b) {
    return {{a[0][0] * b[0][0] + a[0][1] * b[1][0], a[0][0] * b[0][1] + a[0][1] * b[1][1]}, 
            {a[1][0] * b[0][0] + a[1][1] * b[1][0], a[1][0] * b[0][1] + a[1][1] * b[1][1]}};
}
Matrix T(const Matrix & a) {
    return {{a[0][0], a[1][0]}, {a[0][1], a[1][1]}};
}
Frac det(const Matrix & a) {
    return a[0][0] * a[1][1] - a[0][1] * a[1][0];
}
bool Connect(const Matrix & a, const Matrix & b) {
    return det(a + b).fi == 0;
}
bool GL(const Matrix & a) {
    return det(a).fi != 0;
}
ostream& operator<< (ostream& out, const Frac & a) {
    if (a.fi == 0) {
        out << 0;
        return out;
    }
    if (a.se == 1) {
        out << a.fi;
        return out;
    }
    if (a.se < 0) 
        out << -a.fi << "/" << -a.se;
    else
        out << a.fi << "/" << a.se;
    return out;
}
ostream& operator<< (ostream& out, const Matrix & a) {
    if (a.size() == 2 && a[0].size() == 2)
        out << a[0][0] << " " << a[0][1] << endl << a[1][0] << " " << a[1][1] << endl;
    else {
        for (int i = 0; i < a.size(); i++) {
            for (auto el : a[i]) {
                cout << el << " ";
            }
            cout << endl;
        }
    }
    return out; 
}

void Gause(Matrix& A) {
    int n = A.size(), m = A[0].size();

    int curr_string = 0;
    for (int i = 0; i < m && curr_string < n; i++) {
        int curr = curr_string;
        while (curr < n && A[curr][i].fi == 0)
            curr++;
        if (curr == n) continue;
        if (curr != curr_string) {
            for (int j = 0; j < m; j++) {
                swap(A[curr][j], A[curr_string][j]);
            }
        }
        Frac div = A[curr_string][i];
        for (int j = 0; j < m; j++) {
            A[curr_string][j] /= div;
        }
        for (int j = 0; j < n; j++) {
            if (j == curr_string)
                continue;
            Frac mul = A[j][i];
            for (int k = 0; k < m; k++) {
                A[j][k] -= A[curr_string][k] * mul;
            }
        }
        curr_string++;
    }
}

Matrix rev(const Matrix& A) {
    int n = A.size();
    Matrix B(n, vec<Frac> (2 * n, {0, 1}));
    for (int i = 0; i < n; i++) {
        B[i][i + n] = {1, 1};
        for (int j = 0; j < n; j++) {
            B[i][j] = A[i][j];
        }
    }
    Gause(B);
    bool check = true;
    for (int i = 0; i < n; i++) {
        if (B[i][i].fi != B[i][i].se) {
            check = false;
        }
    }
    if (!check) {
        assert(0);
    }
    Matrix res(n, vec<Frac>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            res[i][j] = B[i][j + n];
        }
    }
    return res;
}
};

bool isomorphic_subgraphs_2_by_2(const vector<Matrix>& subg1, const vector<Matrix>& subg2) {
    if (subg1.size() != subg2.size())
        return false;
    int n = subg1.size();
    vector<int> p(n);
    for (int i = 0; i < n; i++)
        p[i] = i;
    bool ok = false;
    // A_i = U B_{p_i} V
    // U = B_{p_i} V (A_i)^-1   (1)

    // V = (x y
    //      z t)
    // Have n wasy, how to express U, in term of x, y, z, t
    set<Matrix> ways; // possible linear equasions on coefficients x, y, z, t
    do {
        vector<Matrix> M; // express U coefficients in term of x, y, z, t 
        for (int i = 0; i < n; i++) {
            Matrix A = subg2[p[i]], B = rev(subg1[i]); 
            Matrix U_coefficients  = {
                {A[0][0] * B[0][0], A[0][0] * B[1][0], A[0][1] * B[0][0], A[0][1] * B[1][0]},
                {A[1][0] * B[0][0], A[1][0] * B[1][0], A[1][1] * B[0][0], A[1][1] * B[1][0]},
                {A[0][0] * B[0][1], A[0][0] * B[1][1], A[0][1] * B[0][1], A[0][1] * B[1][1]},
                {A[1][0] * B[0][1], A[1][0] * B[1][1], A[1][1] * B[0][1], A[1][1] * B[1][1]}
            }; 
            M.push_back(U_coefficients);
        }
        Matrix linear_equasions((n - 1) * M[0].size(), vec<Frac>(4)); // equasions on x, y, z, t
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < M[i].size(); j++) {
                for (int k = 0; k < M[i][j].size(); k++) {
                    // from U_0 - U_i = 0, where U_i = B_{p_i} V (A_i)^-1, we have:
                    linear_equasions[(i - 1) * (M[i].size()) + j][k] = M[0][j][k] - M[i][j][k]; }
            }
        }
        Gause(linear_equasions);
        ways.insert(linear_equasions);
    } while (next_permutation(p.begin(), p.end()));

    // tha same for tansposition
    do {
        vector<Matrix> M;
        for (int i = 0; i < n; i++) {
            Matrix A = subg2[p[i]], B = rev(T(subg1[i]));
            Matrix U_coefficients = {
                {A[0][0] * B[0][0], A[0][0] * B[1][0], A[0][1] * B[0][0], A[0][1] * B[1][0]},
                {A[1][0] * B[0][0], A[1][0] * B[1][0], A[1][1] * B[0][0], A[1][1] * B[1][0]},
                {A[0][0] * B[0][1], A[0][0] * B[1][1], A[0][1] * B[0][1], A[0][1] * B[1][1]},
                {A[1][0] * B[0][1], A[1][0] * B[1][1], A[1][1] * B[0][1], A[1][1] * B[1][1]}
            }; 
            M.push_back(U_coefficients);
        }
        Matrix linear_equasions((n - 1) * M[0].size(), vec<Frac>(4));
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < M[i].size(); j++) {
                for (int k = 0; k < M[i][j].size(); k++) {
                    linear_equasions[(i - 1) * (M[i].size()) + j][k] = M[0][j][k] - M[i][j][k];
                }
            }
        }
        Gause(linear_equasions);
        ways.insert(linear_equasions);
    } while (next_permutation(p.begin(), p.end()));

    // check is there non zero solutions of x, y, z, t, which give V from GL_2

    for (auto lin_eq : ways) {
        if (lin_eq[3][3].fi == lin_eq[3][3].se) {
            // 1 0 0 0
            // 0 1 0 0 
            // 0 0 1 0
            // 0 0 0 1
            continue;
        }
        if (lin_eq[2][2].fi == 0 && lin_eq[2][3].fi == lin_eq[2][3].se) {
            //     0 0           0
            // 0 0 1 0     0 1 * 0
            // 0 0 0 1     0 0 0 1
            if (lin_eq[1][2].fi == lin_eq[1][2].se && lin_eq[1][1].fi == 0)
                continue;
            if (lin_eq[1][2].fi == 0)
                continue;
            return true;
        }
        if (lin_eq[2][2].fi == lin_eq[2][2].se) {
            // 1 0 0 *
            // 0 1 0 *
            // 0 0 1 *
            Frac det = lin_eq[0][3] - lin_eq[1][3] * lin_eq[2][3];
            if (det.fi != 0) {
                return true;
            }
            continue;
        }
        if (lin_eq[0][0].fi == 0) {
            // 0 1 0 *     0 1 * 0     0 1 * *     0 0 1 0     0 0 1 0    0 0 0 1    0 0 0 0
            // 0 0 1 *     0 0 0 1     0 0 0 0     0 0 0 1     0 0 0 0    0 0 0 0    0 0 0 0
            if (lin_eq[1][2].fi == 0 && lin_eq[1][3].fi == lin_eq[1][3].se) {
                if (lin_eq[0][1].fi == lin_eq[0][1].se && lin_eq[0][2].fi != 0) {
                    return true;
                }
                continue;
            }
            return true;
        }
        // 1 0 * *     1 * 0 *     1 * * 0   1 * * *
        // 0 1 * *     0 0 1 *     0 0 0 1   0 0 0 0
        if (lin_eq[1][1].fi == lin_eq[1][1].se) {
            // 1 0 * *
            // 0 1 * * 
            if (lin_eq[0][3].fi == 0 && lin_eq[1][2].fi == 0 && lin_eq[0][2] == lin_eq[1][3]) {
                continue;
            }
            return true;
        }
        if (lin_eq[1][2].fi == lin_eq[1][2].se) {
            // 1 * 0 *
            // 0 0 1 *
            if (lin_eq[0][1].fi == 0 && (lin_eq[0][3] == lin_eq[1][3])) {
                continue;
            }
            return true;
        }
        return true;

    }
    return false;
}


int main() {
    
    vec<Matrix> Akbari = {
        {{{1, 1}, {0, 1}},
         {{0, 1}, {1, 1}}},
        
        {{{0, 1}, {1, 1}},
         {{1, 1}, {0, 1}}},
        
        {{{0, 1}, {-1, 1}},
         {{-1, 1}, {0, 1}}},

        {{{0, 1}, {1, 1}},
         {{-1, 1}, {-2, 1}}},
        
        {{{0, 1}, {-1, 1}},
         {{1, 1}, {-2, 1}}}
    };

    vec<Matrix> new_sample = {
        {{{2, 1}, {2, 1}},
         {{2, 1}, {1, 1}}},
        
        {{{2, 1}, {2, 1}},
         {{2, 1}, {3, 1}}},
        
        {{{-2, 1}, {2, 1}},
         {{-2, 1}, {1, 1}}},

        {{{-2, 1}, {-2, 1}},
         {{2, 1}, {-1, 1}}},
        
        {{{-2, 1}, {-2, 1}},
         {{-4, 1}, {-1, 1}}}
    };

    if (isomorphic_subgraphs_2_by_2(Akbari, new_sample)) 
        cout << "isomorphic subgraphs" << endl;
    else
        cout << "not isomorphic subgraphs" << endl;
    
}
