#include <iostream>
#include <random>
#include <vector>
#include <cassert>
#include <algorithm>
#include <set>
#include "matrices.h"

using namespace std;

bool only_R = false;

vector<Frac> solution(Frac A, Frac B, Frac C, char var) {
    if (A.fi == 0) {
        if (B.fi == 0) {
            return {};
        }
        return {-C / B};
    }
    Frac D = B * B - (make_frac(4, 1) * A * C);
    // var = -B +- sqrt(D) / 2A
    if (D.fi == 0) {
        return {-B / (A + A)};
    }
    if (D.fi < 0) {
        return {};
    }
    long long x = sqrt(D.fi), y = sqrt(D.se);
    if (x * x == D.fi && y * y == D.se) {
        Frac D_ = make_frac(x, y);
        Frac x1 = (-B + D_) / (A + A);
        Frac x2 = (-B - D_) / (A + A);
        return {x1, x2};
    } else {
        cout << "Not in Q, but maybe in R" << endl;
        printf("equasion: %d/%d %c^2 + %d/%d %c + %d/%d = 0\n", A.fi, A.se, var, B.fi, B.se, var, C.fi, C.se);
        only_R = true;
        return {};
    }
}

inline vector<Frac> through_x(Frac y0, Frac y1, Frac z0, Frac z1, Frac t0, Frac t1, const Matrix & a) {
    // x (t0 + t1x) - (y0 + y1x)(z0 + z1x) + a[1][1]x - a[1][0](y0 - y1x) - a[0][1](z0 + z1x) + a[0][0](t0 + t1x) + det(A) = 0
    Frac A = t1 - y1 * z1;
    Frac B = t0 - y1 * z0 - z1 * y0 + a[1][1] - a[1][0] * y1 - a[0][1] * z1 + a[0][0] * t1;
    Frac C = - y0 * z0 - a[1][0] * y0 - a[0][1] * z0 + a[0][0] * t0 + det(a);
    // Ax^2 + Bx + C
    return solution(A, B, C, 'x');
}

inline vector<Frac> through_y(Frac x0, Frac x1, Frac z0, Frac z1, Frac t0, Frac t1, const Matrix & a) {
    // (x0 + x1y) (t0 + t1y) - y(z0 + z1y) + a[1][1](x0 + x1y) - a[1][0]y - a[0][1](z0 + z1y) + a[0][0](t0 + t1y) + det(A) = 0
    Frac A = x1 * t1 - z1;
    Frac B = t1 * x0 + x1 * t0 - z0 + a[1][1] * x1 - a[1][0] - a[0][1] * z1 + a[0][0] * t1;
    Frac C = x0 * t0 + a[1][1] * x0 - a[0][1] * z0 + a[0][0] * t0 + det(a);
    // Ay^2 + By + C
    return solution(A, B, C, 'y');
}

inline vector<Frac> through_z(Frac x0, Frac x1, Frac y0, Frac y1, Frac t0, Frac t1, const Matrix & a) {
    // (x0 + x1z)(t0 + t1z) - (y0 + y1z)z + a[1][1](x0 + x1z) - a[1][0](y0 - y1z) - a[0][1]z + a[0][0](t0 + t1z) + det(A) = 0
    Frac A = -y1 + x1 * t1;
    Frac B = -y0 + x1 * t0 + t1 * x0 + a[1][1] * x1 - a[1][0] * y1 - a[0][1] + a[0][0] * t1;
    Frac C = x0 * t0 + a[1][1] * x0 - a[1][0] * y0 + a[0][0] * t0 + det(a);
    // Az^2 + Bz + C
    return solution(A, B, C, 'z');
}

inline vector<Frac> through_t(Frac x0, Frac x1, Frac y0, Frac y1, Frac z0, Frac z1, const Matrix & a) {
    // (x0 + x1t)t  - (y0 + y1t)(z0 + z1t) + a[1][1](x0 + x1t) - a[1][0](y0 + y1t) - a[0][1](z0 + z1t) + a[0][0]t + det(A) = 0
    Frac A = -y1 * z1 + x1;
    Frac B = -z1 * y0 - y1 * z0 + x0 + a[1][1] * x1 - a[1][0] * y1 - a[0][1] * z1 + a[0][0];
    Frac C = -y0 * z0 + a[1][1] * x0 - a[1][0] * y0 - a[0][1] * z0  + det(a);
    // At^2 + Bt + C = 0
    return solution(A, B, C, 't');
}

inline vector<Matrix> check(vector<Frac> & xs, vector<Frac> & ys, vector<Frac> & zs, vector<Frac> & ts, vector<Matrix> & clique) {
    set<Matrix> answer;
    for (auto x : xs)
        for (auto y : ys)
            for (auto z : zs)
                for (auto t : ts) {
                    Matrix A = {{x, y}, {z, t}};
                    if (det(A).fi == 0)
                        continue;
                    bool ans = true;
                    for (auto X : clique)
                        if (det(A + X).fi != 0) {
                            ans = false;
                            break; 
                        }
                    if (ans)
                        answer.insert(A);
                }
    vector<Matrix> ans(answer.begin(), answer.end());
    return ans;
}

// addition K(1, 1, 1, 1) to K(1, 1, 1, 1, 2)|K(1, 1, 1, 1, 1)|K(1, 1, 1, 1, 0)
vector<Matrix> is_additionable(vector<Matrix> & clique) {
    assert(clique.size() == 4);
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
    // in clique this linear equasions always independent
    Gause(linear_equations);
    if (linear_equations[2][0].fi == 0 && linear_equations[2][1].fi == 0 &&
        linear_equations[2][2].fi == 0 && linear_equations[2][3].fi == 0) {    
        // * * * * *     
        // * * * * *
        // 0 0 0 0 *
        return {};
    }
    
    int val = 3;
    for (int i = 0; i < 3; i++) {
        if (linear_equations[i][i].fi == 0) {
            val = i;
            break;
        }
    }
    vector<Frac> xs, ys, zs, ts;
    if (val == 0) {
        // 0 1 0 0 *
        // 0 0 1 0 *     val = 0
        // 0 0 0 1 *
        xs = through_x(-linear_equations[0][4], make_frac(0, 1), 
                       -linear_equations[1][4], make_frac(0, 1),
                       -linear_equations[2][4], make_frac(0, 1),
                       clique[0]);
        ys = {-linear_equations[0][4]};
        zs = {-linear_equations[1][4]};
        ts = {-linear_equations[2][4]};
    }
    if (val == 1) {     
        // 1 * 0 0 *
        // 0 0 1 0 *     val = 1
        // 0 0 0 1 *
        ys = through_y(-linear_equations[0][4], -linear_equations[0][1], 
                       -linear_equations[1][4], make_frac(0, 1),
                       -linear_equations[2][4], make_frac(0, 1),
                       clique[0]);
        for (auto el : ys) {
            xs.push_back(-linear_equations[0][4] - linear_equations[0][1] * el);
        }
        zs = {-linear_equations[1][4]};
        ts = {-linear_equations[2][4]};
    }
    if (val == 2) {
        // 1 0 -a 0 -c 
        // 0 1 -b 0 -d   t = e; z; y = d + bz; x = c + az;
        // 0 0  0 1 -e
        zs = through_z(-linear_equations[0][4], -linear_equations[0][2],
                       -linear_equations[1][4], -linear_equations[1][2],
                       -linear_equations[2][4], make_frac(0, 1),
                       clique[0]);
        for (auto z : zs) { 
            xs.push_back(-linear_equations[0][4] - linear_equations[0][2] * z);
            ys.push_back(-linear_equations[1][4] - linear_equations[1][2] * z);
        }
        ts.push_back(-linear_equations[2][4]);
    }
    if (val == 3) {
        // 1 0 0 -a -d
        // 0 1 0 -b -e   t; z = f + ct; y = e + bt; x = d + at;  
        // 0 0 1 -c -f
        ts = through_t(-linear_equations[0][4], -linear_equations[0][3],
                       -linear_equations[1][4], -linear_equations[1][3],
                       -linear_equations[2][4], -linear_equations[2][3],
                       clique[0]);
        for (auto t : ts) { 
            xs.push_back(-linear_equations[0][4] - linear_equations[0][3] * t);
            ys.push_back(-linear_equations[1][4] - linear_equations[1][3] * t);
            zs.push_back(-linear_equations[2][4] - linear_equations[2][3] * t);
        }
    }
    vector<Matrix> M = check(xs, ys, zs, ts, clique);
    return M;
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

    for (int i = 0; i < 5; i++) {
        vector<Matrix> a = Akbari;

        a.erase(a.begin() + i);

        vector<Matrix> add = is_additionable(a);
        cout << add.size() << endl;
    }
}
