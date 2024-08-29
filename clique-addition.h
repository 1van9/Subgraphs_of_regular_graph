#pragma once
#include <iostream>
#include <vector>
#include <cassert>
#include <set>
#include <cmath>
#include "matrix.h"


template<typename T>
vec<T> solution(T A, T B, T C, char var) {
    if (!A) {
        if (!B) {
            return {};
        }
        return {-C / B};
    }
    T D = B * B - (A * C * T(4));
    // var = -B +- sqrt(D) / 2A
    if (!D) {
        return {-B / (A + A)};
    }
    if (D < T(0)) {
        return {};
    }
    T D_ = sqrt(D);
    if (D_ * D_ == D) {
        T x1 = (-B + D_) / (A + A);
        T x2 = (-B - D_) / (A + A);
        return {x1, x2};
    } else {
        //printf("Not in Q, but maybe in R\nequasion: %d/%d %c^2 + %d/%d %c + %d/%d = 0\n", A.fi, A.se, var, B.fi, B.se, var, C.fi, C.se);
        return {};
    }
}

template<typename T>
inline vec<T> through_x(T y0, T y1, T z0, T z1, T t0, T t1, const Matrix<T> & a) {
    // x (t0 + t1x) - (y0 + y1x)(z0 + z1x) + a[1][1]x - a[1][0](y0 - y1x) - a[0][1](z0 + z1x) + a[0][0](t0 + t1x) + det(A) = 0
    T A = t1 - y1 * z1;
    T B = t0 - y1 * z0 - z1 * y0 + a[1][1] - a[1][0] * y1 - a[0][1] * z1 + a[0][0] * t1;
    T C = - y0 * z0 - a[1][0] * y0 - a[0][1] * z0 + a[0][0] * t0 + a.det();
    // Ax^2 + Bx + C
    return solution(A, B, C, 'x');
}

template<typename T>
inline vec<T> through_y(T x0, T x1, T z0, T z1, T t0, T t1, const Matrix<T> & a) {
    // (x0 + x1y) (t0 + t1y) - y(z0 + z1y) + a[1][1](x0 + x1y) - a[1][0]y - a[0][1](z0 + z1y) + a[0][0](t0 + t1y) + det(A) = 0
    T A = x1 * t1 - z1;
    T B = t1 * x0 + x1 * t0 - z0 + a[1][1] * x1 - a[1][0] - a[0][1] * z1 + a[0][0] * t1;
    T C = x0 * t0 + a[1][1] * x0 - a[0][1] * z0 + a[0][0] * t0 + a.det();
    // Ay^2 + By + C
    return solution(A, B, C, 'y');
}

template<typename T>
inline vec<T> through_z(T x0, T x1, T y0, T y1, T t0, T t1, const Matrix<T> & a) {
    // (x0 + x1z)(t0 + t1z) - (y0 + y1z)z + a[1][1](x0 + x1z) - a[1][0](y0 - y1z) - a[0][1]z + a[0][0](t0 + t1z) + det(A) = 0
    T A = -y1 + x1 * t1;
    T B = -y0 + x1 * t0 + t1 * x0 + a[1][1] * x1 - a[1][0] * y1 - a[0][1] + a[0][0] * t1;
    T C = x0 * t0 + a[1][1] * x0 - a[1][0] * y0 + a[0][0] * t0 + a.det();
    // Az^2 + Bz + C
    return solution(A, B, C, 'z');
}

template<typename T>
inline vec<T> through_t(T x0, T x1, T y0, T y1, T z0, T z1, const Matrix<T> & a) {
    // (x0 + x1t)t  - (y0 + y1t)(z0 + z1t) + a[1][1](x0 + x1t) - a[1][0](y0 + y1t) - a[0][1](z0 + z1t) + a[0][0]t + det(A) = 0
    T A = -y1 * z1 + x1;
    T B = -z1 * y0 - y1 * z0 + x0 + a[1][1] * x1 - a[1][0] * y1 - a[0][1] * z1 + a[0][0];
    T C = -y0 * z0 + a[1][1] * x0 - a[1][0] * y0 - a[0][1] * z0  + a.det();
    // At^2 + Bt + C = 0
    return solution(A, B, C, 't');
}

template<typename T>
inline vec<Matrix<T>> check(vec<T> & xs, vec<T> & ys, vec<T> & zs, vec<T> & ts, const vec<Matrix<T>> & clique) {
    std::set<Matrix<T>> answer;
    for (auto x : xs)
        for (auto y : ys)
            for (auto z : zs)
                for (auto t : ts) {
                    Matrix<T> A({{x, y}, {z, t}});
                    if (!A.det())
                        continue;
                    bool ans = true;
                    for (auto X : clique)
                        if (!A.Connect(X)) {
                            ans = false;
                            break; 
                        }
                    if (ans)
                        answer.insert(A);
                }
    vec<Matrix<T>> ans(answer.begin(), answer.end());
    return ans;
}

// addition K4 to K(1, 1, 1, 1, 2)|K(1, 1, 1, 1, 1)|K(1, 1, 1, 1, 0)
template<typename T>
vec<Matrix<T>> is_additionable(const vec<Matrix<T>> & clique) {
    assert(clique.size() == 4);
    // a[0][0] a[0][1]
    // a[1][0] a[1][1]
    // equations:
    // xz - yt + a[1][1]x - a[1][0]y - a[0][1]z + a[0][0]t + det(A) = 0
    //              ___ x -    ___ y -    ___ z +    ___ t + ___    = 0
    //              ___ x -    ___ y -    ___ z +    ___ t + ___    = 0
    //              ___ x -    ___ y -    ___ z +    ___ t + ___    = 0
    Matrix<T> linear_equations(3, 5);
    T d = clique[0].det();
    for (int i = 1; i < 4; i++) {
        linear_equations[i - 1][0] = clique[i][1][1] - clique[0][1][1];
        linear_equations[i - 1][1] = clique[0][1][0] - clique[i][1][0];
        linear_equations[i - 1][2] = clique[0][0][1] - clique[i][0][1];
        linear_equations[i - 1][3] = clique[i][0][0] - clique[0][0][0];
        linear_equations[i - 1][4] = clique[i].det() - d;
    }
    // in clique this linear equasions always independent
    linear_equations.Gause();
    if (!linear_equations[2][0] && !linear_equations[2][1] &&
        !linear_equations[2][2] && !linear_equations[2][3]) {    
        // * * * * *     
        // * * * * *
        // 0 0 0 0 *
        return {};
    }
    
    int val = 3;
    for (int i = 0; i < 3; i++) {
        if (!linear_equations[i][i]) {
            val = i;
            break;
        }
    }
    vec<T> xs, ys, zs, ts;
    if (val == 0) {
        // 0 1 0 0 *
        // 0 0 1 0 *     val = 0
        // 0 0 0 1 *
        xs = through_x(-linear_equations[0][4], T(0), 
                       -linear_equations[1][4], T(0),
                       -linear_equations[2][4], T(0),
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
                       -linear_equations[1][4], T(0),
                       -linear_equations[2][4], T(0),
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
                       -linear_equations[2][4], T(0),
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
    vec<Matrix<T>> M = check(xs, ys, zs, ts, clique);
    return M;
}

// how many ways to expand K5 to K(1, 1, 1, 1, 2)
template<typename T>
int addition(const vec<Matrix<T>> & clique) {
    int cnt = 0;
    for (int i = 0; i < 5; i++) {
        vec<Matrix<T>> a = clique;
        a.erase(a.begin() + i);
        vec<Matrix<T>> add = is_additionable(a);
        if (add.size() > 1)
            cnt++;
    }
    return cnt;
}
