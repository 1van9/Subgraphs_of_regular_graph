#pragma once

#include "matrices.h"
#include <set>
#include <algorithm>

bool isomorphic_subgraphs_2_by_2(const vec<Matrix>& subg1, const vec<Matrix>& subg2) {
    if (subg1.size() != subg2.size())
        return false;
    int n = subg1.size();
    vec<int> p(n);
    for (int i = 0; i < n; i++)
        p[i] = i;
    bool ok = false;
    // A_i = U B_{p_i} V
    // U = B_{p_i} V (A_i)^-1   (1)

    // V = (x y
    //      z t)
    // Have n wasy, how to express U, in term of x, y, z, t
    std::set<Matrix> ways; // possible linear equasions on coefficients x, y, z, t
    do {
        vec<Matrix> M; // express U coefficients in term of x, y, z, t 
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

    // the same for tansposition
    do {
        vec<Matrix> M;
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
