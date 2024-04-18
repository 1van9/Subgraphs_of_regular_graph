#pragma once
// header for matrices over field Q

#include <iostream>
#include <vector>
#include <cassert>

#define vec std::vector
#define fi first
#define se second
#define make_frac(x, y) upd(std::make_pair(x, y))

using Frac = std::pair<long long, long long>;
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
    return make_frac(a.fi * b.se + b.fi * a.se, a.se * b.se);
}
Frac operator+= (Frac & a, const Frac & b) {
    return a = make_frac(a.fi * b.se + b.fi * a.se, a.se * b.se);
}
Frac operator- (const Frac & a, const Frac & b) {
    return make_frac(a.fi * b.se - b.fi * a.se, a.se * b.se);
}
Frac operator-= (Frac & a, const Frac & b) {
    return a = make_frac(a.fi * b.se - b.fi * a.se, a.se * b.se);
}
Frac operator* (const Frac & a, const Frac & b) {
    return make_frac(a.fi * b.fi, a.se * b.se);
}
Frac operator*= (Frac & a, const Frac & b) {
    return a = make_frac(a.fi * b.fi, a.se * b.se);
}
Frac operator/ (const Frac & a, const Frac & b) {
    return make_frac(a.fi * b.se, a.se * b.fi);
}
Frac operator/= (Frac & a, const Frac & b) {
    return a = make_frac(a.fi * b.se, a.se * b.fi);
}
bool operator!=(const Frac & a, const Frac & b) {
    return a.fi * b.se != a.se * b.fi;
}
bool operator==(const Frac & a, const Frac & b) {
    return a.fi * b.se == a.se * b.fi;
}
std::ostream& operator<< (std::ostream& out, const Frac & a) {
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

Matrix operator+ (const Matrix & a, const Matrix & b) {
    int n = a.size(), m = a[0].size();
    Matrix c(n, vec<Frac>(m));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            c[i][j] = a[i][j] + b[i][j];
    return c;
}
Matrix operator- (const Matrix & a, const Matrix & b) {
    int n = a.size(), m = a[0].size();
    Matrix c(n, vec<Frac>(m));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            c[i][j] = a[i][j] - b[i][j];
    return c;
}
Matrix operator+= (Matrix & a, const Matrix & b) {
    int n = a.size(), m = a[0].size();
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            a[i][j] += b[i][j];
    return a;
}
Matrix operator-= (Matrix & a, const Matrix & b) {
    int n = a.size(), m = a[0].size();
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            a[i][j] -= b[i][j];
    return a;
}
Matrix operator* (const Matrix & a, const Matrix & b) {
    int n = a.size(), m = a[0].size(), k = b[0].size();
    Matrix c(n, vec<Frac>(k, {0, 1}));
    for (int i = 0; i < n; i++) 
        for (int j = 0; j < k; j++) 
            for (int l = 0; l < m; l++)
                c[i][j] += a[i][l] * b[l][j];
    return c;                
}
Matrix operator*= (Matrix & a, const Matrix & b) {
    int n = a.size(), m = a[0].size(), k = b[0].size();
    Matrix c(n, vec<Frac>(k, {0, 1}));
    for (int i = 0; i < n; i++) 
        for (int j = 0; j < k; j++) 
            for (int l = 0; l < m; l++)
                c[i][j] += a[i][l] * b[l][j];
    return a = c;                
}
std::ostream& operator<< (std::ostream& out, const Matrix & a) {
    for (int i = 0; i < a.size(); i++) {
        for (auto el : a[i]) {
            out << el << " ";
        }
        out << std::endl;
    }
    return out; 
}

Matrix T(const Matrix & a) {
    int n = a.size(), m = a[0].size();
    Matrix b(m, vec<Frac>(n));
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            b[i][j] = a[j][i];
    return b;
}

Frac Gause(Matrix& a) {
    int n = a.size(), m = a[0].size();

    int curr_string = 0;
    Frac det = {1, 1};
    for (int i = 0; i < m && curr_string < n; i++) {
        int curr = curr_string;
        while (curr < n && a[curr][i].fi == 0)
            curr++;
        if (curr == n) continue;
        if (curr != curr_string) {
            for (int j = 0; j < m; j++) {
                swap(a[curr][j], a[curr_string][j]);
            }
        }
        Frac div = a[curr_string][i];
        det /= div;
        for (int j = 0; j < m; j++) {
            a[curr_string][j] /= div;
        }
        for (int j = 0; j < n; j++) {
            if (j == curr_string)
                continue;
            Frac mul = a[j][i];
            for (int k = 0; k < m; k++) {
                a[j][k] -= a[curr_string][k] * mul;
            }
        }
        curr_string++;
    }
    if (curr_string != n || n != m) 
        det = {0, 1};
    return det;
}

Matrix rev(const Matrix& a) {
    int n = a.size();
    Matrix b(n, vec<Frac> (2 * n, {0, 1}));
    for (int i = 0; i < n; i++) {
        b[i][i + n] = {1, 1};
        for (int j = 0; j < n; j++) {
            b[i][j] = a[i][j];
        }
    }
    Gause(b);
    bool check = true;
    for (int i = 0; i < n; i++) {
        if (b[i][i].fi != b[i][i].se) {
            check = false;
        }
    }
    if (!check) {
        assert(0);
    }
    Matrix res(n, vec<Frac>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            res[i][j] = b[i][j + n];
        }
    }
    return res;
}

Frac det(const Matrix & a) {
    if (a.size() == 2 && a[0].size() == 2)
        return a[0][0] * a[1][1] - a[0][1] * a[1][0];
    Matrix copy_a = a;
    return Gause(copy_a);
}

bool Connect(const Matrix & a, const Matrix & b) {
    return det(a + b).fi == 0;
}

bool GL(const Matrix & a) {
    return det(a).fi != 0;
}