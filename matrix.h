#pragma once
// header for matrices over any feild

#include <iostream>
#include <vector>
#include <cassert>
#define vec std::vector

template<typename T>
struct Matrix {
    Matrix() {}
    Matrix(size_t nn) : n(nn), m(nn) {
        a.resize(n, vec<T>(m));
    }
    Matrix(size_t nn, size_t mm) : n(nn), m(mm) {
        a.resize(n, vec<T>(m));
    }
    Matrix(const vec<vec<T>> & ar) : a(ar) {
        n = a.size();
        m = (n) ? a[0].size() : 0;
    }
    Matrix(const Matrix & b) : n(b.n), m(b.m), a(b.a) {}
    Matrix operator = (const Matrix & b) {
        n = b.n, m = b.m;
        a = b.a;
        return b;
    }
    Matrix operator + (const Matrix & b) const {
        Matrix c(n, m);
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < m; j++)
                c.a[i][j] = a[i][j] + b.a[i][j];
        return c;
    }
    Matrix operator - (const Matrix & b) const {
        Matrix c(n, m);
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < m; j++)
                c.a[i][j] = a[i][j] - b.a[i][j];
        return c;
    }
    Matrix operator - () const {
        Matrix c(n, m);
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < m; j++)
                c.a[i][j] -= a[i][j];
        return c;
    }
    Matrix operator += (const Matrix & b) {
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < m; j++)
                a[i][j] += b.a[i][j];
        return a;
    }
    Matrix operator -= (const Matrix & b) {
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < m; j++)
                a[i][j] -= b.a[i][j];
        return a;
    }
    Matrix operator * (const Matrix & b) const {
        Matrix c(n, b.m);
        for (size_t i = 0; i < n; i++) 
            for (size_t j = 0; j < b.m; j++) 
                for (size_t k = 0; k < m; k++)
                    c.a[i][j] += a[i][k] * b.a[k][j];
        return c;                
    }
    Matrix operator * (const T & x) const {
        Matrix c(n, m);
        for (size_t i = 0; i < n; i++) 
            for (size_t j = 0; j < m; j++) 
                c.a[i][j] = a[i][j] * x;
        return c;
    }
    Matrix operator *= (const Matrix & b) {
        Matrix c(n, b.m);
        for (size_t i = 0; i < n; i++) 
            for (size_t j = 0; j < b.m; j++) 
                for (size_t k = 0; k < m; k++)
                    c.a[i][j] += a[i][k] * b.a[k][j];
        return a = c;                
    }
    Matrix operator *= (const T & x) {
        for (size_t i = 0; i < n; i++) 
            for (size_t j = 0; j < m; j++) 
                a[i][j] *= x;
        return a;
    }
    bool operator < (const Matrix & b) const {
        return a < b.a;
    }
    bool operator <= (const Matrix & b) const {
        return a <= b.a;
    }
    bool operator > (const Matrix & b) const {
        return a > b.a;
    }
    bool operator >= (const Matrix & b) const {
        return a >= b.a;
    }
    bool operator == (const Matrix & b) const {
        return a == b.a;
    }
    bool operator != (const Matrix & b) const {
        return a != b.a;
    }
    vec<T> & operator [] (size_t i) {
        return a[i];
    }
    const vec<T> & operator [] (size_t i) const {
        return a[i];
    }
    size_t size() {
        return n;
    }
    Matrix t() {
        Matrix b(m, n);
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                b.a[i][j] = a[j][i];
        return b;
    }
    T Gause() {
        size_t curr_string = 0;
        T det;
        det = 1;
        for (size_t i = 0; i < m && curr_string < n; i++) {
            size_t curr = curr_string;
            while (curr < n && !a[curr][i])
                curr++;
            if (curr == n) 
                continue;
            if (curr != curr_string)
                for (size_t j = 0; j < m; j++)
                    swap(a[curr][j], a[curr_string][j]);
            T div = a[curr_string][i];
            det *= div;
            for (size_t j = 0; j < m; j++)
                a[curr_string][j] /= div;
            for (size_t j = 0; j < n; j++) {
                if (j == curr_string)
                    continue;
                T mul = a[j][i];
                for (size_t k = 0; k < m; k++)
                    a[j][k] -= a[curr_string][k] * mul;
            }
            curr_string++;
        }
        if (curr_string != n) 
            det = 0;
        return det;
    }
    T det() const {
        if (n == 2 && m == 2)
            return a[0][0] * a[1][1] - a[0][1] * a[1][0];
        Matrix copy_a = a;
        return copy_a.Gause();
    }
    Matrix rev() const {
        Matrix b(n, 2 * n);
        for (size_t i = 0; i < n; i++) {
            b.a[i][i + n] = 1;
            for (size_t j = 0; j < n; j++)
                b.a[i][j] = a[i][j];
        }
        T det = b.Gause();
        assert(det);
        Matrix res(n, n);
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                res.a[i][j] = b.a[i][j + n];
        return res;
    }
    bool GL() const {
        return det();
    }
    bool Connect(const Matrix & b) const {
        if (n == 2)
            return !((a[0][0] + b.a[0][0]) * (a[1][1] + b.a[1][1]) - (a[0][1] + b.a[0][1]) * (a[1][0] + b.a[1][0]));
        return !(*this + b).det();
    }

    vec<vec<T>> a;
    size_t n, m;
};

template<typename T>
Matrix<T> operator * (const T & x, const Matrix<T> & a) {
    Matrix<T> c(a.n, a.m);
    for (size_t i = 0; i < a.n; i++) 
        for (size_t j = 0; j < a.m; j++) 
            c.a[i][j] = a[i][j] * x;
    return c;
}

template <typename T>
std::ostream & operator << (std::ostream & out, const Matrix<T> & a) {
    for (size_t i = 0; i < a.n; i++) {
        for (auto el : a.a[i])
            out << el << " ";
        out << std::endl;
    }
    return out; 
}
