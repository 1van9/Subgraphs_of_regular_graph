#pragma once
// for working with polynomials

#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>
#define vec std::vector

template<typename T>
struct Poly {
    void update() {
        while (!coefficients.empty() && !coefficients.back()) 
            coefficients.pop_back();
    }
    size_t degree() const {
        return coefficients.size();
    }
    Poly() {}
    Poly(T x) {
        if (x)
            coefficients.push_back(x);
    } 
    Poly(long long int x) {
        if (x)
            coefficients.emplace_back(x);
    }
    Poly(const vec<T> & c) : coefficients(c) {
        update();
    }
    Poly(const vec<std::pair<size_t, T>> & c) {
        size_t max_degree = std::max_element(c.begin(), c.end()) -> first;
        coefficients.resize(max_degree + 1);
        for (auto [deg, coeff] : c) 
            coefficients[deg] = coeff;
        update();
    }
    Poly(const Poly & p) : coefficients(p.coefficients) {}
    Poly operator = (const Poly & p) {
        coefficients = p.coefficients;
        return p;
    }
    Poly operator + (const Poly & p) const {
        size_t max_degree = std::max(degree(), p.degree());
        Poly c;
        c.coefficients.resize(max_degree);
        for (size_t i = 0; i < max_degree; i++) {
            if (i < degree()) 
                c.coefficients[i] += coefficients[i];
            if (i < p.degree())
                c.coefficients[i] += p.coefficients[i];
        }
        c.update();
        return c;
    }
    Poly & operator += (const Poly & p) {
        for (size_t i = 0; i < degree() && i < p.degree(); i++)
            coefficients[i] += p.coefficients[i];
        while (degree() < p.degree())
            coefficients.push_back(p.coefficients[degree()]);
        update();
        return *this;
    }
    Poly operator - (const Poly & p) const {
        size_t max_degree = max(degree(), p.degree());
        Poly c;
        c.coefficients.resize(max_degree);
        for (size_t i = 0; i < max_degree; i++) {
            if (i < degree()) 
                c.coefficients[i] += coefficients[i];
            if (i < p.degree())
                c.coefficients[i] -= p.coefficients[i];
        }
        c.update();
        return c;
    }
    Poly operator - () const {
        Poly c;
        c.coefficients.resize(degree());
        for (size_t i = 0; i < degree(); i++) {
            c.coefficients[i] = -coefficients[i];
        }
        c.update();
        return c;
    }
    Poly & operator -= (const Poly & p) {
        for (int i = 0; i < degree() && i < p.degree(); i++)
            coefficients[i] -= p.coefficients[i];
        while (degree() < p.degree())
            coefficients.push_back(-p.coefficients[degree()]);
        update();
        return *this;
    }
    Poly operator * (const Poly & p) const {
        Poly c;
        for (size_t i = 0; i < degree(); i++) {
            if (coefficients[i]) {
                Poly add;
                add.coefficients.resize(p.degree() + i);
                for (size_t j = 0; j < p.degree(); j++)
                    add.coefficients[j + i] = p.coefficients[j] * coefficients[i];
                c += add;
            }
        }
        return c;
    }
    Poly & operator *= (const Poly & p) {
        *this = (*this) * p;
        return *this;
    }
    std::pair<Poly, Poly> long_division(const Poly & p) const {
        Poly remainder = *this;
        Poly quotient;
        if (p.degree() == 0)
            assert(0);
        quotient.coefficients.resize(degree());
        for (size_t i = degree(); i >= p.degree(); i--) {
            if (remainder.degree() < i) 
                continue;
            T coeff = remainder.coefficients[i - 1] / p.coefficients[p.degree() - 1];
            quotient.coefficients[i - p.degree()] = coeff;
            remainder -= (p * Poly({make_pair(i - p.degree(), coeff)}));
        }
        quotient.update();
        remainder.update();
        return {quotient, remainder};
    }
    Poly operator / (const Poly & p) const {
        return long_division(p).first;
    }
    Poly & operator /= (const Poly & p) {
        *this = long_division(p).first;
        return *this;
    }
    Poly operator % (const Poly & p) const {
        return long_division(p).second;
    }
    Poly & operator %= (const Poly & p) {
        (*this) = long_division(p).second;
        return *this;
    }
    T operator () (T x) const {
        T result = T(0);
        T pw = T(1);
        for (auto coeff : coefficients) {
            result += coeff * pw;
            pw *= x;
        }
        return result;
    }
    T operator [] (size_t i) const {
        if (i >= degree() || !coefficients[i])
            return T(0);
        return coefficients[i];
    }
    T & operator [] (size_t i) {
        return coefficients[i];
    }
    bool operator == (const Poly & p) const {
        return coefficients == p.coefficients;
    }
    bool operator != (const Poly & p) const {
        return coefficients != p.coefficients;
    }
    bool operator == (const T & x) const {
        if (x)
            return degree() == 1 && coefficients[0] == x;
        return !degree();
    }
    bool operator != (const T& x) const {
        if (x)
            return degree() != 1 || coefficients[0] != x;
        return degree();
    }
    bool operator > (const T & x) const {
        return (degree() ? coefficients.back() > x : x < T(0));
    }
    bool operator >= (const T & x) const {
        return (degree() ? coefficients.back() >= x : x <= T(0));
    }
    bool operator < (const T & x) const {
        return (degree() ? coefficients.back() < x : x > T(0));
    }
    bool operator <= (const T & x) const {
        return (degree() ? coefficients.back() <= x : x >= T(0));
    }
    bool operator > (const Poly & p) const {
        T x = (degree() ? coefficients.back() : T(0));
        return p < x;
    }
    bool operator >= (const Poly & p) const {
        T x = (degree() ? coefficients.back() : T(0));
        return p <= x;
    }
    bool operator < (const Poly & p) const {
        T x = (degree() ? coefficients.back() : T(0));
        return p > x;
    }
    bool operator <= (const Poly & p) const {
        T x = (degree() ? coefficients.back() : T(0));
        return p >= x;
    }
    operator bool() const {
        return degree();
    }

    vec<T> coefficients;
};

template<typename T>
char PSign(T a) {
    if (a >= T(0))
        return '+';
    return '-';
}
template<typename T>
T abs(T a) {
    if (a >= T(0))
        return a;
    return -a;
}

template<typename T>
std::ostream & operator << (std::ostream & out, const Poly<T> & p) {
    if (p.degree() == 0) {
        out << "0";
        return out;
    }
    bool first = true;
    if (p.coefficients[0])
        first = false, out << p.coefficients[0];
    for (size_t i = 1; i < p.degree(); i++) {
        if (!p.coefficients[i]) 
            continue;
        if (first) {
            if (p.coefficients[i] == -1ll)
                out << '-';                
            if (abs(p.coefficients[i]) != 1ll)
                out << p.coefficients[i];
            first = false;
        } else {
            out << " " << PSign(p.coefficients[i]) << " ";
            if (abs(p.coefficients[i]) != 1ll)
                out << abs(p.coefficients[i]);
        }
        out << 'x';
        if (i > 1) 
            out << '^' << i; 
    }
    return out;
}

template<typename T>
void swap(Poly<T> & p, Poly<T> & q) {
    swap(q.coefficients, p.coefficients);
}

template<typename T>
Poly<T> sqrt(const Poly<T> & a) {
    if (a.degree() % 2 == 0)
        return a;
    if (!a.degree())
        return a;
    if (a.degree() == 1)
        return Poly<T>(sqrt(a[0]));
    size_t len = a.degree() / 2;
    vec<T> b(len + 1);
    T last = sqrt(a.coefficients.back());
    b[len] = last;
    Poly<T> res(b);
    Poly<T> curr = a - res * res;
    for (size_t i = len; i > 0; i--) {
        res[i - 1] = curr[i + len - 1] / (last + last);
        curr = a - res * res;
    }
    return res;
}
