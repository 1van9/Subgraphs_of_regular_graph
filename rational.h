#pragma ones
#include <iostream>
#include <cassert>
#include <cmath>

template<typename T>
T gcd(T a, T b) {
    if (!b) return a;
    return gcd(b, a % b);
}


template<typename T>
struct Frac {

    void upd() {
        if (!fi && !se)
            return;
        T d = gcd(abs(fi), se);
        fi /= d, se /= d;
        if (se < T(0)) 
            fi = -fi, se = -se;
    }
    Frac() : fi(0), se(1) {}
    Frac(const T& x, const T& y) : fi(x), se(y) {
        upd();
    }
    Frac(const T& x) : fi(x), se(1) {
        upd();
    }
    Frac(const Frac & x) : fi(x.fi), se(x.se) {}
    Frac operator=(const Frac & x) {
        fi = x.fi, se = x.se;
        return x;
    }
    Frac operator= (const T & x) {
        fi = x;
        se = T(1);
        return x;
    }

    Frac operator+ (const Frac & b) const {
        return Frac(fi * b.se + b.fi * se, se * b.se);
    }
    Frac operator+= (const Frac & b) {
        fi = fi * b.se + b.fi * se;
        se *= b.se;
        upd();
        return *this;
    }
    Frac operator- (const Frac & b) const {
        return Frac(fi * b.se - b.fi * se, se * b.se);
    }
    Frac operator- () const {
        return Frac(-fi, se);
    }
    Frac operator-= (const Frac & b) {
        fi = fi * b.se - b.fi * se;
        se *= b.se;
        upd();
        return *this;
    }
    Frac operator* (const Frac & b) const {
        return Frac(fi * b.fi, se * b.se);
    }
    Frac operator*= (const Frac & b) {
        fi *= b.fi, se *= b.se;
        upd();
        return *this;
    }
    Frac operator/ (const Frac & b) const {
        assert(b);
        return Frac(fi * b.se, se * b.fi);
    }
    Frac operator/= (const Frac & b) {
        assert(b);
        fi *= b.se, se *= b.fi;
        upd();
        return *this;
    }
    bool operator==(const Frac & b) const {
        return fi == b.fi && se == b.se;
    }
    bool operator!=(const Frac & b) const {
        return fi != b.fi || se != b.se;
    }
    bool operator==(const T & x) const {
        return fi == x && se == T(1);
    }
    bool operator!=(const T & x) const {
        return fi != x || se != T(1);
    }
    bool operator<(const Frac & b) const {
        return fi * b.se < se * b.fi;
    }
    bool operator<=(const Frac & b) const {
        return fi * b.se <= se * b.fi;
    }
    bool operator>(const Frac & b) const {
        return fi * b.se > se * b.fi;
    }
    bool operator>=(const Frac & b) const {
        return fi * b.se >= se * b.fi;
    }
    bool operator<(const T & x) const {
        return fi < se * x;
    }
    bool operator<=(const T & x) const {
        return fi <= se * x;
    }
    bool operator>(const T & x) const {
        return fi > se * x;
    }
    bool operator>=(const T & x) const {
        return fi >= se * x;
    }
    operator bool() const {
        return fi;
    }

    T fi, se;
};

template<typename T>
std::ostream& operator<< (std::ostream& out, const Frac<T> & a) {
    if (!a.fi) {
        out << 0;
        return out;
    }
    if (a.se == T(1)) {
        out << a.fi;
        return out;
    }
    out << a.fi << "/" << a.se;
    return out;
}

void swap(long long int & a, long long int & b) {
    a ^= b;
    b ^= a;
    a ^= b;
}

template<typename T>
void swap(Frac<T> & a, Frac<T> & b) {
    swap(a.fi, b.fi);
    swap(a.se, b.se);
    return;
}

template<typename T>
Frac<T> sqrt(const Frac<T> & a) {
    return Frac<T>(sqrt(a.fi), sqrt(a.se));
}


template<typename T>
Frac<T> abs(const Frac<T> & a) {
    return Frac<T>(abs(a.fi), a.se);
}

