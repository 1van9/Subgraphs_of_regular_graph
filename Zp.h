#include<iostream>


template<int P>
struct Zp {
    Zp(): x(0) {}
    Zp(int i) : x(i % P) {
        if (x < 0)
            x += P;
    }
    Zp(const Zp & a) : x(a.x) {}
    Zp operator=(const Zp & a) {
        x = a.x;
        return a;
    }
    Zp operator + (const Zp& other) const {
        return Zp(x + other.x);
    }
    Zp operator - (const Zp& other) const {
        return Zp(x - other.x);
    }
    Zp operator += (const Zp& other) {
        x += other.x;
        if (x >= P)
            x -= P;
        return *this;
    }
    Zp operator -= (const Zp& other) {
        x -= other.x;
        if (x < 0)
            x += P;
        return *this;
    }
    Zp operator * (const Zp& other) const {
        return Zp(x * other.x);
    }
    Zp operator *= (const Zp& other) {
        x *= other.x;
        x %= P;
        return *this;
    }
    Zp operator^(int k) const {
        if (k == 0)
            return Zp(1);
        Zp val = ((*this) * (*this))^(k >> 1);
        if (k & 1)
            return val * (*this);
        else
            return val; 
    }
    Zp operator / (const Zp& other) const {
        return (*this) * (other^(P - 2));
    }
    Zp operator /= (const Zp& other) {
        *this = (*this / other);
        return (*this);
    }
    operator bool() const {
        return x;
    }
    
    int x;
};

template<int P>
std::istream& operator >>(std::istream & in, Zp<P>& a) {
    in >> a.x;
    a.x %= P;
    if (a.x < 0)
        a.x += P;
    return in;
} 

template<int P>
std::ostream& operator << (std::ostream & out, const Zp<P>& a) {
    out << a.x;
    return out;
}

template<int P>
void swap(Zp<P> & a, Zp<P> & b) {
    a.x ^= b.x;
    b.x ^= a.x;
    a.x ^= b.x;
}
