#pragma ones
// header for working with integer numbers

#include<iostream>
#include<vector>
#include<algorithm>
#define vec std::vector

inline vec<long long> nm(long long x, long long sys) {
    vec<long long> res;
    while (x > 0) {
        res.push_back(x % sys);
        x /= sys;
    }
    return res;
}

inline vec<long long> sm(const vec<long long> & x, 
                         const vec<long long> & y, 
                         long long sys) {
    vec<long long> res;
    size_t i = 0;
    long long shift = 0;
    while (true) {
        long long currd = shift;
        if (i < x.size()) 
            currd += x[i];
        if (i < y.size()) 
            currd += y[i];
        if (i >= std::max(x.size(), y.size()) && currd == 0) 
            break;
        res.push_back(currd - sys * (currd >= sys));
        shift = (currd >= sys);
        i++;
    }
    return res;
}

inline vec<long long> df(const vec<long long> & x, 
                         const vec<long long> & y, 
                         long long sys) {
    vec<long long> res;
    size_t i = 0;
    long long shift = 0;
    while (true) {
        long long currd = shift;
        if (i < x.size()) 
            currd += x[i];
        if (i < y.size()) 
            currd -= y[i];
        if (i >= std::max(x.size(), y.size()) && currd <= 0) 
            break;
        shift = -1 * (currd < 0);
        res.push_back(currd + sys * (currd < 0));
        i++;
    }
    return res;
}

inline bool ls(const vec<long long> & x, 
               const vec<long long> & y, 
               long long sys) {
    size_t i = 0;
    long long shift = 0;
    while (true) {
        long long currd = shift;
        if (i < x.size()) 
            currd += x[i];
        if (i < y.size()) 
            currd -= y[i];
        if (i >= std::max(x.size(), y.size()) && currd <= 0) 
            break;
        shift = -1 * (currd < 0);
        i++;
    }
    return shift < 0;
}

inline vec<long long> ml(const vec<long long> & x, 
                         const vec<long long> & y, 
                         long long sys) {
    vec<long long> res(x.size() + y.size());
    for (size_t i = 0; i < x.size(); i++) {
        for (size_t j = 0; j < y.size(); j++) {
            res[i + j] += x[i] * y[j];
        }
    }
    long long shift = 0;
    for (auto & d : res) {
        d += shift;
        shift = d / sys;
        d %= sys;
    }
    return res;
}

inline std::pair<vec<long long>, vec<long long>> dv(const vec<long long> & x, 
                                               const vec<long long> & y, 
                                               long long sys) {
    vec<long long> res;
    if (x.size() < y.size()) {
        return make_pair(res, x);
    }
    vec<long long> rem = x;
    res.resize(x.size() - y.size() + 1);
    for (int i = x.size() - y.size(); i >= 0; i--) {
        vec<long long> num(i + 1);
        long long l = 0, r = sys;
        while (r - l > 1) {
            long long mid = (l + r) / 2;
            num[i] = mid;
            if (ls(rem, ml(y, num, sys), sys)) {
                r = mid;
            } else {
                l = mid;
            }
        }
        res[i] = l;
        num[i] = l;
        rem = df(rem, ml(y, num, sys), sys);
    }
    return std::make_pair(res, rem);
}

vec<long long> conv_seg(size_t l, size_t r, long long sys1, 
                        const vec<long long> & dig, 
                        long long sys2, vec<long long>& pw) {
    if (l + 1 == r) {
        pw = nm(sys1, sys2);
        return nm(dig[l], sys2);
    }
    int mid = (l + r) / 2;
    vec<long long> pw1, pw2;
    vec<long long> x1 = conv_seg(l, mid, sys1, dig, sys2, pw1);
    vec<long long> x2 = conv_seg(mid, r, sys1, dig, sys2, pw2);
    pw = ml(pw1, pw2, sys2);
    return sm(ml(x2, pw1, sys2), x1, sys2);
}

inline vec<long long> conv(long long sys1, const vec<long long>& dig, 
                           long long sys2) {
    if (sys1 == sys2)
        return dig;
    vec<long long> tmp;
    return conv_seg(0, dig.size(), sys1, dig, sys2, tmp);
}

struct Integer {
    void update() {
        while (!digits.empty() && !digits.back()) 
            digits.pop_back();
        if (digits.empty()) 
            sign = 1;
    }
    Integer() {}
    Integer(long long x) {
        if (x < 0) {
            sign = -1;
            x = -x;
        }
        digits = nm(x, system);
        update();
    } 
    Integer(const vec<long long> & dig, long long sg) : digits(dig), sign(sg) {
        update();
    }
    Integer(const std::string & s) {
        size_t start = 0;
        if (s[start] == '-')
            sign = -1, start++;
        if (s[start] == '+')
            start++;
        vec<long long> num(s.size() - start);
        for (size_t i = 0; i < num.size(); i++)
            num[i] = s[s.size() - i - 1] - '0';
        digits = conv(10, num, system);
        update();
    }
    Integer(const Integer & n) : digits(n.digits), sign(n.sign) {}
    Integer operator = (const Integer & n) {
        digits = n.digits;
        sign = n.sign;
        return *this;
    }
    Integer operator + (const Integer & n) const {
        if (sign == n.sign)
            return Integer(sm(digits, n.digits, system), sign);
        if (ls(n.digits, digits, system))
            return Integer(df(digits, n.digits, system), sign);
        return Integer(df(n.digits, digits, system), n.sign);
    }
    Integer & operator += (const Integer & n) {
        *this = (*this) + n;
        return *this;
    }
    Integer operator - (const Integer & n) const {
        return (*this) + Integer(n.digits, -n.sign);
    }
    Integer & operator -= (const Integer & n) {
        (*this) = (*this) - n;
        return *this;
    }
    Integer operator * (const Integer & n) const {
        return Integer(ml(n.digits, digits, system), sign * n.sign);
    }
    Integer & operator *= (const Integer & n) {
        *this = (*this) * n;
        return *this;
    }
    Integer operator / (const Integer & n) const {
        return Integer(dv(digits, n.digits, system).first, sign * n.sign);
    }
    Integer & operator /= (const Integer & n) {
        *this = (*this) / n;
        return *this;
    }
    Integer operator % (const Integer & n) const {
        vec<long long> res = dv(digits, n.digits, system).second;
        if (sign == -1 && res.size())
            return Integer(df(n.digits, res, system), 1);
        return Integer(res, 1);
    }
    Integer & operator %= (const Integer & n) {
        (*this) = (*this) % n;
        return *this;
    }
    Integer operator ^ (Integer x) const {
        if (x.digits.empty())
            return 1;
        if (x.digits[0] & 1)
            return ((*this) ^ (x - Integer(1))) * (*this);
        return ((*this) * (*this)) ^ (x / Integer(2));
    }
    bool operator == (const Integer & n) const {
        return digits == n.digits && sign == n.sign;
    }
    bool operator != (const Integer & n) const {
        return !(digits == n.digits);
    }
    bool operator < (const Integer & n) const {
        if (sign == 1) {
            if (n.sign == -1) 
                return false;
            return ls(digits, n.digits, system);
        }
        if (n.sign == 1)
            return true;
        return ls(n.digits, digits, system);
    }
    bool operator <= (const Integer & n) const {
        return !(n < (*this));
    }
    bool operator > (const Integer & n) const {
        return n < (*this);
    }
    bool operator >= (const Integer & n) const {
        return !((*this) < n);
    }
    operator bool() const {
        return digits.size();
    }
    
    vec<long long> digits;
    int sign = 1;
    const long long system = 1 << 16;
};

std::istream & operator>>(std::istream & in, Integer & n) {
    std::string s;
    in >> s;
    n = Integer(s);
    return in;
}

std::ostream & operator<<(std::ostream & out, const Integer & n) {
    if (n.digits.empty()) {
        out << 0;
        return out;
    }
    if (n.sign == -1) out << '-';
    vec<long long> dig = conv(n.system, n.digits, 10);
    reverse(dig.begin(), dig.end());
    size_t i = 0;
    while (i < dig.size() && dig[i] == 0)
        i++;
    for (; i < dig.size(); i++)
        out << dig[i];
    return out;
}

Integer abs(const Integer & n) {
    return Integer(n.digits, 1);
}

void swap(Integer & n1, Integer & n2) {
    std::swap(n1.digits, n2.digits);
    std::swap(n1.sign, n2.sign);
}

Integer min(const Integer & a, const Integer & b) {
    if (a < b)
        return a;
    return b;
}

Integer max(const Integer & a, const Integer & b) {
    if (a > b)
        return a;
    return b;
}
