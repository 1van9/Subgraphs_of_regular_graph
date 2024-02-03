#include <iostream>
#include <random>
#include <vector>
using namespace std;

#define vec vector
#define fi first
#define se second
#define mfrac make_pair


long long gcd(long long a, long long b) {
    if (b == 0) return a;
    return gcd(b, a % b);
}


namespace {
    
using Frac = pair<long long, long long>;
using Matrix = vec<vec<Frac>>;

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
Frac det(const Matrix & a) {
    return a[0][0] * a[1][1] - a[0][1] * a[1][0];
}
bool operator- (const Matrix & a, const Matrix & b) {
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

};

vector<Matrix> small;

void gen_small_matrix(int i, int j, Matrix &curr, vec<Frac>& nums) {
    if (j == 2) {   //shift
        i++;
        j = 0;
    }
    if (i == 2) {
        if (GL(curr)) {
            small.push_back(curr);
        }
        return;
    }
    for (auto num : nums) {
        curr[i][j] = num;
        gen_small_matrix(i, j + 1, curr, nums);
    }
}

bool is_tl = false;
int frequency = 1000;
double tl = 10;
double start;

bool is_tle() {
    static int cnt = 0;
    cnt++;
    if (cnt == frequency) {
        cnt = 0;
        is_tl |= (clock() - start > (tl * CLOCKS_PER_SEC));
    }
    return is_tl;
}

Frac gen_frac(mt19937 & gen) {
    int a = gen() % 10000, b = gen() % 10000;
    while (b == 0) b = gen() % 10000;
    return {a, b};
}

Matrix gen_matrix1(mt19937& gen) {
    Matrix a(2, vec<Frac>(2, {0, 1}));
    while (det(a).fi == 0) {
        a[0][0] = gen_frac(gen);
        a[0][1] = gen_frac(gen);
        a[1][0] = gen_frac(gen);
        a[1][1] = gen_frac(gen);
    }
    return a;
}

Matrix gen_matrix2(mt19937 & gen) {
    return small[gen() % small.size()];
}

int number_of_steps = 100'000, desperation = 500;
int big_number_of_steps = 200'000, big_desperation = 5000;
bool mode = true; // all small matrix


void gen_clique(vector<Matrix> & clique, mt19937 & gen, int last_ind = -1) {
    if (clique.size() < 4) {
        int num_steps;
        if (mode)
            num_steps = (last_ind == -1) ? small.size() : last_ind;
        else
            num_steps = number_of_steps;

        while (num_steps--) {
            Matrix new_matrix;
            if (mode)
                new_matrix = small[num_steps];
            else
                new_matrix = (num_steps > desperation) ? gen_matrix2(gen) : gen_matrix1(gen); 
            bool ok = true;
            for (auto el : clique) {
                if (!(el - new_matrix)) {
                    ok = false; 
                    break;
                }
            }
            if (ok) {
                clique.push_back(new_matrix);
                gen_clique(clique, gen, num_steps);
                clique.pop_back();
            }
        }
    } else {
        cout << "minimal success" << endl;
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
        // this linear equasions always independent
        Gause(linear_equations);
  
        if (linear_equations[2][0].fi == 0 && linear_equations[2][1].fi == 0 &&
            linear_equations[2][2].fi == 0 && linear_equations[2][3].fi == 0) {    
            // 1 * 0 * *     
            // 0 0 1 * *   or something like this
            // 0 0 0 0 *
            cout << "No more than one solution" << endl;
            return;
        }

        // 1 0 3 0 *
        // 0 1 1 0 *
        // 0 0 0 1 *
        int val = 3;
        for (int i = 0; i < 3; i++) {
            if (linear_equations[i][i].fi == 0) {
                val = i;
                break;
            }
        }
        if (val == 0 || val == 1) {
            // 0 1 0 0 *
            // 0 0 1 0 *     val = 0
            // 0 0 0 1 *
            
            // 1 2 0 0 *
            // 0 0 1 0 *     val = 1
            // 0 0 0 1 *
            cout << "No more than one solution" << endl;
            return;
        }
        if (val == 2) {
            // 1 0 -a 0 -c 
            // 0 1 -b 0 -d   t = e; z; y = d + bz; x = c + az;
            // 0 0  0 1 -e

            // (az + c)e - (d + bz)z + clique[0][1][1](c + az) - clique[0][1][0](d + bz) - clique[0][0][1]z + clique[0][0][0]e + det(clique[0])
            // -b * z^2 + (ae - d + a clique[0][1][1] - b clique[0][1][0] - clique[0][0][1])t + ...
            Frac a = mfrac(0, 1)-linear_equations[0][3];
            Frac b = mfrac(0, 1)-linear_equations[1][3];
            Frac c = mfrac(0, 1)-linear_equations[0][4];
            Frac d = mfrac(0, 1)-linear_equations[1][4];
            Frac e = mfrac(0, 1)-linear_equations[2][4];
            Frac A = mfrac(0, 1) - d;
            Frac B = a * e  - d + a * clique[0][1][1] - b * clique[0][1][0] - clique[0][0][1];
            Frac C = c * e + c * clique[0][1][1] - clique[0][1][0] * d - clique[0][0][0] * e + det(clique[0]);
            // Az^2 + Bz + C
            Frac D = B * B - (mfrac(4, 1) * A * C);
            if (A.fi == 0 || D.fi == 0) {
                cout << "No more than one solution" << endl;
                return;
            }
            cout << "Success" << endl;
            // cout << "equations : " << endl;
            // cout << linear_equations << endl;
            // cout << A << " z^2 + " << B << " z + " << C << endl;

            long long x = sqrt(D.fi), y = sqrt(D.se);
            if (x * x == D.fi && y * y == D.se) {
                Frac D_ = mfrac(x, y);
                Frac z1 = (D_ - B) / (mfrac(2, 1) * A);
                Frac z2 = (mfrac(0, 1) - D_ - B) / (mfrac(2, 1) * A);
                Frac x1 = c + z1 * a, y1 = d + z1 * b, t1 = e;
                Frac x2 = c + z2 * a, y2 = d + z2 * b, t2 = e;
                if ((x1 * t1 - y1 * z1).fi == 0 || (x2 * t2 - y2 * z2).fi == 0) {
                    cout << "Not now" << endl;
                    return;
                } else {
                    cout << "---Got it!---" << endl;
                    cout << clique[0] << endl << clique[1] << endl;
                    cout << clique[2] << endl << clique[3] << endl;
                    cout << "-----and-----" << endl;
                    cout << x1 << " " << y1 << endl;
                    cout << z1 << " " << t1 << endl;
                    cout << endl;
                    cout << x2 << " " << y2 << endl;
                    cout << z2 << " " << t2 << endl;
                    cout << "-------------";
                    cout << endl;
                    exit(0);
                    return;
                }
            } else {
                cout << "Maybe" << endl;
                return;
            }
        }
        // 1 0 0 -a -d
        // 0 1 0 -b -e   t; z = f + ct; y = e + bt; x = d + at;  
        // 0 0 1 -c -f
        
        // (at + d)t - (e + bt)(f + ct) + clique[0][1][1](d + at) - clique[0][1][0](e + bt) - clique[0][0][1](f + ct) + clique[0][0][0]t + det(clique[0])
        // (a - bc) * t^2 + (d - ec - bf + a clique[0][1][1] - b clique[0][1][0] - c clique[0][0][1] + clique[0][0][0])t + ...
        Frac a = mfrac(0, 1)-linear_equations[0][3];
        Frac b = mfrac(0, 1)-linear_equations[1][3];
        Frac c = mfrac(0, 1)-linear_equations[2][3];
        Frac d = mfrac(0, 1)-linear_equations[0][4];
        Frac e = mfrac(0, 1)-linear_equations[1][4];
        Frac f = mfrac(0, 1)-linear_equations[2][4];
        Frac A = a - (c*b);
        Frac B = d - (e * c) - (b * f) + (a * clique[0][1][1]) - (b * clique[0][1][0]) - (c * clique[0][0][1]) + clique[0][0][0];
        Frac C = mfrac(0, 1)-(e * f) + (clique[0][1][1] * d) - (clique[0][1][0] * e) - (clique[0][0][1] * f) + det(clique[0]);
        //At^2 + Bt + C
        Frac D = B * B - (mfrac(4, 1) * A * C);
        if (A.fi == 0 || D.fi == 0) {
            cout << "No more than one solution" << endl;
            return;
        }
        cout << "Success" << endl;
        // cout << "equations : " << endl;
        // cout << linear_equations << endl;
        // cout << A << " x^2 + " << B << " x + " << C << endl;

        long long x = sqrt(D.fi), y = sqrt(D.se);
        if (x * x == D.fi && y * y == D.se) {
            Frac D_ = mfrac(x, y);
            Frac t1 = (D_ - B) / (mfrac(2, 1) * A);
            Frac t2 = (mfrac(0, 1) - D_ - B) / (mfrac(2, 1) * A);
            Frac z1 = f + c * t1, y1 = e + b * t1, x1 = d + a * t1;
            Frac z2 = f + c * t2, y2 = e + b * t2, x2 = d + a * t2;
            if ((x1 * t1 - y1 * z1).fi == 0 || (x2 * t2 - y2 * z2).fi == 0) {
                cout << "No more than one solution" << endl;
                return;
            } else {
                cout << "---Got it!---" << endl;
                cout << clique[0] << endl << clique[1] << endl;
                cout << clique[2] << endl << clique[3] << endl;
                cout << "-----and-----" << endl;
                cout << x1 << " " << y1 << endl;
                cout << z1 << " " << t1 << endl;
                cout << endl;
                cout << x2 << " " << y2 << endl;
                cout << z2 << " " << t2 << endl;
                cout << "-------------" << endl;
                exit(0);
                return;
            }
        } else {
            cout << "Maybe" << endl;
            return;
        }
        return;

    }
}

int main() {
    mt19937 gen(1786);
    vector<Frac> small_nums;
    int C = 2;
    for (int i = -C*4; i <= C*4; i++) {
        small_nums.push_back(upd(mfrac(i, 4)));
    }
    for (int i = -C*3; i <= C*3; i++) {
        if (i % 3 != 0)
            small_nums.push_back(upd(mfrac(i, 3)));
    }
    // for (int i = -C*5; i <= C*5; i++) {
    //     if (i % 5 != 0)
    //         small_nums.push_back(mfrac(i, 5));
    // }
    Matrix curr(2, vec<Frac>(2));
    gen_small_matrix(0, 0, curr, small_nums);
    int n = small.size();
    cout << n << endl;
    
    start = clock();
    vector<Matrix> clique;
    gen_clique(clique, gen);
    //result:

    // Matrix A1 = {{{5, 1}, {5, 1}}, 
    //              {{5, 1}, {4, 1}}};
    
    // Matrix A2 = {{{5, 1}, {5, 1}}, 
    //              {{5, 1}, {6, 1}}};
    
    // Matrix A3 = {{{-5, 1}, {5, 1}}, 
    //              {{-5, 1}, {4, 1}}};
    
    // Matrix A4 = {{{-5, 1}, {-5, 1}}, 
    //              {{5, 1}, {-4, 1}}};

    // Matrix X1 = {{{-5, 1}, {-5, 1}}, 
    //              {{-13, 1}, {-4, 1}}};
    
    // Matrix X2 = {{{-5, 1}, {3, 1}}, 
    //              {{-5, 1}, {4, 1}}};
}
