#include "matrices.h"
#include "isomorphism.h"
using namespace std;

vec<Matrix> get_K11112(Frac n, Frac m) { 
    Matrix A1 = {{n, n}, 
                 {n, m}};
    
    Matrix A2 = {{n, n}, 
                 {n, n + n - m}};
    
    Matrix A3 = {{-n, n}, 
                 {-n, m}};
    
    Matrix A4 = {{-n, -n}, 
                 {n, -m}};

    Matrix X1 = {{-n, -n}, 
                 {-m - m - n, -m}};

    Matrix X2 = {{-n, m + m - n}, 
                 {-n, m}};
    
    return {A1, A2, A3, A4, X1, X2};
}


int main() {
    // int maxn = 50, maxm = 5;
    // int maxn = 100, maxm = 2;
    int maxn = 20, maxm = 20;
    bool ok = true;
    for (int n1 = 1; n1 < maxn; n1++) {
        for (int m1 = 1; m1 < maxm; m1++) {
            if (m1 == n1)
                continue;
            Frac f1 = make_frac(n1, 1);
            Frac f2 = make_frac(m1, 1);
            vec<Matrix> k1 = get_K11112(f1, f2);
            k1.pop_back();
            for (int n2 = n1 + 1; n2 < maxn; n2++) {
                for (int m2 = 1; m2 < maxm; m2++) {
                    if (n2 == m2)
                        continue;
                    if (make_frac(n1, m1) == make_frac(n2, m2) || make_frac(n1, m1) == make_frac(m2, n2))
                        continue;
                    Frac f3 = make_frac(n2, 1);
                    Frac f4 = make_frac(m2, 1);
                    vec<Matrix> k2 = get_K11112(f3, f4);
                    k2.pop_back();
                    
                    if (isomorphic_subgraphs_2_by_2(k1, k2)) {
                        cout << "Not coresopnds with hypothesis" << endl;
                        cout << n1 << " " << m1 << ", " << n2 << " " << m2 << endl; 
                        ok = false;
                    }
                }
            }
            if (ok)
                cout << "ok for " << n1 << " " << m1 << endl;
        }
    }
    if (ok)
        cout << "Evrything ok" << endl;
}
