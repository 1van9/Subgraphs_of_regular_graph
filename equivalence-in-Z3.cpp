#include "any-matrix.h"

using namespace std;

void PrintClique(const Clique3 & a) {
    for (int i = 0; i < a[0].size(); i++) {
        for (auto el : a) {
            for (int j = 0; j < el[i].size(); j++) cout << el[i][j] << " ";
            cout << "   ";
        }
        cout << endl;
    }
    cout << endl;
}

vector<Mat3> Gl;

vector<Clique3> all_max_cliques;
int max_ans = 5;


void GetCliques(Clique3 & curr, int i) {
    if (curr.size() == max_ans) {
        all_max_cliques.push_back(curr);
        return;
    }
    for (; i < Gl.size(); i++) {
        bool connect = true;
        for (auto el : curr) {
            if (!Gl[i].Connect(el)) {
                connect = false;
            }
        }
        if (connect) {
            curr.push_back(Gl[i]);
            GetCliques(curr, i + 1);
            curr.pop_back();
        }
    }
}




int main() {
    vec<Z3> all_nums({0, 1, 2});
    Gl = get_matrixes(all_nums, 2);    
    cout << "Number of matrix : " <<  Gl.size() << endl;
    
    Clique3 curr;
    GetCliques(curr, 0);
    cout << "Number of cliques : " << all_max_cliques.size() << endl;
  
    for (int i = 1; i < all_max_cliques.size(); i++) {
        if (!isomorphic_subgraphs_2_by_2(all_max_cliques[0], all_max_cliques[i])) {
            cout << "Not equivalent cliques:" << endl;
            PrintClique(all_max_cliques[0]);
            PrintClique(all_max_cliques[i]);
            return 0;
        }
    }
    cout << "All cliques are equivalent" << endl;
    
    return 0;
}
