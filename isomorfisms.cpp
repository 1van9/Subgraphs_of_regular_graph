#include "matrices.h"
#include "isomorfisms.h"
using namespace std;

int main() {
    
    vec<Matrix> Akbari = {
        {{{1, 1}, {0, 1}},
         {{0, 1}, {1, 1}}},
        
        {{{0, 1}, {1, 1}},
         {{1, 1}, {0, 1}}},
        
        {{{0, 1}, {-1, 1}},
         {{-1, 1}, {0, 1}}},

        {{{0, 1}, {1, 1}},
         {{-1, 1}, {-2, 1}}},
        
        {{{0, 1}, {-1, 1}},
         {{1, 1}, {-2, 1}}}
    };

    vec<Matrix> new_sample = {
        {{{2, 1}, {2, 1}},
         {{2, 1}, {1, 1}}},
        
        {{{2, 1}, {2, 1}},
         {{2, 1}, {3, 1}}},
        
        {{{-2, 1}, {2, 1}},
         {{-2, 1}, {1, 1}}},

        {{{-2, 1}, {-2, 1}},
         {{2, 1}, {-1, 1}}},
        
        {{{-2, 1}, {-2, 1}},
         {{-4, 1}, {-1, 1}}}
    };

    if (isomorphic_subgraphs_2_by_2(Akbari, new_sample)) 
        cout << "isomorphic subgraphs" << endl;
    else
        cout << "not isomorphic subgraphs" << endl;
    
}