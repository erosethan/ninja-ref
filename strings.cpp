#include <bits/stdc++.h>
using namespace std;

void BucketSort(vector<int>& sa,
    const vector<int>& rank, int ranks) {

    vector<int> bucket(ranks, 0);
    vector<int> tmp_sa(sa.size());
    for (int i = 0; i < sa.size(); ++i)
        ++bucket[rank[sa[i]]];

    for (int i = 0, sum = 0; i < ranks; ++i) {
        swap(bucket[i], sum); sum += bucket[i];
    } for (int i = 0; i < sa.size(); ++i)
        tmp_sa[bucket[rank[sa[i]]]++] = sa[i];
    swap(sa, tmp_sa);
}

// Recuerden poner '$' al final de la cadena.

vector<int> SuffixArray(const string& str) {
    int ranks = 255;
    vector<int> sa(str.size());
    vector<int> nrank(str.size());
    vector<int> rank(str.size(), 0);
    vector<int> tmp_rank(str.size());
    for (int i = 0; i < str.size(); ++i)
        nrank[i] = str[i], sa[i] = i;

    for (int p = 0; true; ++p) {
        BucketSort(sa, nrank, ranks + 1);
        BucketSort(sa,  rank, ranks + 1);

        tmp_rank[0] = ranks = 0;
        for (int i = 1; i < str.size(); ++i)
            if (rank[sa[i]] != rank[sa[i - 1]] ||
                nrank[sa[i]] != nrank[sa[i - 1]])
                 tmp_rank[i] = ++ranks;
            else tmp_rank[i] = ranks;

        if (ranks + 1 == str.size()) break;

        for (int i = 1; i <= 1 << p; ++i)
            nrank[str.size() - i] = 0;
        for (int i = 0; i < str.size(); ++i) {
            int prv = sa[i] - (1 << p);
            if (prv >= 0) nrank[prv] = tmp_rank[i];
            rank[sa[i]] = tmp_rank[i];
        }
    }
    return sa;
}

vector<int> LCP(const string& str,
    const vector<int>& sa) {

    vector<int> lcp(str.size());
    vector<int> plcp(str.size());
    vector<int> phi(str.size(), -1);
    for(int i = 1; i < str.size(); ++i)
        phi[sa[i]] = sa[i - 1];
    
    int len = 0;
    for(int i = 0; i < str.size(); ++i) {
        if (phi[i] == -1) continue;
        for (; str[phi[i] + len] ==
               str[i + len]; ++len) {}
        plcp[i] = len; len = max(len-1, 0);
    }
    for (int i = 0; i < str.size(); ++i)
        lcp[i] = plcp[sa[i]];
    return lcp;
}

int main() {
    return 0;
}
