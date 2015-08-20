#include <bits/stdc++.h>
using namespace std;

typedef vector< vector<double> > Matriz;
const double ERROR = 1e-9;

bool Igual(double a, double b) {
    return fabs(a - b) < ERROR;
}

void EliminaGaussiana(Matriz& m) {
    for (int i = 0; i < m.size(); ++i) {
        int fila_mayor = i;
        for (int j = i + 1; j < m.size(); ++j)
            if (fabs(m[fila_mayor][i]) <
                fabs(m[j][i])) fila_mayor = j;
        swap(m[i], m[fila_mayor]);
        if (Igual(m[i][i], 0)) continue;
        for (int j = m[i].size() - 1; j >= i; --j)
            m[i][j] /= m[i][i];
        for (int j = 0; j < m.size(); ++j) {
            if (i == j || Igual(m[j][i], 0)) continue;
            double ratio = m[j][i] / m[i][i];
            for (int k = i; k < m[j].size(); ++k)
                m[j][k] -= m[i][k] * ratio;
        }
    }
}

int main() {
    return 0;
}