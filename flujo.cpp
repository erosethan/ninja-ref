#include <bits/stdc++.h>
using namespace std;

// Definiciones iniciales
typedef pair<int, int> Arista;

const int INF = 1 << 30;
const int MAXN = 100;

vector<int> grafo[MAXN];

// Maximo emparejamiento en grafo bipartito. Nodos indexados
// de 0 a n - 1. Recibe dos parametros: Un vector con los indices
// del conjunto izquierdo y otro con los indices del derecho.

int pareja[MAXN];
bool visitado[MAXN];

int CaminoIncremental(int u) {
    visitado[u] = true;
    for (int i = 0; i < grafo[i].size(); ++i) {
        int v = grafo[u][i];
        if (pareja[v] == -1) 
            return pareja[v] = u;
        if (visitado[pareja[v]]) continue;
        if (v == CaminoIncremental(
            pareja[v])) return u;
    }
    return -1;
}

vector<Arista> MaxEmparejamientoBipartito(
    const vector<int>& l, const vector<int>& r) {
    int n = l.size() + r.size();
    fill(pareja, pareja + n, -1);
    for (int i = 0; i < l.size(); ++i) {
        fill(visitado, visitado + n, false);
        CaminoIncremental(l[i]);
    }
    vector<Arista> parejas;
    for (int i = 0; i < r.size(); ++i)
        if (pareja[r[i]] != -1) parejas.push_back(
            Arista(pareja[r[i]], r[i]));
    return parejas;
}

// Emparejamiento de costo maximo en grafo bipartito ponderado.
// Nodos con indice del 0 al n - 1. Recibe los mismos parametros
// que MaxEmparejamientoBipartito. Utiliza algoritmo hungaro.
// Cuidado: Se ocupa una matriz de costo entre nodos.

int slack[MAXN];
int retorno[MAXN];
int etiqueta[MAXN];
int costo[MAXN][MAXN];

vector<Arista> EmparejaCostoMaxBipartito(
    const vector<int>& l, const vector<int>& r) {
    int n = l.size() + r.size();
    fill(pareja, pareja + n, -1);
    fill(etiqueta, etiqueta + n, 0);
    for (int i = 0; i < l.size(); ++i)
        for (int j = 0; j < r.size(); ++j)
            etiqueta[l[i]] = max(etiqueta[l[i]],
                costo[l[i]][r[j]]);
            
    for (int i = 0; i < l.size(); ++i) {
        for (int j = 0; j < r.size(); ++j)
            slack[r[j]] = -costo[l[i]][r[j]] +
                etiqueta[l[i]] + etiqueta[r[j]];
        fill(visitado, visitado + n, false);
        fill(retorno, retorno + n, l[i]);
        visitado[l[i]] = true;
        
        bool emparejado = false;
        for (int j = 0; !emparejado; ++j) {
            int t = 0; for (; t < r.size(); ++t) {
                if (visitado[r[t]]) continue;
                if (!slack[r[t]]) break;
            }
            if (t < r.size()) {
                visitado[t = r[t]] = true;
                if (pareja[t] == -1) {
                    emparejado = true;
                    for (int p; ; t = p) {
                        pareja[t] = retorno[t];
                        p = pareja[retorno[t]];
                        pareja[retorno[t]] = t;
                        if (retorno[t] == l[i]) break;
                    }
                } else {
                    visitado[t = pareja[t]] = true;
                    for (int k = 0; k < r.size(); ++k) {
                        int new_slack = etiqueta[t] +
                            etiqueta[r[k]] - costo[t][r[k]];
                        if (new_slack < slack[r[k]]) {
                            slack[r[k]] = new_slack;
                            retorno[r[k]] = t;
                        }
                    }
                }
            } else {
                int d = INF;
                for (int k = 0; k < r.size(); ++k)
                    if (slack[r[k]]) d = min(d, slack[r[k]]);
                for (int k = 0; k < l.size(); ++k)
                    if (visitado[l[k]]) etiqueta[l[k]] -= d;
                for (int k = 0; k < r.size(); ++k)
                    if (!visitado[r[k]]) slack[r[k]] -= d;
                    else etiqueta[r[k]] += d;
            }
        }
    }
    vector<Arista> parejas;
    for (int i = 0; i < l.size(); ++i)
        parejas.push_back(Arista(l[i],pareja[l[i]]));
    return parejas;
}

int main() {
    return 0;
}
