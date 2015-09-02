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
    for (int i = 0; i < grafo[u].size(); ++i)
        if (pareja[grafo[u][i]] == -1) 
            return pareja[grafo[u][i]] = u;
    for (int i = 0; i < grafo[u].size(); ++i) {
        int v = grafo[u][i];
        if (visitado[pareja[v]]) continue;
        if (CaminoIncremental(pareja[v]) != -1)
            return pareja[v] = u;
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
    // Si l.size() != r.size() KABOOM!
    assert(l.size() == r.size());

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
        parejas.push_back(Arista(
            l[i], pareja[l[i]]));
    return parejas;
}

// Flujo maximo en un grafo dirigido mediante Edmonds-Karp O(VE^2).
// Los nodos están indexados del 0 al n - 1. Para cada arista (u, v)
// en el grafo debe existir la arista residual (v, u). Se debe cumplir
// cap[u][v] = cap[v][u], flujo[u][v] = 0 y flujo[v][u] = cap[u][v].

int padre[MAXN];
int cap[MAXN][MAXN];
int flujo[MAXN][MAXN]; 

int ActualizarFlujo(int u, int f) {
    int p = padre[u];
    if (p == u) return f;
    f = ActualizarFlujo(p, min(
        f, cap[p][u] - flujo[p][u]));
    flujo[p][u] += f;
    flujo[u][p] -= f;
    return f;
}

int AumentarFlujo(int s, int t, int n) {
    fill(padre, padre + n, -1);
    queue<int> q; q.push(s); padre[s] = s;
    while (!q.empty()) {
        int u = q.front();
        q.pop(); if (u == t) break;
        for (int i = 0; i < grafo[u].size(); ++i) {
            int v = grafo[u][i];
            if (flujo[u][v] == cap[u][v] ||
                padre[v] != -1) continue;
            padre[v] = u, q.push(v);
        }
    }
    if (padre[t] == -1) return 0;
    return ActualizarFlujo(t, INF);
}

int EdmondsKarp(int s, int t, int n) {
    int flujo_maximo = 0, f;
    while (f = AumentarFlujo(s, t, n))
        flujo_maximo += f;
    return flujo_maximo;
}

// Flujo maximo en un grafo dirigido mediante Dinic O(V^2E).
// Los nodos están indexados del 0 al n - 1. Para cada arista (u, v)
// en el grafo debe existir la arista residual (v, u). Se debe cumplir
// cap[u][v] = cap[v][u], flujo[u][v] = 0 y flujo[v][u] = cap[u][v].

int dist[MAXN];

int FlujoBloqueante(int u, int t, int f) {
    if (u == t) return f; int fluido = 0;
    for (int i = 0; i < grafo[u].size(); ++i) {
        int v = grafo[u][i];
        if (dist[u] + 1 > dist[v]) continue;
        int fv = FlujoBloqueante(v, t,
            min(f - fluido, cap[u][v] - flujo[u][v]));
        flujo[u][v] += fv; flujo[v][u] -= fv;
        fluido += fv; if (fluido == f) break;
    }
    return fluido;
}

int Dinic(int s, int t, int n) {
    int flujo_maximo = dist[t] = 0;
    while (dist[t] < INF) {
        fill(dist, dist + n, INF);
        queue<int> q; q.push(s); dist[s] = 0;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int i = 0; i < grafo[u].size(); ++i) {
                int v = grafo[u][i];
                if (flujo[u][v] == cap[u][v] ||
                    dist[v] <= dist[u] + 1) continue;
                dist[v] = dist[u] + 1, q.push(v);
            }
        }
        if (dist[t] < INF) flujo_maximo +=
            FlujoBloqueante(s, t, INF);
    }
    return flujo_maximo;
}

int main() {
    return 0;
}
