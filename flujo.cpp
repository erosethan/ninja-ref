#include <bits/stdc++.h>
using namespace std;

// Definiciones iniciales.

typedef long long Dato; // Ajustable.
typedef vector<Dato> Vec;
typedef vector<Vec> Mat;
typedef vector<int> Lista;
typedef pair<int, int> Par;
typedef pair<Dato, Dato> Datos;

const Dato INF = 1 << 30;

// EMPAREJAMIENTO BIPARTITO
// Nodos indexados de 0 a n - 1.

struct Bipartito {

    int n;
    Lista pareja;
    vector<bool> lado;
    vector<bool> visitado;
    vector<Lista> aristas;

    Bipartito(int N) : lado(N), pareja(N),
        visitado(N), aristas(N), n(N) {}

    void AgregarArista(int u, int v) {
        aristas[u].push_back(v);
        aristas[v].push_back(u);
    }

    void AgregarIzq(int u) { lado[u] = true; }
    void AgregarDer(int u) { lado[u] = false; }

    int CaminoIncremental(int u) {
        visitado[u] = true;
        for (int i = 0; i < aristas[u].size(); ++i)
            if (pareja[aristas[u][i]] == -1) 
                return pareja[aristas[u][i]] = u;
        for (int i = 0; i < aristas[u].size(); ++i) {
            int v = aristas[u][i];
            if (visitado[pareja[v]]) continue;
            if (CaminoIncremental(pareja[v]) != -1)
                return pareja[v] = u;
        }
        return -1;
    }

    // Maximo emparejamiento en grafo bipartito.

    vector<Par> MaxEmparejamiento() {
        fill(pareja.begin(), pareja.end(), -1);
        for (int i = 0; i < n; ++i) {
            if (!lado[i]) continue; CaminoIncremental(i);
            fill(visitado.begin(), visitado.end(), false);
        }
        vector<Par> pares;
        for (int i = 0; i < n; ++i)
            if (!lado[i] && pareja[i] != -1)
                pares.push_back(Par(pareja[i], i));
        return pares; // Cardinalidad = pares.size()
    }
};


// EMPAREJAMIENTO BIPARTITO DE COSTO MAX/MIN
// Nodos indexados de 0 a n - 1, diferencia
// entre nodos en el conjunto izquierdo y derecho.
// Es posible que alguna variable se desborde y se
// cicle, para evitarlo cambien Dato a long long.

struct BipartitoCosto {

    int n, s;
    Mat costo;
    Vec slack, etiqueta;
    Lista pareja, retorno;
    vector<bool> visitado;

    // Emparejamiento de costo maximo S =  1
    // Emparejamiento de costo minimo S = -1
    
    BipartitoCosto(int N, int S = 1)
        : costo(N, Vec(N, -INF * S)), s(S),
          slack(2*N), etiqueta(2*N), pareja(2*N),
          retorno(2*N), visitado(2*N), n(N) {}

    void AgregarArista(
        int u, int v, int c) {
        costo[u][v] = c * s;
    }

    vector<Par> EmparejamientoOptimo() {

        fill(pareja.begin(), pareja.end(), -1);
        fill(etiqueta.begin(), etiqueta.end(), 0);
        for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j)
            etiqueta[i] = max(etiqueta[i], costo[i][j]);
                
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j)
                slack[j + n] = etiqueta[i] +
                    etiqueta[j + n] - costo[i][j];
            fill(visitado.begin(), visitado.end(), false);
            fill(retorno.begin(), retorno.end(), i);
            visitado[i] = true;
            
            bool emparejado = false;
            for (int j = 0; !emparejado; ++j) {
                int t = n; for (; t < 2 * n; ++t) {
                    if (visitado[t]) continue;
                    if (!slack[t]) break;
                }
                if (t < 2 * n) {
                    visitado[t] = true;
                    if (pareja[t] == -1) {
                        emparejado = true;
                        for (int p; ; t = p) {
                            pareja[t] = retorno[t];
                            p = pareja[retorno[t]];
                            pareja[retorno[t]] = t;
                            if (retorno[t] == i) break;
                        }
                    } else {
                        visitado[t = pareja[t]] = true;
                        for (int k = 0; k < n; ++k) {
                            Dato new_slack = etiqueta[t] +
                                etiqueta[k + n] - costo[t][k];
                            if (new_slack < slack[k + n]) {
                                slack[k + n] = new_slack;
                                retorno[k + n] = t;
                            }
                        }
                    }
                } else {
                    Dato d = INF;
                    for (int k = n; k < 2 * n; ++k)
                        if (slack[k]) d = min(d, slack[k]);
                    for (int k = 0; k < n; ++k)
                        if (visitado[k]) etiqueta[k] -= d;
                    for (int k = n; k < 2 * n; ++k)
                        if (!visitado[k]) slack[k] -= d;
                        else etiqueta[k] += d;
                }
            }
        }
        vector<Par> pares;
        for (int i = 0; i < n; ++i)
            if (costo[i][pareja[i] - n] != s * -INF)
                pares.push_back(Par(i, pareja[i] - n));
        return pares; // Emparejamiento optimo.
    }
};

// FLUJO MAXIMO
// Nodos indexados de 0 a n - 1.

struct Flujo {

    int n;
    Lista padre;
    Mat cap, flujo;
    vector<Dato> dist;
    vector<Lista> aristas;

    Flujo(int N) : dist(N), padre(N), aristas(N),
        cap(N, Vec(N)), flujo(N, Vec(N)), n(N) {}

    void AgregarArista(int u, int v, Dato c) {
        flujo[v][u] += c; // Solo dirigidas!
        cap[u][v] += c, cap[v][u] += c;
        aristas[u].push_back(v);
        aristas[v].push_back(u);
    }

    // Flujo maximo mediante Edmonds-Karp O(VE^2).

    Dato ActualizarFlujo(int u, Dato f) {
        int p = padre[u];
        if (p == u) return f;
        f = ActualizarFlujo(p, min(
            f, cap[p][u] - flujo[p][u]));
        flujo[p][u] += f;
        flujo[u][p] -= f;
        return f;
    }

    Dato AumentarFlujo(int s, int t) {
        fill(padre.begin(), padre.end(), -1);
        queue<int> q; q.push(s); padre[s] = s;
        while (!q.empty()) {
            int u = q.front();
            q.pop(); if (u == t) break;
            for (int i = 0; i < aristas[u].size(); ++i) {
                int v = aristas[u][i];
                if (flujo[u][v] == cap[u][v] ||
                    padre[v] != -1) continue;
                padre[v] = u, q.push(v);
            }
        }
        if (padre[t] == -1) return 0;
        return ActualizarFlujo(t, INF);
    }

    Dato EdmondsKarp(int s, int t) {
        Dato flujo_maximo = 0, f;
        while (f = AumentarFlujo(s, t))
            flujo_maximo += f;
        return flujo_maximo;
    }

    // Flujo maximo mediante Dinic O(V^2E).

    Dato FlujoBloqueante(int u, int t, Dato f) {
        if (u == t) return f; Dato fluido = 0;
        for (int i = 0; i < aristas[u].size(); ++i) {
            if (fluido == f) break; int v = aristas[u][i];
            if (dist[u] + 1 > dist[v]) continue;
            Dato fv = FlujoBloqueante(v, t,
                min(f - fluido, cap[u][v] - flujo[u][v]));
            flujo[u][v] += fv, fluido += fv;
            flujo[v][u] -= fv;
        }
        return fluido;
    }

    Dato Dinic(int s, int t) {
        Dato flujo_maximo = dist[t] = 0;
        while (dist[t] < INF) {
            fill(dist.begin(), dist.end(), INF);
            queue<int> q; q.push(s); dist[s] = 0;
            while (!q.empty()) {
                int u = q.front(); q.pop();
                for (int i = 0; i < aristas[u].size(); ++i) {
                    int v = aristas[u][i];
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
};

// FLUJO MEMORIA OPTIMIZADA Y
// FLUJO MAXIMO DE COSTO MINIMO
// Nodos indexados de 0 a n - 1.
// No utiliza matrices de adyacencia.

struct FlujoCosto {

    struct AristaFlujo {
        int dst; Dato cap,
        flujo, peso, npeso;
        AristaFlujo* residual;

        AristaFlujo(int d, Dato f, Dato c)
            : dst(d), flujo(f), cap(c) {}

        Dato AumentarFlujo(Dato f) {
            residual->flujo -= f;
            this->flujo += f;
            return peso * f;
        }
    };

    vector< vector<AristaFlujo*> > aristas;
    int n; vector<Par> prv; vector<Dato> dist;

    FlujoCosto(int N) : aristas(N),
        prv(N), dist(N), n(N) {}

    // No omitan el siguiente destructor!
    
    ~FlujoCosto() { for (int i = 0; i < n; ++i)
        for (int j = 0; j < aristas[i].size(); ++j)
            delete aristas[i][j];
    }

    // Para aristas bidireccionales agreguen dos aristas
    // dirigidas. Si las aristas no son ponderadas dejen
    // el ultimo parametro con el valor por defecto.

    void AgregarArista(int u, int v, Dato c, Dato p = 0) {
        AristaFlujo* uv = new AristaFlujo(v, 0, c);
        AristaFlujo* vu = new AristaFlujo(u, c, c);
        uv->residual = vu, vu->residual = uv;
        uv->peso = uv->npeso =  p;
        vu->peso = vu->npeso = -p;
        aristas[u].push_back(uv);
        aristas[v].push_back(vu);
    }

    // Dinic para flujo maximo con memoria optimizada.
    // Prefieran esta version solo cuando n > 5,000.

    Dato FlujoBloqueante(int u, int t, Dato f) {
        if (u == t) return f; Dato fluido = 0;
        for (int i = 0; i < aristas[u].size(); ++i) {
            if (fluido == f) break;
            AristaFlujo* v = aristas[u][i];
            if (dist[u] + 1 == dist[v->dst]) {
                Dato fv = FlujoBloqueante(v->dst, t,
                    min(f - fluido, v->cap - v->flujo));
                v->AumentarFlujo(fv), fluido += fv;
            }
        }
        return fluido;
    }

    Dato Dinic(int s, int t) {
        Dato flujo_maximo = dist[t] = 0;
        while (dist[t] < INF) {
            fill(dist.begin(), dist.end(), INF);
            queue<int> q; q.push(s); dist[s] = 0;
            while (!q.empty()) {
                int u = q.front(); q.pop();
                for (int i = 0; i < aristas[u].size(); ++i) {
                    AristaFlujo* v = aristas[u][i];
                    if (dist[v->dst] < INF) continue;
                    if (v->flujo == v->cap) continue;
                    dist[v->dst] = dist[u] + 1;
                    q.push(v->dst);
                }
            }
            if (dist[t] < INF) flujo_maximo +=
                FlujoBloqueante(s, t, INF);
        }
        return flujo_maximo;
    }

    // Flujo de costo minimo en O(VElogV * flow). Si dejan el
    // valor por defecto del parametro k saca el flujo maximo.

    void RecalcularCosto(const vector<Dato>& pi) {
        for (int u = 0; u < n; ++u) {
            for (int i = 0; i < aristas[u].size(); ++i) {
                AristaFlujo* v = aristas[u][i];
                v->npeso = v->npeso + pi[u] - pi[v->dst];
            }
        }
    }

    Datos ActualizarFlujo(int u, Dato f) {
        int p = prv[u].first;
        int i = prv[u].second;
        if (p == -1) return Datos(f, 0);
        AristaFlujo* pu = aristas[p][i];

        Datos res = ActualizarFlujo(
            p, min(f, pu->cap - pu->flujo));
        res.second += pu->AumentarFlujo(
            res.first); return res;
    }

    Datos AumentarFlujo(int s, int t, Dato f) {
        vector<Dato> dist(n, INF);
        fill(prv.begin(), prv.end(), Datos(-1, -1));
        priority_queue<Datos, vector<Datos>,
                       greater<Datos> > pq;
        pq.push(Datos(0, s)); dist[s] = 0;
        
        while (!pq.empty()) {
            int u = pq.top().second;
            int p = pq.top().first; pq.pop();
            if (dist[u] < p) continue;        
            for (int i = 0; i < aristas[u].size(); ++i) {
                AristaFlujo* v = aristas[u][i];
                if (v->flujo == v->cap) continue;
                if (dist[u] + v->npeso < dist[v->dst]) {
                    dist[v->dst] = dist[u] + v->npeso;
                    pq.push(Datos(dist[v->dst], v->dst));
                    prv[v->dst].second = i;
                    prv[v->dst].first = u;
                }
            }
        }
        if (dist[t] == INF)
            return Datos(0, 0);
        RecalcularCosto(dist);
        return ActualizarFlujo(t, f);
    }

    Datos FlujoCostoMin(int s, int t, int k = -1) {
        vector<Dato> dist(n, INF); dist[s] = 0;
        for (int i = 0; i < n; ++i) {
            for (int u = 0; u < n; ++u) {
                if (dist[u] == INF) continue;
                for (int j = 0; j < aristas[u].size(); ++j) {
                    AristaFlujo* v = aristas[u][j];
                    if (v->flujo < v->cap)  dist[v->dst] = min(
                        dist[v->dst], dist[u] + v->npeso);
                }
            }
        }
        RecalcularCosto(dist);

        Datos flujo_costo(0, 0);
        while (flujo_costo.first < k || k == -1) {
            Datos fc = AumentarFlujo(s, t, (k == -1)?
                INF: k - flujo_costo.first);
            flujo_costo.second += fc.second;
            flujo_costo.first += fc.first;
            if (!fc.first) break;
        }
        return flujo_costo;
    }
};

int main() {
    return 0;
}
