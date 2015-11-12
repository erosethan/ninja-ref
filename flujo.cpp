#include <bits/stdc++.h>
using namespace std;

// Definiciones iniciales.

typedef int Flujo; // Ajustable.

typedef vector<int> Lista;
typedef pair<int, int> Par;
typedef vector<Flujo> Flujo1D;
typedef vector<Flujo1D> Flujo2D;

const Flujo FINF = 1 << 30;

// EMPAREJAMIENTO BIPARTITO
// Nodos indexados de 0 a n - 1.

struct Bipartito {

    int n; Lista pareja;
    vector<Lista> aristas;
    vector<bool> lado, visitado;

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

// FLUJO MAXIMO
// Nodos indexados de 0 a n - 1.

struct GrafoFlujo {

    int n; vector<Lista> aristas;
    Flujo2D cap, flujo; Lista padre, dist;
    
    GrafoFlujo(int N) : dist(N), padre(N), aristas(N),
        cap(N, Flujo1D(N)), flujo(N, Flujo1D(N)), n(N) {}

    void AgregarArista(int u, int v, Flujo c) {
        flujo[v][u] += c; // Solo dirigidas!
        cap[u][v] += c, cap[v][u] += c;
        aristas[u].push_back(v);
        aristas[v].push_back(u);
    }

    // Flujo maximo mediante Edmonds-Karp O(VE^2).

    Flujo ActualizarFlujo(int u, Flujo f) {
        int p = padre[u];
        if (p == u) return f;
        f = ActualizarFlujo(p, min(f,
            cap[p][u] - flujo[p][u]));
        flujo[p][u] += f;
        flujo[u][p] -= f;
        return f;
    }

    Flujo AumentarFlujo(int s, int t) {
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
        return ActualizarFlujo(t, FINF);
    }

    Flujo EdmondsKarp(int s, int t) {
        Flujo flujo_maximo = 0, f;
        while (f = AumentarFlujo(s, t))
            flujo_maximo += f;
        return flujo_maximo;
    }

    // Flujo maximo mediante Dinic O(V^2E).

    Flujo FlujoBloqueante(int u, int t, Flujo f) {
        if (u == t) return f; Flujo fluido = 0;
        for (int i = 0; i < aristas[u].size(); ++i) {
            if (fluido == f) break; int v = aristas[u][i];
            if (dist[u] + 1 > dist[v]) continue;
            Flujo fv = FlujoBloqueante(v, t,
                min(f - fluido, cap[u][v] - flujo[u][v]));
            flujo[u][v] += fv, fluido += fv;
            flujo[v][u] -= fv;
        }
        return fluido;
    }

    Flujo Dinic(int s, int t) {
        Flujo flujo_maximo = dist[t] = 0;
        while (dist[t] < INT_MAX) {
            fill(dist.begin(), dist.end(), INT_MAX);
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
            if (dist[t] < INT_MAX) flujo_maximo +=
                FlujoBloqueante(s, t, FINF);
        }
        return flujo_maximo;
    }
};

// Definiciones adicionales.

typedef int Costo; // Ajustable.

typedef vector<Costo> Costo1D;
typedef vector<Costo1D> Costo2D;
typedef pair<Costo, int> CostoNodo;
typedef pair<Flujo, Costo> FlujoCosto;

const double ERROR = 1e-9;
const Costo CINF = 1 << 30;

// Tolerancia en flotantes.

bool Igual(double a, double b) {
    return fabs(a - b) < ERROR;
}

// EMPAREJAMIENTO BIPARTITO DE COSTO MAX/MIN
// Nodos indexados de 0 a n - 1, diferencia
// entre nodos en el conjunto izquierdo y derecho.
// Es posible que alguna variable se desborde y se
// cicle, para evitarlo cambien Dato a long long.

struct BipartitoCosto {

    Lista pareja, retorno; vector<bool> visitado;
    int n, s; Costo1D slack, etiqueta; Costo2D costo;

    // Emparejamiento de costo maximo S =  1
    // Emparejamiento de costo minimo S = -1
    
    BipartitoCosto(int N, int S = 1)
        : costo(N, Costo1D(N, S * -CINF)), s(S),
          slack(2 * N), etiqueta(2 * N), pareja(2 * N),
          retorno(2 * N), visitado(2 * N), n(N) {}

    void AgregarArista(int u, int v,
        Costo c) { costo[u][v] = c * s; }

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
                    if (Igual(slack[t], 0)) break;
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
                            Costo new_slack = etiqueta[t] +
                                etiqueta[k + n] - costo[t][k];
                            if (!Igual(new_slack, slack[k + n])
                                && new_slack < slack[k + n]) {
                                slack[k + n] = new_slack;
                                retorno[k + n] = t;
                            }
                        }
                    }
                } else {
                    Costo d = CINF;
                    for (int k = n; k < 2 * n; ++k)
                        if (!Igual(slack[k], 0))
                            d = min(d, slack[k]);
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
            if (!Igual(costo[i][pareja[i] - n], s * -CINF))
                pares.push_back(Par(i, pareja[i] - n));
        return pares; // Emparejamiento optimo.
    }
};

// FLUJO MEMORIA OPTIMIZADA Y
// FLUJO MAXIMO DE COSTO MINIMO
// Nodos indexados de 0 a n - 1.
// No utiliza matrices de adyacencia.

struct GrafoFlujoCosto {

    struct AristaFlujo {

        int dst; AristaFlujo* residual;
        Flujo cap, flujo; Costo peso, npeso;

        AristaFlujo(int d, Flujo f, Flujo c)
            : dst(d), flujo(f), cap(c) {}

        Costo AumentarFlujo(Flujo f) {
            residual->flujo -= f;
            this->flujo += f;
            return peso * f;
        }
    };

    int n; vector<Par> prv; Lista dist;
    vector< vector<AristaFlujo*> > aristas;

    GrafoFlujoCosto(int N) : n(N),
        aristas(N), prv(N), dist(N) {}
    
    ~GrafoFlujoCosto() { for (int i = 0; i < n; ++i)
        for (int j = 0; j < aristas[i].size(); ++j)
            delete aristas[i][j]; // NO OMITIR!!!
    }

    // Para aristas bidireccionales agreguen dos aristas
    // dirigidas. Si las aristas no son ponderadas dejen
    // el ultimo parametro con el valor por defecto.

    void AgregarArista(int u, int v, Flujo c, Costo p = 0) {
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

    Flujo FlujoBloqueante(int u, int t, Flujo f) {
        if (u == t) return f; Flujo fluido = 0;
        for (int i = 0; i < aristas[u].size(); ++i) {
            if (fluido == f) break;
            AristaFlujo* v = aristas[u][i];
            if (dist[u] + 1 == dist[v->dst]) {
                Flujo fv = FlujoBloqueante(v->dst, t,
                    min(f - fluido, v->cap - v->flujo));
                v->AumentarFlujo(fv), fluido += fv;
            }
        }
        return fluido;
    }

    Flujo Dinic(int s, int t) {
        Flujo flujo_maximo = dist[t] = 0;
        while (dist[t] < INT_MAX) {
            fill(dist.begin(), dist.end(), INT_MAX);
            queue<int> q; q.push(s); dist[s] = 0;
            while (!q.empty()) {
                int u = q.front(); q.pop();
                for (int i = 0; i < aristas[u].size(); ++i) {
                    AristaFlujo* v = aristas[u][i];
                    if (dist[v->dst] < INT_MAX) continue;
                    if (v->flujo == v->cap) continue;
                    dist[v->dst] = dist[u] + 1;
                    q.push(v->dst);
                }
            }
            if (dist[t] < INT_MAX) flujo_maximo +=
                FlujoBloqueante(s, t, FINF);
        }
        return flujo_maximo;
    }

    // Flujo de costo minimo en O(VElogV * flow). Si dejan el
    // valor por defecto del parametro k saca el flujo maximo.

    void RecalcularCosto(const Costo1D& pi) {
        for (int u = 0; u < n; ++u) {
            for (int i = 0; i < aristas[u].size(); ++i) {
                AristaFlujo* v = aristas[u][i];
                v->npeso = v->npeso + pi[u] - pi[v->dst];
            }
        }
    }

    FlujoCosto ActualizarFlujo(int u, Flujo f) {
        int p = prv[u].first, i = prv[u].second;
        if (p == -1) return FlujoCosto(f, 0);
        AristaFlujo* pu = aristas[p][i];

        FlujoCosto res = ActualizarFlujo(
            p, min(f, pu->cap - pu->flujo));
        res.second += pu->AumentarFlujo(
            res.first); return res;
    }

    FlujoCosto AumentarFlujo(int s, int t, Flujo f) {
        Costo1D dist(n, CINF);
        fill(prv.begin(), prv.end(), Par(-1, -1));
        priority_queue<CostoNodo, vector<CostoNodo>,
                       greater<CostoNodo> > pq;
        pq.push(FlujoCosto(0, s)); dist[s] = 0;
        
        while (!pq.empty()) {
            int u = pq.top().second;
            Costo p = pq.top().first; pq.pop();
            if (!Igual(dist[u], p)) continue;
            for (int i = 0; i < aristas[u].size(); ++i) {
                AristaFlujo* v = aristas[u][i];
                if (v->flujo == v->cap) continue;
                Costo ndist = dist[u] + v->npeso;
                if (!Igual(ndist, dist[v->dst]) &&
                           ndist < dist[v->dst]) {
                    dist[v->dst] = dist[u] + v->npeso;
                    pq.push(CostoNodo(ndist, v->dst));
                    prv[v->dst].second = i;
                    prv[v->dst].first = u;
                }
            }
        }
        if (Igual(dist[t], CINF))
            return FlujoCosto(0, 0);
        RecalcularCosto(dist);
        return ActualizarFlujo(t, f);
    }

    FlujoCosto FlujoCostoMin(int s, int t, Flujo k = FINF) {
        Costo1D dist(n, CINF); dist[s] = 0;
        for (int i = 0; i < n; ++i) {
            for (int u = 0; u < n; ++u) {
                if (Igual(dist[u], CINF)) continue;
                for (int j = 0; j < aristas[u].size(); ++j) {
                    AristaFlujo* v = aristas[u][j];
                    if (v->flujo < v->cap) dist[v->dst] =
                        min(dist[v->dst], dist[u] + v->npeso);
                }
            }
        }
        RecalcularCosto(dist);

        FlujoCosto flujo_costo(0, 0);
        while (flujo_costo.first < k) {
            FlujoCosto fc = AumentarFlujo(
                s, t, k - flujo_costo.first);
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
