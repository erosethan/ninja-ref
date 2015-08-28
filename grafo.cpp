#include <bits/stdc++.h>
using namespace std;

// Definiciones iniciales
typedef pair<int, int> Arista;
typedef pair<int, Arista> PesoArista;

const int INF = 1 << 30;
const int MAXN = 1000;

int in_degree[MAXN];
vector<int> grafo[MAXN];
vector<Arista> grafo_peso[MAXN];

void AgregarArista(int u, int v) {
    grafo[u].push_back(v);
    ++in_degree[v];
}

// Detecta ciclos en una componente conexa bidireccional.
// Devuelve -1 si no existen ciclos. En caso de existir,
// los nodos del ciclo se guardan en el arreglo ciclo.

int top_ciclo;
int ciclo[MAXN];
bool ciclo_activo;
bool visitado[MAXN];

int EncontrarCiclo_(int u, int p) {
    visitado[u] = true;
    for (int i = 0; i < grafo[u].size(); ++i) {
        if (grafo[u][i] == p) continue;
        int v = grafo[u][i], cic;
        if (visitado[v]) return ciclo[top_ciclo++] = v; 
        if ((cic = EncontrarCiclo_(v, u)) >= 0) {
            if (ciclo_activo) ciclo[top_ciclo++] = v;
            if (cic == u) ciclo_activo = false;
            return cic;
        }
    }
    return -1;
}

int EncontrarCiclo(int u) {
    ciclo_activo = true, top_ciclo = 0;
    return EncontrarCiclo_(u, -1);
}

// Orden topologico de un grafo dirigido aciclico.
// Los indices de los nodos se asumen de 0 a n - 1.

vector<int> orden_topo; // Resultado

void OrdenTopologico_(int u) {
    visitado[u] = true;
    for (int i = 0; i < grafo[u].size(); ++i)
        if (!visitado[grafo[u][i]])
            OrdenTopologico_(grafo[u][i]);
    orden_topo.push_back(u);
}

void OrdenTopologico(int n) {
    orden_topo.clear();
    fill(visitado, visitado + n, false);
    for (int i = 0; i < n; ++i)
        if (!visitado[i]) OrdenTopologico_(i);
}

// Obtener componentes fuertemente conexas en un grafo
// dirigido. Asume los nodos indexados de 0 a n - 1.

int numeracion, top_activo;
int low[MAXN], num[MAXN], activo[MAXN];
vector< vector<int> > CFCs; // Resultado

void ObtenerCFCs_(int u) {
    activo[top_activo++] = u;
    low[u] = num[u] = ++numeracion;
    for (int i = 0; i < grafo[u].size(); ++i) {
        int v = grafo[u][i];
        if (!num[v]) ObtenerCFCs_(v);
        low[u] = min(low[u], low[v]);
    }
    if (low[u] == num[u]) {
        vector<int> CFC;
        while (activo[top_activo - 1] != u) {
            CFC.push_back(activo[--top_activo]);
            low[activo[top_activo]] = INF;
        }
        CFC.push_back(activo[--top_activo]);
        low[u] = INF; CFCs.push_back(CFC);
    }
}

void ObtenerCFCs(int n) {
    CFCs.clear();
    numeracion = 0;
    fill(num, num + n, 0);
    fill(low, low + n, 0);
    for (int i = 0; i < n; ++i)
        if (!num[i]) ObtenerCFCs_(i);
}

// Detecta los puentes y puntos de articulacion en
// un grafo bidireccional. Indices de 0 a n - 1.

vector<int> puentes[MAXN]; // Resultado
bool punto_art[MAXN]; // Resultado

void PuntosArtPuentes_(int u, int p) {
    int hijos = 0;
    low[u] = num[u] = ++numeracion;
    for (int i = 0; i < grafo[u].size(); ++i) {
        int v = grafo[u][i];
        if(v == p) continue;
        if (!num[v]) {
            ++hijos;
            PuntosArtPuentes_(v, u);
            if (low[v] > num[u]) {
                puentes[u].push_back(v);
                puentes[v].push_back(u);
            }
            low[u] = min(low[u], low[v]);
            punto_art[u] |= low[v] >= num[u];
        } else low[u] = min(low[u], num[v]);
    }
    if (p == -1) punto_art[u] = hijos > 1;
}

void PuntosArtPuentes(int n) {
    numeracion = 0;
    fill(num, num + n, 0);
    fill(low, low + n, 0);
    fill(punto_art, punto_art + n, false);
    for (int i = 0; i < n; ++i)
        puentes[i].clear();
    for (int i = 0; i < n; ++i)
        if (!num[i]) PuntosArtPuentes_(i, -1);
}

// Estructura de conjuntos disjuntos.
// Conjuntos indexados de 0 a n - 1.

struct UnionFind {
    int nconjuntos;
    vector<int> padre;
    vector<int> tamano;

    UnionFind(int n) : nconjuntos(n),
        padre(n), tamano(n, 1) {
        for(int i = 0; i < n; ++i)
            padre[i] = i;
    }

    int Encontrar(int u) {
        if (padre[u] == u) return u;
        return padre[u] = Encontrar(padre[u]);
    }

    void Unir(int u, int v) {
        int Ru = Encontrar(u);
        int Rv = Encontrar(v);
        if (Ru == Rv) return;
        -- nconjuntos, padre[Ru] = Rv;
        tamano[Rv] += tamano[Ru];
    }

    bool MismoConjunto(int u, int v) {
        return Encontrar(u) == Encontrar(v);
    }

    int TamanoConjunto(int u) {
        return tamano[Encontrar(u)];
    }
};

// Saber si un grafo es bicoloreable o bipartito.
// Se asumen los nodos indexados de 0 a n - 1.

char color[MAXN];

bool Bicolorear_(int u, int c) {
    color[u] = c;
    for (int i = 0; i < grafo[u].size(); ++i) {
        int v = grafo[u][i];
        if (color[v] == 1 - c) continue;
        if (color[v] == c) return false;
        if (!Bicolorear_(v, 1 - c)) return false;
    }
}

bool Bicolorear(int n) {
    fill(color, color + n, -1);
    for (int i = 0; i < n; ++i)
        if (color[i] == -1 &&
            !Bicolorear_(i, 0))
            return false;
    return true;
}

// Busqueda en amplitud. Nodos con indice del 0 al n - 1.

vector<int> BFS(int o, int n) {
    vector<int> dist(n, INF);
    queue<int> q; q.push(o); dist[o] = 0;
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int i = 0; i < grafo[u].size(); ++i) {
            int v = grafo[u][i];
            if (dist[u] + 1 < dist[v]) {
                dist[v] = dist[u] + 1;
                q.push(v);
            }
        }
    }
    return dist;
}

// Algoritmo de Dijkstra. Nodos indexados del 0 al n - 1.

vector<int> Dijkstra(int o, int n) {
    vector<int> dist(n, INF);
    priority_queue<Arista, vector<Arista>,
                   greater<Arista> > pq;
    pq.push(Arista(0, o)); dist[o] = 0;
    
    while (!pq.empty()) {
        int u = pq.top().second;
        int p = pq.top().first; pq.pop();
        if (dist[u] < p) continue;        
        for (int i = 0; i < grafo_peso[u].size(); ++i) {
            p = grafo_peso[u][i].second;
            int v = grafo_peso[u][i].first;
            if (dist[u] + p < dist[v]) {
                dist[v] = dist[u] + p;
                pq.push(Arista(dist[v], v));
            }
        }
    }
    return dist;
}

// Dijkstra version lineal. Nodos indexados del 0 al n - 1.
// ¡Peligro! Recuerden cuidar el peso maximo de las aristas.

const int MAXP = 100; // Peso maximo

vector<int> DijkstraLineal(int o, int n) {
    vector<int> dist(n, INF);
    vector<bool> proc(n, false);
    vector< queue<int> > q(MAXP);
    q[0].push(o); dist[o] = 0;
    int qi = 0, total = 1;
    
    while (total) {
        if (!q[qi].empty()) {
            int u = q[qi].front(); q[qi].pop(), --total;
            if (proc[u]) continue; proc[u] = true;
            for (int i = 0; i < grafo_peso[u].size(); ++i) {
                int v = grafo_peso[u][i].first;
                int p = grafo_peso[u][i].second;
                if (dist[u] + p < dist[v])
                    q[(qi + p) % MAXP].push(v), ++total;
            }
        } else qi = (qi + 1) % MAXP;
    }
    return dist;
}

// Algoritmo de Floyd-Warshall. Nodos con indice del 0 al n - 1.

int dist[MAXN][MAXN]; // Matriz de adyacencia

void FloydWarshall(int n) {
    for (int k = 0; k < n; ++k)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                dist[i][j] = min(dist[i][j],
                    dist[i][k] + dist[k][j]);
}

// Arbol de expansion minima por Kruskal. Nodos del 0 al n - 1. 

vector<PesoArista> Kruskal(int n) {
    vector<PesoArista> aristas;
    for (int u = 0; u < n; ++u) {
        for (int i = 0; i < grafo_peso[u].size(); ++i) {
            int v = grafo_peso[u][i].first;
            int p = grafo_peso[u][i].second;
            aristas.push_back(PesoArista(
                p, Arista(u, v)));
        }
    }
    sort(aristas.begin(), aristas.end());
    
    UnionFind uf(n);
    vector<PesoArista> mst;
    for (int i = 0; i < aristas.size(); ++i) {
        int u = aristas[i].second.first;
        int v = aristas[i].second.second;
        if (uf.MismoConjunto(u, v)) continue;
        uf.Unir(u, v); mst.push_back(aristas[i]);
    }
    return mst;
}

// Algoritmo de Bellman Ford. Nodos indexados de 0 a n - 1.

vector<int> Bellmanford(int o, int n) {
    vector<int> dist(n, INF); dist[o] = 0;
    for (int i = 0; i < n; ++i) {
        for (int u = 0; u < n; ++u) {
            if (dist[u] == INF) continue;
            for (int j = 0; j < grafo_peso[u].size(); ++j) {
                int v = grafo_peso[u][j].first;
                int p = grafo_peso[u][j].second;
                dist[v] = min(dist[v], dist[u] + p);
            }
        }
    }
    bool ciclo_neg = false;
    for (int u = 0; u < n; ++u) {
        if (dist[u] == INF) continue;
        for (int j = 0; j < grafo_peso[u].size(); ++j) {
            int v = grafo_peso[u][j].first;
            int p = grafo_peso[u][j].second;
            ciclo_neg |= dist[u] + p < dist[v];
        }
    }
    if (!ciclo_neg) return dist;
    for (int u = 0; u < n; ++u)
        if (dist[u] < INF) dist[u] = -INF;
    return dist;
}

int main() {
    return 0;
}
