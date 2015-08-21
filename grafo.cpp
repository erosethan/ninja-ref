#include <bits/stdc++.h>
using namespace std;

// Definiciones iniciales
typedef int Arista;
typedef vector<Arista> Nodo;
typedef pair<int, int> AristaPeso;
typedef vector<AristaPeso> NodoPeso;

const int INF = 1 << 30;
const int MAXN = 100000;
const int MAXM = 100000;

Nodo grafo[MAXN];
NodoPeso grafo_peso[MAXN];
int in_degree[MAXN];

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
            low[activo[top_activo]] = MAXN;
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

Nodo puentes[MAXN]; // Resultado
bool punto_art[MAXN]; // Resultado

void PuntosArtPuentes_(int u, int p) {
    low[u] = num[u] = ++numeracion;
    for (int i = 0; i < grafo[u].size(); ++i) {
        if (grafo[u][i] == p) continue;
        int v = grafo[u][i];
        if (!num[v]) PuntosArtPuentes_(v, u);
        if (low[v] > num[u]) {
            puentes[u].push_back(v);
            puentes[v].push_back(u);
        }
        punto_art[u] |= low[v] >= num[u];
        low[u] = min(low[u], low[v]);
    }
}

void PuntosArtPuentes(int n) {
    numeracion = 0;
    fill(num, num + n, 0);
    fill(low, low + n, 0);
    fill(punto_art, punto_art + n, false);
    for(int i = 0; i < n; ++i)
        puentes[i].clear();
    for(int i = 0; i < n; ++i) {
        if(num[i]) continue;
        PuntosArtPuentes_(i, -1);
        punto_art[i] = grafo[i].size() > 1;
    }
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
    priority_queue<AristaPeso,
                   vector<AristaPeso>,
                   greater<AristaPeso> > pq;
    pq.push(AristaPeso(0, o)); dist[o] = 0;
    
    while (!pq.empty()) {
        int u = pq.top().second;
        int p = pq.top().first; pq.pop();
        if (dist[u] < p) continue;        
        for (int i = 0; i < grafo_peso[u].size(); ++i) {
            p = grafo_peso[u][i].second;
            int v = grafo_peso[u][i].first;
            if (dist[u] + p < dist[v]) {
                dist[v] = dist[u] + p;
                pq.push(AristaPeso(dist[v], v));
            }
        }
    }
    return dist;
}

const int MAXP = 100;
vector<int> DijkstraLineal(int o, int n) {
    vector < queue<int> > q(MAXP);
    vector <bool> procesado(n, false);
    vector<int > dist(n, INF);
    q[0].push(o); dist[0] = 0;
    int qi = 0, total = 1;
    while(total) {
        if(!q[qi].empty()) {
            int u = q[qi].front();
            q[qi].pop();
            total--;
            if(procesado[u])continue;
            procesado[u] = true;
            for(int i = 0; i < grafo_peso[u].size(); ++i) {
                int v = grafo_peso[u][i].first;
                int p = grafo_peso[u][i].second;
                if(dist[u] + p < dist[v]) {
                    q[(qi + p) % MAXP].push(v);
                    ++total;
                }
            }
        }
        else qi = (qi + 1) % MAXP;
    }
}
// O(n³)
//dist matriz de adyacencia
/*
void floy(int n){
    for(int k = 0; k < n; k++)
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
}
*/

//por que gama y alichos no sabian exp binaria
long long Exp(long long a, long long n, long long m) {
    long long res = 1, p = a;
    for(long long i = 0; 1ll << i <= n; i++) {
        if(n & 1ll << i) res = (res * p) % m;
        p = (p * p) % m;
    }
    return res;
}

//multiplicacion binaria
long long Mult(long long a, long long n, long long m) {
    long long res = 0, p = a;
    for(long long i = 0; 1ll << i <= n; i++) {
        if(n & 1ll << i) res = (res + p) % m;
        p = (p + p) % m;
    }
    return res;
}

//kruskal
typedef pair <int, AristaPeso> KArista;
vector <KArista> kruskal(int n) {
    vector <KArista> mst;
    vector <KArista> aristas;
    for(int u = 0; u < n; u++)
        for(int i = 0; i < grafo_peso[u].size(); i++) {
            int v = grafo_peso[u][i].first;
            int p = grafo_peso[u][i].second;
            aristas.push_back(KArista(p, AristaPeso(u, v)));
        }
    sort(aristas.begin(), aristas.end());
    UnionFind uf(n);
    for(int i = 0; i < aristas.size(); ++i) {
        int u = aristas[i].second.first;
        int v = aristas[i].second.second;
        if(uf.MismoConjunto(u, v)) continue;
        uf.Unir(u, v);
        mst.push_back(aristas[i]);
    }
    return mst;
}



int main() {
    return 0;
}
vector<int> Bellmanford(int o, int n) {
    vector<int >dist(n, INF);
    dist[o] = 0;
    for(int i = 0; i < n; i++) {
        for(int u = 0; u < n; u++) {
            if(dist[u] == INF) continue;
            for(int j = 0; j < grafo_peso[u].size(); j++) {
                int v = grafo_peso[u][j].first;
                int p = grafo_peso[u][j].second;
                dist[v] = min(dist[v], dist[u] + p);
            }
        }
    }
    bool ciclo_neg = false;
    for(int u = 0; u < n; u++) {
        if(dist[u] == INF) continue;
        for(int j = 0; j < grafo_peso[u].size(); j++) {
            int v = grafo_peso[u][j].first;
            int p = grafo_peso[u][j].second;
            ciclo_neg |= dist[u] + p < dist[v];
        }
    }
    if(!ciclo_neg) return dist;
    for(int u = 0; u < n; ++u)
        if(dist[u] < INF)dist[u]= -INF;
    return dist;
}