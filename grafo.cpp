#include <bits/stdc++.h>
using namespace std;

typedef int Costo;
typedef pair<int, int> Arista;

const Costo INF = 1 << 30;

// Grafos no ponderados.
// Nodos indexados de 0 a n - 1.
// Grafo(n,  true) -> Bidireccional.
// Grafo(n, false) -> Dirigido.

struct Grafo {

    int n; bool bi;
    vector<vector<int>> ady;
    Grafo(int N, bool B = true)
        : n(N), bi(B), ady(N) {}

    void AgregarArista(int u, int v) {
        if (bi) ady[v].push_back(u);
        ady[u].push_back(v);
    }

    // Detecta ciclos en un grafo o multigrafo.
    // Llamar a DetectarCiclos() devuelve un
    // vector de vectores; cada vector interno es
    // una lista que representa un ciclo del grafo.
    // NOTA: Solo detecta un ciclo por componente.

    vector<int> ciclo;
    vector<char> color;

    void DetectarCiclo(int u, int p) {
        int retorno = bi? 0: 1;
        color[u] = ciclo.empty()? 'G': 'N';
        for (int v : ady[u]) {
            if (v == p && !retorno++) continue;
            if (ciclo.empty() && color[v] == 'G') {
                color[v] = 'A', ciclo.push_back(v);
                if (u != v) color[u] = 'R',
                    ciclo.push_back(u);
            }
            if (color[v] != 'B') continue;

            DetectarCiclo(v, u);
            if (color[u] == 'G' && color[v] == 'R')
                color[u] = 'R', ciclo.push_back(u);
        }
        if (color[u] == 'G') color[u] = 'N';
    }

    vector<vector<int>> DetectarCiclos() {
        vector<vector<int>> ciclos;
        color = vector<char>(n, 'B');
        for (int u = 0; u < n; ++u) {
            if (color[u] != 'B') continue;
            ciclo.clear(); DetectarCiclo(u, n);
            reverse(ciclo.begin(), ciclo.end());
            if (!ciclo.empty())
                ciclos.push_back(ciclo);
        }
        return ciclos;
    }

    // Deteccion de puentes y puntos de articulacion
    // en un grafo o multigrafo bidireccional. Los puentes
    // quedan guardados en un vector de aristas. Lo puntos
    // de articulacion son marcados como true o false
    // respectivamente en un vector de booleanos.

    int tiempo;
    vector<int> label, low;
    vector<Arista> puentes;    // <- Resultado
    vector<bool> articulacion; // <- Resultado

    int PuentesArticulacion(int u, int p) {
        label[u] = low[u] = ++tiempo;

        int hijos = 0, retorno = 0;
        for (int v : ady[u]) {
            if (v == p && !retorno++) continue;
            if (!label[v]) { ++hijos;
                PuentesArticulacion(v, u);
                if (label[u] <= low[v])
                    articulacion[u] = true;
                if (label[u] <  low[v])
                    puentes.push_back(Arista(u, v));
                low[u] = min(low[u], low[v]);
            }
            low[u] = min(low[u], label[v]);
        }
        return hijos;
    }

    void PuentesArticulacion() {
        low = vector<int>(n);
        label = vector<int>(n);
        tiempo = 0, puentes.clear();
        articulacion = vector<bool>(n);
        for (int u = 0; u < n; ++u)
            if (!label[u]) articulacion[u] =
                PuentesArticulacion(u, n) > 1;
    }

    // Deteccion de componentes fuertemente conexas
    // en un grafo dirigido. Las componentes quedan
    // guardadas en un vector de vectores, donde
    // cada vector interno contiene los nodos
    // de una componente fuertemente conexa.

    vector<vector<int>> scc; // <- Resultado
    int top; vector<int> pila;

    void FuertementeConexo(int u) {
        label[u] = low[u] = ++tiempo;

        pila[++top] = u;
        for (int v : ady[u]) {
            if (!label[v]) FuertementeConexo(v);
            low[u] = min(low[u], low[v]);
        }
        if (label[u] == low[u]) {
            vector<int> componente;
            while (pila[top] != u) {
                componente.push_back(pila[top]);
                low[pila[top--]] = n + 1;
            }
            componente.push_back(pila[top--]);
            scc.push_back(componente);
            low[u] = n + 1;
        }
    }

    void FuertementeConexo() {
        low = vector<int>(n);
        label = vector<int>(n);
        tiempo = 0, scc.clear();
        top = -1, pila = vector<int>(n);
        for (int u = 0; u < n; ++u)
            if (!label[u]) FuertementeConexo(u);
    }
    
    // Determina si un grafo bidireccional
    // es bipartito (o bien, es bicoloreable).

    bool EsBipartito() {
        vector<char> color(n, -1);
        for (int u = 0; u < n; ++u) {
            if (color[u] >= 0) continue;

            color[u] = 0;
            queue<int> q; q.push(u);
            while (!q.empty()) {
                int u = q.front(); q.pop();
                for (int v : ady[u]) {
                    if (color[v] < 0) q.push(v),
                        color[v] = 1 - color[u];
                    if (color[u] == color[v])
                        return false;
                }
            }
        }
        return true;
    }

    // REVISIT
    // Obtiene el orden topologico de los nodos
    // de un grafo dirigido. Orden ascendente
    // respecto al numero de dependencias.

    vector<bool> vis;
    vector<int> ordenados;

    void OrdenTopologico(int u) {
        vis[u] = true;
        for (int v : ady[u])
            if (!vis[v]) OrdenTopologico(v);
        ordenados.push_back(u);
    }

    void OrdenTopologico() {
        ordenados.clear();
        vis = vector<bool>(n);
        for (int u = 0; u < n; ++u)
            if (!vis[u]) OrdenTopologico(u);
    }

    // Busqueda en amplitud desde el nodo s.
    // Devuelve el vector de distancias a todos
    // los nodos desde s. Un valor INF indica que
    // no es posible ir de s al respectivo nodo.

    vector<Costo> BFS(int s) {
        queue<int> q;
        vector<Costo> d(n, INF);
        q.push(s), d[s] = 0;

        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int v : ady[u])
                if (d[u] + 1 < d[v])
                    d[v] = d[u] + 1,
                    q.push(v);
        }
        return d;
    }
};

// Conjuntos disjuntos con Union-Find.

struct UnionFind {

    int n; vector<int> padre, tam;

    UnionFind(int N) : n(N),
        tam(N, 1), padre(N) {
        while (--N) padre[N] = N;
    }

    int Raiz(int u) {
        if (padre[u] == u) return u;
        return padre[u] = Raiz(padre[u]);
    }

    bool SonConexos(int u, int v) {
        return Raiz(u) == Raiz(v);
    }

    void Unir(int u, int v) {
        int Ru = Raiz(u);
        int Rv = Raiz(v);
        if (Ru == Rv) return;
        --n, padre[Ru] = Rv;
        tam[Rv] += tam[Ru];
    }

    int Tamano(int u) {
        return tam[Raiz(u)];
    }
};

typedef pair<Costo, int> CostoNodo;
typedef pair<Costo, Arista> Ponderada;

// Grafos con ponderacion.
// Nodos indexados de 0 a n - 1.
// GrafoPonderado(n,  true) -> Bidireccional.
// GrafoPonderado(n, false) -> Dirigido.

struct GrafoPonderado {

    int n; bool bi;
    vector<vector<CostoNodo>> ady;
    GrafoPonderado(int N, bool B = true)
        : n(N), bi(B), ady(N) {}

    void AgregarArista(int u, int v, Costo c) {
        if (bi) ady[v].push_back(CostoNodo(c, u));
        ady[u].push_back(CostoNodo(c, v));
    }

    // Obtiene el arbol de expansion minima de un
    // grafo bidireccional. Para obtener el arbol
    // de expansion maxima descomentar el reverse.
    // En caso de tener varias componentes conexas,
    // obtiene el bosque de expansion minima.

    vector<Ponderada> Kruskal() {
        vector<Ponderada> todas;
        for (int u = 0; u < n; ++u)
            for (CostoNodo arista : ady[u])
                todas.push_back(
                    Ponderada(arista.first,
                    Arista(u, arista.second)));
        sort(todas.begin(), todas.end());
        // reverse(todas.begin(), todas.end());

        vector<Ponderada> mst;
        UnionFind componentes(n);
        for (Ponderada arista : todas) {
            int u = arista.second.first;
            int v = arista.second.second;
            if (!componentes.SonConexos(u, v))
                componentes.Unir(u, v),
                mst.push_back(arista);
        }
        return mst;
    }

    // Algoritmo de dijkstra desde el nodo s.
    // Devuelve el vector de distancias a todos
    // los nodos desde s. Un valor INF indica que
    // no es posible ir de s al respectivo nodo.

    vector<Costo> Dijkstra(int s) {
        vector<Costo> dist(n, INF);
        priority_queue<CostoNodo> pq;
        pq.push(CostoNodo(0, s)), dist[s] = 0;

        while (!pq.empty()) {
            Costo p = -pq.top().first;
            int u = pq.top().second; pq.pop();
            if (dist[u] < p) continue;

            for (CostoNodo arista : ady[u]) {
                int v = arista.second;
                p = dist[u] + arista.first;
                if (p < dist[v]) dist[v] = p,
                    pq.push(CostoNodo(-p, v));
            }
        }
        return dist;
    }

    // Algoritmo de Bellman-Ford optimizado, desde
    // el nodo s. Devuelve un booleano indicando si
    // existe un ciclo negativo en un digrafo.
    // Obtiene el vector de distancias a todos.

    vector<Costo> dist; // <- Resultado

    bool BellmanFerrari(int s) {
        queue<int> q; 
        vector<bool> en_cola(n);
        vector<int> procesado(n);
        dist = vector<Costo>(n, INF);
        q.push(s), dist[s] = 0;

        while (!q.empty()) {
            int u = q.front();
            q.pop(), en_cola[u] = false;
            if (++procesado[u] == n) return true;
            for (CostoNodo arista : ady[u]) {
                int v = arista.second;
                Costo p = arista.first;
                if (dist[u] + p < dist[v]) {
                    if (!en_cola[v]) q.push(v);
                    dist[v] = dist[u] + p;
                    en_cola[v] = true;
                }
            }
        }
        return false;
    }
};

int main() {
    return 0;
}
