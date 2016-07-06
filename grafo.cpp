#include <bits/stdc++.h>
using namespace std;

typedef int Costo;
typedef pair<int, int> Arista;

const Costo INF = 1 << 30;

//

struct Grafo {

	int n; bool bi;
	vector<vector<int>> ady;
	Grafo(int N, bool B = true)
		: n(N), bi(B), ady(N) {}

	void AgregarArista(int u, int v) {
		if (bi) ady[v].push_back(u);
		ady[u].push_back(v);
	}

	//

	vector<int> ciclo;
	vector<char> color;

	void DetectarCiclo(int u, int p) {
		color[u] = ciclo.empty()? 'G': 'N';
		for (int v : ady[u]) {
			if (bi && v == p) continue;
			if (ciclo.empty() && color[v] == 'G')
				color[v] = 'A', ciclo.push_back(v),
				color[u] = 'R', ciclo.push_back(u);
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
			ciclo.clear(); DetectarCiclo(u, u);
			reverse(ciclo.begin(), ciclo.end());
			ciclos.push_back(ciclo);
		}
		return ciclos;
	}

	//

	int tiempo;
	vector<int> label, low;
	vector<Arista> puentes;
	vector<bool> articulacion;

	int PuentesArticulacion(int u, int p) {
		label[u] = low[u] = ++tiempo;
		int hijos = 0;
		for (int v : ady[u]) {
			if (v == p) continue;
			if (!label[v]) { ++hijos;
				PuentesArticulacion(v, u);
				if (label[u] < low[v])
					puentes.push_back(Arista(u, v));
				if (label[u] <= low[v])
					articulacion[u] = true;
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
				PuentesArticulacion(u, u) > 1;
	}

	//

	vector<vector<int>> scc;
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

	//

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

	//

	vector<Costo> BFS(int s) {
		queue<int> q;
	    vector<Costo> d(n, INF);
	    d[s] = 0, q.push(s);

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

//

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

//

struct GrafoPonderado {

	int n; bool bi;
	vector<vector<CostoNodo>> ady;
	GrafoPonderado(int N, bool B = true)
		: n(N), bi(B), ady(N) {}

	void AgregarArista(int u, int v, Costo c) {
		if (bi) ady[v].push_back(CostoNodo(c, u));
		ady[u].push_back(CostoNodo(c, v));
	}

	//

	vector<Ponderada> Kruskal() {
		vector<Ponderada> todas;
		for (int u = 0; u < n; ++u)
			for (CostoNodo arista : ady[u])
				todas.push_back(
					Ponderada(arista.second,
					Arista(u, arista.first)));
		sort(todas.begin(), todas.end());

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

	//

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
					pq.push(CostoNodo(p, v));
			}
		}
		return dist;
	}
};
