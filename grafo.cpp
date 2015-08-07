#include <bits/stdc++.h>
using namespace std;

const int MAXN = 100000;
const int MAXM = 100000;

typedef int Arista;
typedef pair<int, int> AristaPeso;
typedef vector<Arista> Nodo;

Nodo grafo[MAXN];

int in_degree[MAXN];

void AgregarArista(int u, int v) {
	grafo[u].push_back(v);
	++in_degree[v];
}

int top_ciclo;
int ciclo[MAXN];
bool ciclo_activo;
bool visitado[MAXN];

int EncontrarCiclo_(int u, int p) {
	visitado[u] = true;
	for(int i = 0; i < grafo[u].size(); ++i) {
		int v = grafo[u][i];
		if (v == p) continue;
		if (visitado[v]) return ciclo[top_ciclo++] = v; 
		if ((int cic = EncontrarCiclo_(v, u)) >= 0) {
			if (ciclo_activo) ciclo[top_ciclo++] = v;
			if (cic == u) ciclo_activo = false;
			return cic;
		}
	}
	return -1;
}

int EncontrarCiclo(int u){
	ciclo_activo = true, top_ciclo = 0;
	return EncontrarCiclo_(u, -1);
}

//menor a mayor dependencia
//indexado de 0 a n-1
vector <int> ordentopo;

void OrdenTopologico_(int u) {
	visitado[u] = true;
	for(int i = 0; i < grafo[u].size(); ++i) {
		int v =  grafo[u][i];
		if(!visitado[v])
			OrdenTopologico_(v);
	}
	ordentopo.push_back(u);
}
void OrdenTopologico(int n){
	fill(visitado, visitado+n, false);
	ordentopo.clear();
	for(int i = 0; i < n; i++) {
		if(!visitado[i])
			OrdenTopologico_(i);
	}
}

int numeracion;
int top_activo;
int low[MAXN];
int num[MAXN];
int activos[MAXN];

void ObtenerCFC_(int u, vector <vector<int> >& componentes) {
	low[u]=num[u]=++numeracion;
	activo[top_activo++] = u;
	for(int i = 0; i < grafo[u].size(); ++i) {
		int v = grafo[u][i];
		if(!num[v])
			ObtenerCFC_(v, componentes);
		low[u] = min(low[u], low[v]);
	}
	if(low[u] == num[u]) {
		vector<int> cFC;
		while( activo[top_activo - 1] != u ) {
			cFC.push_back(activo[--top_activo]);
			low[activo[top_activo]] = INF;
		}
		cFC.push_back(activo[--top_activo]);
		low[activo[u]] = INF;
		componentes.push_back(cFC);
	}
}
vector <vector<int> > ObternerCFC(int n) {
	vector <vector<int> > cfcs;
	fill(num, num + n, 0);
	fill(low, low + n, 0);
	numeracion = 0;
	for(int i = 0; i < n; ++i) {
		if(!num[i]) ObternerCFC_(i, cfcs);
	}
	return cfcs;
}

Nodo puentes[MAXN];
bool punto_art[MAXN];

void PuntosArtPuntes_(int u, int p) {
 	low[u] = num[u] = ++numeracion;
 	for(int i = 0; i < grafo[u].size(); i++) {
 		int v = grafo[u][i];
 		if(v == p) continue;
 		if(!num[v]) PuntosArtPuntes_(v, u);
 		if(low[v] > num[u]) {
 			puentes[u].push_back(v);
 			puentes[v].push_back(u);
 		}
 		punto_art[u] |=low[v] >= num[u];
 		low[u] = min(low[u], low[v]);
 	}
 	low[u] = INF;
}

void PuntosArtPuntes(int n) {
	fill(num, num + n, 0);
	fill(low, low + n, 0);
	numeracion = 0;
	fill(punto_art, punto_art + n, false);
	for(int i = 0; i < n; i++) puentes[i].clear();
	for(int i = 0; i < n; i++) {
		if(num[i]) continue;
		PuntosArtPuntes_(i, -1);
		punto_art[i]= grafo[i].size() > 1;
	}
}


int main() {
	return 0;
}