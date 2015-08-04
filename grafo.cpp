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

int main() {
	return 0;
}