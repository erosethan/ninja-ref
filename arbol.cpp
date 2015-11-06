#include <bits/stdc++.h>
using namespace std;

// Definiciones iniciales.

typedef vector<int> Lista;

// Arbol con Heavy-Light Decomposition.
// Los nodos estan indexados de 0 a n - 1.

struct HeavyLight {

    int n, conteo;
    Lista nivel, tamano, up;
    Lista indice, super, top;
    vector<Lista> aristas;

    HeavyLight(int N) : n(N), conteo(),
        top(N), nivel(N), tamano(N), up(N),
        indice(N), super(N), aristas(N) {}

    void AgregarArista(int u, int v) {
        aristas[u].push_back(v);
        aristas[v].push_back(u);
    }

    void CalcularNivel(int u, int p) {
        for (int i = 0; i < aristas[u].size(); ++i) {
            int v = aristas[u][i]; if (p == v) continue;
            if (super[u] == super[v]) nivel[v] = nivel[u];
            else nivel[v] = nivel[u] + 1;
            CalcularNivel(v, u);
        }
    }

    // Construir realiza todas las operaciones para
    // trabajar con Heavy-Light. Por defecto, la raiz del
    // arbol se establece como el nodo 0. Si quieren definir
    // una raiz diferente, llamen Construir(r) donde el
    // parametro r indica cual sera la raiz del arbol.

    int Construir(int u = 0, int p = -1) {
        int tam_subarbol = 0;
        up[u] = p, super[u] = -1;

        for (int i = 0; i < aristas[u].size(); ++i) {
            int v = aristas[u][i]; if (p == v) continue;
            tam_subarbol += Construir(v, u);
        }
        for (int i = 0; i < aristas[u].size(); ++i) {
            int v = aristas[u][i]; if (p == v) continue;
            if (tamano[v] > tam_subarbol / 2)
                indice[u] = indice[v] + 1,
                super[u] = super[v],
                top[super[v]] = u;
        }
        if (super[u] == -1) super[u] = conteo,
                            top[conteo++] = u;
        if (p == -1) CalcularNivel(u, p);
        return tamano[u] = tam_subarbol + 1;
    }

    int LCA(int u, int v) {
        if (nivel[v] > nivel[u]) swap(u, v);
        while (nivel[u] > nivel[v]) u = up[top[super[u]]];
        while (super[u] != super[v]) u = up[top[super[u]]],
                                     v = up[top[super[v]]];
        return (indice[u] > indice[v])? u: v;
    }
};

int main() {
    return 0;
}
