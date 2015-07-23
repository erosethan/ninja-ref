#include <bits/stdc++.h>
using namespace std;

// Definiciones iniciales
typedef double Coord;
const double ERROR = 1e-9;

// Tolerancia en flotantes
bool Igual(Coord a, Coord b) { 
    return fabs(a - b) < ERROR;
}

// Punto en 2D
struct Punto {
    Coord x, y;

    Punto() : x(0), y(0) {}
    Punto(Coord x_, Coord y_)
        : x(x_), y(y_) {}

    // Izquierda a derecha, abajo a arriba
    bool operator<(const Punto& cmp) const {
        if (Igual(x, cmp.x))
            return y < cmp.y;
        return x < cmp.x;
    }
};

// Puntos: p, q
// Vectores: v, w

// Distancia entre dos puntos p y q
double Distancia(const Punto& p, const Punto& q) {
    return hypot(p.x - q.x, p.y - q.y);
}

// Magnitud de un vector v
double Magnitud(const Punto& v) {
    return hypot(v.x, v.y);
}

// Producto punto entre vectores v y w
double Dot(const Punto& v, const Punto& w) {
    return v.x * w.x + v.y * w.y;
}

// Producto cruz entre vectores v y w
double Cruz(const Punto& v, const Punto& w) {
    return v.x * w.y - v.y * w.x;
}

// Conversion de grados a radianes
double GradARad(double grd) {
    return (grd * M_PI) / 180.0;
}

// Conversion de radianes a grados
double RadAGrad(double rad) {
    return (rad * 180.0) / M_PI;
}

// Rotar un punto respecto al origen
// La rotación se hace en orden CCW, para
// rotar en CW llamar Rotar(p, 360 - alpha)
Punto Rotar(const Punto& p, double alpha) {
    double rad = GradARad(alpha);
    return Punto(p.x*cos(rad) - p.y*sin(rad),
                 p.x*sin(rad) + p.y*cos(rad));
}

// Trasladar p tomando como origen al punto o 
Punto Trasladar(const Punto& o, const Punto& p) {
    return Punto(p.x - o.x, p.y - o.y);
}

// Escalar un vector v por un factor s
Punto Escalar(const Punto& v, double s) {
    return Punto(v.x * s, v.y * s);
}

// Obtener vector opuesto
Punto Opuesto(const Punto& v) {
    return Punto(-v.x, -v.y);
}

// Angulo entre vectores v y w
double Angulo(const Punto& v, const Punto& w) {
    return RadAGrad(acos(Dot(v, w) / (Magnitud(v) * Magnitud(w))));
}

// Test de mano derecha: CCW = 1, CW = -1, Colineal = 0
int ManoDerecha(const Punto& o, const Punto& p, const Punto& q) {
    double ccw = Cruz(Trasladar(o, p), Trasladar(o, q));
    return Igual(ccw, 0)? 0: (ccw < 0)? -1: 1;
}

// Linea en 2D
struct Linea {
    Punto p, q;
    long long a, b, c;

    Linea() : p(), q(), a(0), b(0), c(0) {}
    Linea(Coord a_, Coord b_, Coord c_)
        : p(), q(), a(a_), b(b_), c(c_) {}

    Linea(const Punto& p_, const Punto& q_)
        : p(p_), q(q_), a(), b(), c() {
        // Asegura p como punto menor
        if (q < p) swap(p, q);
        a = q.y - p.y;
        b = p.x - q.x;
        if (!a) c = -p.y, b = 1;
        else if (!b) c = -p.x, a = 1;
        else {
            c = abs(__gcd(a, b));
            if (a < 0) c = -c;
            a /= c, b /= c;
            c = Cruz(q, p);
        }
    }
};

// Saber si un punto p esta en la recta r
bool PuntoEnRecta(const Punto& p, const Linea& r) {
    return Igual(r.a*p.x + r.b*p.y + r.c, 0);
}

// Saber si un punto o esta en el segmento s
bool PuntoEnSegmento(const Punto& p, const Linea& s) {
    return PuntoEnRecta(p, s) && !(p < s.p || s.q < p);
}

// Saber si dos lineas son paralelas
bool LineasParalelas(const Linea& l, const Linea& m) {
    return l.a == m.a && l.b == m.b;
}

// Saber si dos lineas son iguales
bool LineasIguales(const Linea& l, const Linea& m) {
    return LineasParalelas(l, m) && l.c == m.c;
}

// Saber si dos lineas son perpendiculares
bool LineasPerpendiculares(const Linea& l, const Linea& m) {
    return l.a == m.b && l.b == m.a;
}

// Obtener una linea paralela a l que pase por p
Linea ParalelaEnPunto(const Linea& l, const Punto& p) {
    return Linea(p, Punto(p.x - l.b, p.y + l.a));
}

// Obtener una linea perpendicular a l que pase por p
Linea PerpendicularEnPunto(const Linea& l, const Punto& p) {
    return Linea(p, Punto(p.x + l.a, p.y + l.b));
}

// Saber si dos rectas r y s se intersectan
bool InterseccionRectas(const Linea& r, const Linea& s) {
    return !LineasParalelas(r, s);
}

// Saber si dos segmentos s y t se intersectan
bool InterseccionSegmentos(const Linea& s, const Linea& t) {
    if (ManoDerecha(s.p, s.q, t.p) ==
        ManoDerecha(s.p, s.q, t.q)) return false;
    if (ManoDerecha(t.p, t.q, s.p) ==
        ManoDerecha(t.p, t.q, s.p)) return false;
    return true;
}

// Obtener punto de interseccion entre lineas l y m
Punto PuntoInterseccion(const Linea& l, const Linea& m) {
    if (LineasParalelas(l, m)) return Punto();
    if (!l.a) return Punto((double)(-l.c*m.b - m.c) / m.a, -l.c);
    double y = (double)(m.a*l.c - l.a*m.c) / (m.b*l.a - m.a*l.b);
    double x = (double)(l.c + l.b * y) / -l.a;
    return Punto(x, y);
}

// Obtener proyeccion del punto p en la recta r
Punto ProyeccionEnRecta(const Linea& r, const Punto& p) {
    Punto a = Trasladar(r.p, p), b = Trasladar(r.p, r.q);
    return Trasladar(Opuesto(r.p), Escalar(
        b, Dot(a, b) / pow(Magnitud(b), 2)));
}

// Distancia entre un punto p y una recta r
double DistanciaPuntoRecta(const Punto& p, const Linea& r) {
    return Distancia(ProyeccionEnRecta(r, p), p);
}

// Distancia entre un punto p y un segmento s
double DistanciaPuntoSegmento(const Punto& p, const Linea& s) {
    Punto proy = ProyeccionEnRecta(s, p);
    if (proy < s.p) return Distancia(s.p, p);
    if (s.q < proy) return Distancia(s.q, p);
    return Distancia(proy, p);
}

// Distancia entre dos lineas l y m
double DistanciaRectaRecta(const Linea& l, const Linea& m) {
    return LineasParalelas(l, m)? DistanciaPuntoRecta(l.p, m): 0;
}

// Un poligono es una serie de
// vertices conectados por aristas
// P = p1, p2, p3, ..., pn, p1
typedef vector<Punto> Poligono;

// Saber si un punto esta en el perimetro de un poligono
bool PuntoEnPerimetro(const Punto& p, const Poligono& P) {
    for (int i = 1; i < P.size(); ++i) {
        Punto l = P[i - 1], r = P[i];
        if (r < l) swap(l, r);
        if (!ManoDerecha(l, r, p) &&
            !(p < l || r < p)) return true;
    }
    return false;
}

// Saber si un punto esta dentro de un poligono convexo
bool PuntoEnConvexo(const Punto& p, const Poligono& P) {
    int dir = ManoDerecha(P[0], P[1], p);
    for (int i = 2; i < P.size(); ++i)
        if (2 == abs(dir - ManoDerecha(
            P[i - 1], P[i], p))) return false;
    return true;
}

// Saber si un punto esta dentro de un poligono concavo
bool PuntoEnConcavo(const Punto& p, const Poligono& P) {
    double angulo=0;
    int tam=P.size()-1;
    for(int i=0;i<tam;i++)
        angulo+=Angulo(Trasladar(p,P[i]),Trasladar(p,P[i+1]))*ManoDerecha(p,P[i],P[i+1]);
    return (angulo>180)?true:false;
}

int main() {
}
