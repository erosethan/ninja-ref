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
    //Igualdad de puntos 
    bool operator==(const Punto& cmp) const {
        return Igual(x, cmp.x) && Igual(y, cmp.y);
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
        if (!a) c = p.y, b = -1;
        else if (!b) c = p.x, a = -1;
        else {
            c = abs(__gcd(a, b));
            a /= c, b /= c;
            c = -a*p.x - b*p.y;
        }
    }

    bool operator<(const Linea& cmp) const {
        if (a != cmp.a) return a < cmp.a;
        if (b != cmp.b) return b < cmp.b;
        return c < cmp.c;
    }
};

// Rectas: r, s
// Segmentos: s, t

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
    return (l.a == m.b && l.b == -m.a) ||
           (m.a == l.b && m.b == -l.a);
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

// Saber si una recta r y un segmento s se intersectan
bool IntersecRectaSegmen(const Linea& r, const Linea& s) {
    return ManoDerecha(r.p, r.q, s.p) !=
           ManoDerecha(r.p, r.q, s.q);
}

// Saber si dos segmentos s y t se intersectan
bool InterseccionSegmentos(const Linea& s, const Linea& t) {
    if (ManoDerecha(s.p, s.q, t.p) ==
        ManoDerecha(s.p, s.q, t.q)) return false;
    if (ManoDerecha(t.p, t.q, s.p) ==
        ManoDerecha(t.p, t.q, s.q)) return false;
    return true;
}

// Obtener punto de interseccion entre lineas l y m
Punto PuntoInterseccion(const Linea& l, const Linea& m) {
    if (LineasParalelas(l, m)) return Punto();
    if (!l.a) return Punto((double)(l.c*m.b + m.c) / -m.a, l.c);
    double y = (double)(m.a*l.c - l.a*m.c) / (m.b*l.a - m.a*l.b);
    double x = (double)(l.c + l.b * y) / -l.a;
    return Punto(x, y);
}

// Obtener proyeccion del punto p en la recta r
Punto ProyeccionEnRecta(const Punto& p, const Linea& r) {
    Punto a = Trasladar(r.p, p), b = Trasladar(r.p, r.q);
    return Trasladar(Opuesto(r.p), Escalar(
        b, Dot(a, b) / pow(Magnitud(b), 2)));
}

// Distancia entre un punto p y una recta r
double DistanciaPuntoRecta(const Punto& p, const Linea& r) {
    return Distancia(ProyeccionEnRecta(p, r), p);
}

// Distancia entre un punto p y un segmento s
double DistanciaPuntoSegmento(const Punto& p, const Linea& s) {
    Punto proy = ProyeccionEnRecta(p, s);
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
    if (PuntoEnPerimetro(p, P)) return true;
    int dir = ManoDerecha(P[0], P[1], p);
    for (int i = 2; i < P.size(); ++i)
        if (dir != ManoDerecha(P[i - 1], P[i], p))
            return false;
    return true;
}

// Punto en poligono concavo por ray casting
bool RayCasting(const Punto& p, const Poligono& P) {
    if (PuntoEnPerimetro(p, P)) return true;
    Punto o = *min_element(P.begin(), P.end());
    Linea rayo(p, Punto(o.x - M_PI, o.y - M_E));
    int cruces = 0;
    for (int i = 1; i < P.size(); ++i)
        if (InterseccionSegmentos(rayo,
            Linea(P[i - 1], P[i]))) ++cruces;
    return cruces & 1;
}

// Punto en poligono concavo por angle summation
bool AngleSummation(const Punto& p, const Poligono& P) {
    if (PuntoEnPerimetro(p, P)) return true;
    double angulo = 0;
    for (int i = 1; i < P.size(); ++i)
        angulo += ManoDerecha(p, P[i - 1],P[i]) *
            Angulo(Trasladar(p, P[i - 1]), Trasladar(p, P[i]));
    return (fabs(angulo) > 180)? true: false;
}

// Area de un poligono
double Area(const Poligono& P) {
    double area = 0;
    for (int i = 1; i < P.size(); ++i)
        area += Cruz(P[i - 1], P[i]);
    return fabs(area) / 2.0;
}

// Perimetro de un poligono
double Perimetro(const Poligono& P) {
    double perimetro = 0;
    for (int i = 1; i < P.size(); ++i)
        perimetro += Distancia(P[i - 1], P[i]);
    return perimetro;
}

// Cerco convexo de un conjunto de puntos
Poligono CercoConvexo(vector<Punto> P){
    // Si P ya esta ordenado no usar sort
    sort(P.begin(), P.end());
    Poligono arriba, abajo;
    for (int i = 0; i < P.size(); ++i) {
        while (arriba.size() > 1) {
            int p = arriba.size() - 1;
            // Permitir colineales: <=
            if (ManoDerecha(arriba[p - 1],
                arriba[p], P[i]) < 0) break;
            arriba.pop_back();
        }
        arriba.push_back(P[i]);
    }
    arriba.pop_back();
    for (int i = P.size() - 1; i >= 0; --i) {
        while (abajo.size() > 1) {
            int p = abajo.size() - 1;
            // Permitir colineales: <=
            if(ManoDerecha(abajo[p - 1],
                abajo[p], P[i]) < 0) break;
            abajo.pop_back();
        }
        abajo.push_back(P[i]);
    }
    arriba.insert(arriba.end(),
        abajo.begin(), abajo.end());
    return arriba;
}

Punto Centroide(const Poligono& P) {
    double x = 0, y = 0, k = 0;
    for (int i = 1; i < P.size(); ++i) {
        double cruz = Cruz(P[i - 1], P[i]);
        x += cruz * (P[i - 1].x + P[i].x);
        y += cruz * (P[i - 1].y + P[i].y);
        k += cruz * 3;
    }
    return Punto(x / k, y / k);
}

bool RectaCortaPoligono( const Linea& r, const Poligono& P) {
    int dif = 0;
    for(int i = 0; i < P.size(); ++i) {
        int a = ManoDerecha(l.p, l.q, P[i] );
        if(dif == 0) dif = a;
        if(a != dif && a && dif)
            return true;
    }
    return false;
}

vector<Poligono> CortarPoligono(
    const Poligono& P, const Linea& r) {
    if (!RectaCortaPoligono(r, P))
        return vector<Poligono>(1, P);
    int ind = 0;
    vector<Poligono> Ps(2);
    for (int i = 1; i < P.size(); ++i) {
        Linea s(P[i - 1], P[i]);
        if (IntersecRectaSegmen(r, s)) {
            Punto p = PuntoInterseccion(r, s);
            if (P[i - 1] == p) continue;
            Ps[ind].push_back(P[i - 1]);
            Ps[1 - ind].push_back(p);
            Ps[ind].push_back(p);
            ind = 1 - ind;
        }
        else Ps[ind].push_back(P[i - 1]);
    }
    Ps[0].push_back(Ps[0][0]);
    Ps[1].push_back(Ps[1][0]);
    return Ps;
}


int main() {
    return 0;
}
