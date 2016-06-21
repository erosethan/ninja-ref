#include <bits/stdc++.h>
using namespace std;

// Definiciones iniciales.

typedef long long Long;
typedef pair<int, int> Factor;

vector<Factor> FactoresPrimos(int a) {
	int conteo = 0;
	vector<Factor> factores;
	while (!(a % 2)) ++conteo, a /= 2;
	if (conteo) factores.push_back(
		Factor(2, conteo)), conteo = 0;

	int raiz = sqrt(a);
	for (int i = 3; i <= raiz; i += 2) {
		while (!(a % i)) ++conteo, a /= i;
		if (conteo) factores.push_back(
			Factor(i, conteo)), conteo = 0;
	}
	if (a > 1) factores.push_back(
		Factor(a, 1));
	return factores;
}

vector<int> Criba(int n) {
	int raiz = sqrt(n); vector<int> criba(n + 1);
	for (int i = 4; i <= n; i += 2) criba[i] = 2;
	for (int i = 3; i <= raiz; i += 2) if (!criba[i])
		for (int j = i * i; j <= n; j += i)
			if (!criba[j]) criba[j] = i;
	return criba;
}

vector<Factor> FactoresFactorial(
	int n, const vector<int>& primos) {

	vector<Factor> factores(primos.size());
	for (int i = 0; i < primos.size(); ++i) {
		int p = primos[i], reps = n / p;
		while (primos[i] <= n / p)
			p *= primos[i], reps += n / p;
		factores[i] = Factor(primos[i], reps);
	}
	return factores;
}

// Exponenciacion binaria a^n mod m.

Long Exponenciar(Long a, Long n, Long m) {
    Long res = 1, p = a;
    for (; n; n >>= 1) {
        if (n & 1) res =
            (res * p) % m;
        p = (p * p) % m;
    }
    return res;
}

// Multiplicacion binaria a*b mod m.

Long Multiplicar(Long a, Long b, Long m) {
    Long res = 0, p = a;
    for (; b; b >>= 1) {
        if (b & 1) res =
            (res + p) % m;
        p = (p + p) % m;
    }
    return res;
}

Long Euclides(Long a, Long b,
    Long& x, Long& y, Long mod = 0) {
    if (!b) { x = 1, y = 0; return a; }
    Long gcd = Euclides(b, a % b, x, y, mod);

    x = !mod? x - y * (a / b):
    	(mod + x - y * (a / b)) % mod;
    swap(x, y); return gcd;
}

// Definiciones para eliminacion Gaussiana.

typedef vector<double> Vector;
typedef vector<Vector> Matriz;

const double ERROR = 1e-9;
const double M_2PI = 2 * M_PI;

// Tolerancia en flotantes.

bool Igual(double a, double b) {
    return fabs(a - b) < ERROR;
}

// Realiza eliminacion Gaussiana en una matriz.

void EliminaGaussiana(Matriz& m) {
    for (int i = 0; i < m.size(); ++i) {
        int fila_mayor = i;
        for (int j = i + 1; j < m.size(); ++j)
            if (fabs(m[fila_mayor][i]) <
                fabs(m[j][i])) fila_mayor = j;
        swap(m[i], m[fila_mayor]);
        if (Igual(m[i][i], 0)) continue;
        for (int j = m[i].size() - 1; j >= i; --j)
            m[i][j] /= m[i][i];
        for (int j = 0; j < m.size(); ++j) {
            if (i == j || Igual(m[j][i], 0)) continue;
            for (int k = m[j].size()-1; k >= i; --k)
                m[j][k] -= m[i][k] * m[j][i];
        }
    }
}

// Tipo de dato para operar con numeros complejos.

struct Complejo {
    double real, imag;
    Complejo() : real(), imag() {}
    Complejo(double R, double I) : real(R), imag(I) {}

    Complejo operator+(const Complejo& c) {
        return Complejo(real + c.real, imag + c.imag); }
    Complejo operator-(const Complejo& c) {
        return Complejo(real - c.real, imag - c.imag); }
    Complejo operator*(const Complejo& c) {
        return Complejo(real * c.real - imag * c.imag,
                        real * c.imag + imag * c.real); }
};

// Transformada rapida de Fourier.
// Se tiene que garantizar que el numero de
// elementos en el vector sea una potencia de 2.

vector<Complejo> FastAndFourier(
    const vector<Complejo>& a, int k = 1) {

    int n = a.size();
    if (n == 1) return a;
    vector<Complejo> par, impar;
    for (int i = 0; i < n; ++i)
        if (i & 1) par.push_back(a[i]);
        else impar.push_back(a[i]);
    impar = FastAndFourier(impar, k);
    par = FastAndFourier(par, k);

    vector<Complejo> fourier(n);
    Complejo w(1, 0), wn(cos(-k * M_2PI/n),
                         sin(-k * M_2PI/n));
    for (int i = 0; i < n/2; w = w * wn, ++i) {
        fourier[i + n/2] = impar[i] - w * par[i];
        fourier[i] = impar[i] + w * par[i];
    }
    return fourier;
}

// Transformada inversa de Fourier.
// Se tiene que garantizar que el numero de
// elementos en el vector sea una potencia de 2.

vector<Complejo> InvFastAndFourier(
    const vector<Complejo>& a) {

    vector<Complejo> ifft = FastAndFourier(a, -1);
    for (int i = 0; i < ifft.size(); ++i)
        ifft[i].real /= ifft.size(),
        ifft[i].imag /= ifft.size();
    return ifft;
}

// Convolucion discreta de dos vectores usando
// transformada rapida de Fourier O(n log n).
// Multiplica eficientemente dos polinomios.

Vector ConvolucionDiscreta(
    const Vector& a, const Vector& b) {

    int n = a.size() + b.size() - 1;
    int p = 1; while (p < n) p <<= 1;

    vector<Complejo> A(p), B(p), C(p);
    for (int i = 0; i < a.size(); ++i)
        A[i] = Complejo(a[i], 0);
    for (int i = 0; i < b.size(); ++i)
        B[i] = Complejo(b[i], 0);

    A = FastAndFourier(A);
    B = FastAndFourier(B);
    for (int i = 0; i < p; ++i)
        C[i] = A[i] * B[i];
    C = InvFastAndFourier(C);

    Vector convolucion(n);
    for (int i = 0; i < n; ++i)
        convolucion[i] = C[i].real;
    return convolucion;
}

int main() {
    return 0;
}
