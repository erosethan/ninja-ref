#include <bits/stdc++.h>
using namespace std;

// Definiciones iniciales.

typedef long long Long;
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
            double ratio = m[j][i] / m[i][i];
            for (int k = i; k < m[j].size(); ++k)
                m[j][k] -= m[i][k] * ratio;
        }
    }
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

//Factores Primos

#include <cmath>
typedef pair<int, int> pii;

vector<pii> factors(int a){
	int count = 0, root;
	while(!(a%2)){
		a /= 2;
		count++;
	}
	vector<pii> facts;
	if(count) facts.push_back(pii(2,count));
	root = sqrt(a);
	for(int i = 3; i <= root; i++){
		count = 0;
		while(!(a%i)){
			a /= i;
			count ++;
		}
		if(count) facts.push_back(pii(i,count));
	}
	if(a > 1) facts.push_back(pii(a,1));
	return facts;
}

//Sieve
const int MAXN = ;

bitset<MAXN + 1> primes;
void Sieve() {
	primes.set();
	primes[0] = primes[1] = 0;
	long long sroot = sqrt(MAXN);
	for (long long j = 2LL * 2LL; j <= MAXN; j += 2LL)
		primes[j] = 0;
	long long sroot = sqrt(MAXN);
	for (long long i = 3LL; i <= sroot; i += 2LL) {
		if (primes[i]) {
			for (long long j = i * i; j <= MAXN; j += i)
				primes[j] = 0;
		} 
	}
	return;
}

long long sieve[MAXN + 1];
void Sieve() {
	for (long long j = 2LL * 2LL; j <= MAXN; j += 2LL)
		sieve[j] = 2LL;
	long long sroot = sqrt(MAXN);
	for (long long i = 3LL; i <= sroot; i += 2LL) {
		if (sieve[i]) {
			for (long long j = i * i; j <= MAXN; j += i)
			if (!sieve[j])
				sieve[j] = i;
		} 
	}
	return;
}

int exp[MAXN + 1];
void AcumuladoFactorial() {
	fill(exp, exp + MAXN + 1, 1);
	for (int i = n; i > 3; --i) {
		if (!sieve[i]) continue;
		exp[sieve[i]] += exp[i];
		exp[i / sieve[i]] += exp[i];
		exp[i] = 0;
	}
	return;
}


#define ll long long
using namespace std;

ll gcd(ll a, ll b, ll* x, ll* y, ll mod){
    if(b==0){
        if (x && y) {
            *x=1;
            *y=0;
        }
        return a;
    }
    ll GCD = gcd(b, a%b, x, y, mod);
    if (x && y) {
        ll x_1, y_1;
        x_1 = *y;
        y_1 = (*x - ((*y)*(a/b))%mod + mod)%mod;
        *x = x_1;
        *y = y_1;
    }
    return GCD;
}

int main(){
    ll a, b, x, y, GCD;
    cin >> a >> b;
    GCD = gcd(a, b, &x, &y, b);
    cout << "Inverso de : " << a << " Modulo: " << b << " Es: " << x;
}



int main() {
    return 0;
}
