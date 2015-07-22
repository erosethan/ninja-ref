#include <math.h>
#include <bits/stdc++.h>
using namespace std;
const double ERR=1e-9;
typedef double coord;
struct Punto{
	//2d
	coord x,y;
	Punto():x(0),y(0){}
	Punto(coord x_,coord y_)
		:x(x_),y(y_){}
	bool operator<(const Punto&cmp)const{
		if(x==cmp.x)return y<cmp.y;
		return x<cmp.x;
	}
};
//puntos p,q 
//vectores v,w

typedef vector<Punto> Poligono;
//p1,p2,p3,.....,pn,p1;

bool Igual(coord a, coord b){
	//ERR margen de error permitido 
	return fabs(a-b)<ERR;
}

//Distancia entre dos puntos
double Distancia( const Punto& p, const Punto& q){
	return hypot(p.x-q.x,p.y-q.y);
}
//Magnitud de un vector
double Magnitud(const Punto& v){
	return hypot(v.x,v.y);
}

//producto punto de vectores
double Dot(const Punto& v,const Punto& w){
	return v.x*w.y+v.y*w.x;
}

//producto cruz regresa magnitud
double Cruz(const Punto& v,const Punto& w){
	return (v.x*w.y)-(v.y*w.x);
}

//conversion de Grados a radianes 
double GradARad(double grados){
	return (grados*M_PI)/180.0;
}

//conversion de Radianes a grados
double RadAGrad(double radianes){
	return (radianes*180.0)/M_PI;
}
//rotar un punto en ccw
//rotar en cw 360-alpha
Punto Rotar(const Punto& p,double grados){
	double r=GradARad(grados);
	return Punto(p.x*cos(r)-p.y*sin(r),p.x*sin(r)+p.y*cos(r));
}

//Trasladar un punto tomando a o como el origen 
Punto Trasladar(const Punto& o,const Punto& p){
	return Punto(p.x-o.x,p.y-o.y);
}

//vector opuesto
Punto Opuesto(const Punto& v){
	return Punto(-v.x,-v.y);
}

double Angulo(const Punto& v, const Punto& w){
	//regreso en radianes
	double r=acos(Dot(v,w)/(Magnitud(v)*Magnitud(w)));
	//convertimos a grados
	return RadAGrad(r);
}
Punto Escalar(const Punto& v, double s){
	return Punto(v.x*s,v.y*s);
}

//1 izq 0 colineales -1 derecha
int ManoDerecha(const Punto& o,const Punto& v,const Punto& w){
	double ccw=Cruz(Trasladar(o,v),Trasladar(o,w));
	if(Igual(ccw,0))return 0;
	return ccw<0?-1:1;
}

bool PuntoEnPConvexo(const Poligono& P,const Punto& a){
	int dir=ManoDerecha(P[0],P[1],a);
	int tam=P.size()-1;
	for(int i=1;i<tam;i++)
		if((abs(ManoDerecha(P[i],P[i+1],a))-dir)==2)
			return false;
	return true;
}

bool PuntoEnConcavo(const Poligono& P,const Punto& p){
	double angulo=0;
	int tam=P.size()-1;
	for(int i=0;i<tam;i++)
		angulo+=Angulo(Trasladar(p,P[i]),Trasladar(p,P[i+1]))*ManoDerecha(p,P[i],P[i+1]);
	return (angulo>180)?true:false;
}

struct Linea{
	Punto p,q;
	long long a,b,c;
	Linea():p(),q(),a(0),b(1),c(0){}
	Linea(coord a_,coord b_,coord c_):p(),q(),a(a_),b(b_),c(c_){}
	Linea(const Punto& p_,const Punto& q_):p(p_),q(q_),a(),b(),c(){
		//solo para enteros
		if(q<p)swap(p,q);//punto p como el menor
		a=q.y-p.y;
		b=p.x-q.x;
		c=__gcd(a,b);
		a/=c;
		b/=c;
		if(a==0)c=-p.y;
		else if(b==0)c=-p.x;
		else c=Cruz(q,p);
	}
};
bool PuntoenRecta(const Linea& r,const Punto& p){
	return Igual(p.x*r.a+p.y*r.b+r.c,0);
}
bool PuntoEnSegmento(const Punto&p,const Linea& s){
	return PuntoenRecta(s,p)&&!(p<s.p||s.q<p);
}

bool PuntoEnPerimetro(const Poligono& P,const Punto& p){
	int tam=P.size();
	for(int i=1;i<tam;i++){
		Punto p1=P[i-1],p2=P[i];
		if(p2<p1)swap(p1,p2);
		if(!ManoDerecha(p1,p2,p)&&!(p<p1||p2<p))return true;
		return false;
	}
	/*mayor costo gcd en linea
	for(int i=1;i<tam;i++)
		if(PuntoEnSegmento(p,Linea(P[i-1],P[i])))return true;
	return false;*/
}
bool Parelelas(const Linea& l1,const Linea& l2){
	return (l1.a==l2.a&&l2.b==l1.b);
}
bool lineasiguales(const Linea& l1, const Linea&l2){
	return (Igual(l1.a/l1.b,l2.a/l2.b)&&Igual(l1.c/l1.b,l2.c/l2.b));
}
bool lineaperpendicular(){}
bool lineaperpendicularEnpunto(){}

int main(){
	Linea a(Punto(-2,1),Punto(2,3));
	cout<<a.a<<' '<<a.b<<' '<<a.c<<endl;
	Linea b(Punto(-2,-1),Punto(2,1));
	cout<<b.a<<' '<<b.b<<' '<<b.c<<endl;
	cout<<Parelelas(a,b)<<endl;
}
