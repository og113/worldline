/*
	definitions for program to generate unit loops in D dimensions for use in worldline programs.
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_monte.h>
//#include <gsl/gsl_monte_plain.h>
//#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include "simple.h"
#include "genloop.h"
#include "genfloop.h"

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	contents:
		1 - FCoeff class
		2 - functions acting on FCoeffs
		3 - FLoop class
		4 - FLoop functions
		5 - Metropolis
		6 - explicit template instantiation
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	1 - FCoeff class
		- initializers and destructors
		- indexing
		- algebraic operations
----------------------------------------------------------------------------------------------------------------------------*/

// initializers and destructors
template<uint Dim>
FCoeff<Dim>::FCoeff(): Coeffs(2*Dim,0.0) {}

template <uint Dim>
FCoeff<Dim>::FCoeff(const FCoeff<Dim>& p): Coeffs(p.Coeffs) {}

template <uint Dim>
FCoeff<Dim>::FCoeff(const vector<number>& c): Coeffs(2*Dim) {
	if (c.size()==2*Dim)
		Coeffs = c;
	else
		cerr << "FCoeff error: cannot initialize Dim=" << Dim << " FCoeff from vector of size " << c.size() << endl;
}

template<uint Dim>
void FCoeff<Dim>::copy(const FCoeff<Dim>& p) {
	Coeffs = p.Coeffs;
}

template <uint Dim>
FCoeff<Dim>& FCoeff<Dim>::operator=(const FCoeff<Dim>& p) {
	copy(p);
	return *this;
}

template <uint Dim>
void FCoeff<Dim>::operator()(const FCoeff& p) {
	copy(p);
}

// set to zero
template<uint Dim>
void FCoeff<Dim>::zero() {
	fill(Coeffs.begin(), Coeffs.end(), 0.0);
}

// indexing
template <uint Dim>
const number& FCoeff<Dim>::operator[](const uint& loc) const {
	return Coeffs.at(loc);
}

template <uint Dim>
number& FCoeff<Dim>::operator[](const uint& loc) {
	return Coeffs.at(loc);
}

// +=
template <uint Dim>
FCoeff<Dim>& FCoeff<Dim>::operator+=(const FCoeff<Dim>& rhs) {
	for (uint j=0; j<2*Dim; j++)
		Coeffs[j] += rhs.Coeffs[j];
	return *this;
}

// -=
template <uint Dim>
FCoeff<Dim>& FCoeff<Dim>::operator-=(const FCoeff<Dim>& rhs) {
	for (uint j=0; j<2*Dim; j++)
		Coeffs[j] -= rhs.Coeffs[j];
	return *this;
}

// *=
template <uint Dim>
FCoeff<Dim>& FCoeff<Dim>::operator*=(const number& rhs) {
	for (uint j=0; j<2*Dim; j++)
		Coeffs[j] *= rhs;
	return *this;
}

// /=
template <uint Dim>
FCoeff<Dim>& FCoeff<Dim>::operator/=(const number& rhs) {
	for (uint j=0; j<2*Dim; j++)
		Coeffs[j] /= rhs;
	return *this;
}

// operator +
template <uint Dim>
FCoeff<Dim> operator+(const FCoeff<Dim>& lhs,const FCoeff<Dim>& rhs) {
	FCoeff<Dim> p(lhs);
	for (uint j=0; j<2*Dim; j++)
		p.Coeffs[j] += rhs.Coeffs[j];
	return p;
}

// operator -
template <uint Dim>
FCoeff<Dim> operator-(const FCoeff<Dim>& lhs,const FCoeff<Dim>& rhs) {
	FCoeff<Dim> p(lhs);
	for (uint j=0; j<2*Dim; j++)
		p.Coeffs[j] -= rhs.Coeffs[j];
	return p;
}

// operator *
template <uint Dim>
FCoeff<Dim> operator*(const number& lhs,const FCoeff<Dim>& rhs) {
	FCoeff<Dim> p(rhs);
	for (uint j=0; j<2*Dim; j++)
		p.Coeffs[j] *= lhs;
	return p;
}

// operator /
template <uint Dim>
FCoeff<Dim> operator/(const FCoeff<Dim>& lhs,const number& rhs) {
	FCoeff<Dim> p(lhs);
	for (uint j=0; j<2*Dim; j++)
		p.Coeffs[j] /= rhs;
	return p;
}


// operator<<
template <uint Dim>
ostream& operator<<(ostream& os,const FCoeff<Dim>& p) {
	os << left << setprecision(16);
	for (uint j=0; j<2*Dim; j++)
		os << setw(24) << p[j];
	return os;
}

// operator>>
template <uint Dim>
istream& operator>>(istream& is,FCoeff<Dim>& p) {
	for (uint j=0; j<2*Dim; j++)
		is >> p[j];
	return is;
}

// writeBinary
template <uint Dim>
ostream& FCoeff<Dim>::writeBinary(ostream& os) const {
	for (uint j=0; j<2*Dim; j++)
		os.write(reinterpret_cast<const char*>(&(Coeffs[j])),sizeof(number));
	return os;
}

// readBinary
template <uint Dim>
istream& FCoeff<Dim>::readBinary(istream& is) {
	for (uint j=0; j<2*Dim; j++)
		is.read(reinterpret_cast<char*>(&(Coeffs[j])),sizeof(number));
	return is;
}	

// operator==
template <uint Dim>
bool operator==(const FCoeff<Dim>& lhs, const FCoeff<Dim>& rhs) {
	bool equal = true;
	for (uint j=0; j<2*Dim; j++) {
		if (abs(lhs[j]-rhs[j])>MIN_NUMBER) {
			equal = false;
			break;
		}
	}
	return equal;
}

// operator^=
template <uint Dim>
bool operator^=(const FCoeff<Dim>& lhs, const FCoeff<Dim>& rhs) {
	number closeness = MIN_NUMBER*1.0e4;
	bool equal = true;
	for (uint j=0; j<2*Dim; j++) {
		if (abs(lhs[j]-rhs[j])>closeness) {
			equal = false;
			break;
		}
	}
	return equal;
}


/*----------------------------------------------------------------------------------------------------------------------------
	2 - functions acting on FCoeffs

----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	3 - FLoop class
		- initialisers and destructor
		- steps
		- save
		- <<
----------------------------------------------------------------------------------------------------------------------------*/

// minimal initialiser
template <uint Dim>
FLoop<Dim>::FLoop(const uint& k, const uint& seed): 
		K(k), Seed(seed), Size(pow(2,k)), Grown(false) {
	FCoeffs.resize(Size);
	(FCoeffs[0]).zero();
	Generator = gsl_rng_alloc(gsl_rng_taus); // could also use gsl_rng_mt19937 (mersener twist)
}

// destructor
template <uint Dim>
FLoop<Dim>::~FLoop() {
	//gsl_rng_free(Generator);
	delete Generator;
}

// size
template <uint Dim>
uint FLoop<Dim>::size() const {
	return Size;
}	

// clear
template <uint Dim>
void FLoop<Dim>::clear() {
	for (uint j=0; j<Size; j++)
		FCoeffs[j].zero();
}

// set seed
template <uint Dim>
void FLoop<Dim>::setSeed(const uint& s) {
	Seed = s;
}

// indexing
template <uint Dim>
const FCoeff<Dim>& FLoop<Dim>::operator[](const uint& loc) const {
	return FCoeffs.at(loc);
}

// indexing
template <uint Dim>
FCoeff<Dim>& FLoop<Dim>::operator[](const uint& loc) {
	return FCoeffs.at(loc);
}

// save
template <uint Dim>
void FLoop<Dim>::save(const string& f) const {
	ofstream os;
	os.open(f.c_str(),ios::binary);
	if (os.good()) {
		for (uint j=0; j<Size; j++) {
			(FCoeffs[j]).writeBinary(os);
		}
		os.close();
	}
	else {
		cerr << "FLoop::save error: cannot write to " << f << endl;
		os.close();
		return;
	}
}

// saveAscii
template <uint Dim>
void FLoop<Dim>::saveAscii(const string& f) const {
	ofstream os;
	os.open(f.c_str());
	if (os.good()) {
		os << *this;
		os.close();
	}
	else {
		cerr << "FLoop::save error: cannot write to " << f << endl;
		os.close();
		return;
	}
}

// load
template <uint Dim>
void FLoop<Dim>::load(const string& f) {
	if (!fileExists(f)) {
		cerr << "FLoop::load error: file, " << f << ", doesn't exist" << endl;
		return;
	}
	uint fd = countDoubles(f);
	if (fd!=2*Dim*Size) {
		cerr << "FLoop::load error: " << f << " contains " << fd << " doubles, floop requires " << 2*Size*Dim << endl;
		return;
	}
	else {	
		ifstream is;
		is.open(f.c_str(),ios::binary);
		if (is.good()) {
			for (uint j=0; j<Size; j++) {
				(FCoeffs[j]).readBinary(is);
			}
			is.close();
		}
		else {
			cerr << "FLoop::load error: cannot read from " << f << endl;
			is.close();
			return;
		}
	}
	Grown = true;
}

// loadAscii
template <uint Dim>
void FLoop<Dim>::loadAscii(const string& f) {
	uint fl = countLines(f);
	uint fc = countColumns(f);
	if (fl!=Size) {
		cerr << "FLoop::load error: " << f << " contains " << fl << " lines, floop requires " << Size << endl;
		return;
	}
	else if (fc!=2*Dim) {
		cerr << "FLoop::load error: " << f << " contains " << fc << " columns, floop requires " << 2*Dim << endl;
		return;
	}
	else {	
		ifstream is;
		is.open(f.c_str());
		if (is.good()) {
			is >> *this;
			is.close();
		}
		else {
			cerr << "FLoop::load error: cannot read from " << f << endl;
			is.close();
			return;
		}
	}
	Grown = true;
}

// grow
template <uint Dim>
void FLoop<Dim>::grow() {
	Generator = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(Generator,Seed);
	double sigma = 1.0/PI;
		
	for (uint i=0; i<Size; i++) {
		for (uint j=0; j<Dim; j++) {
	  		(FCoeffs[i])[2*j] = gsl_ran_gaussian(Generator, sigma/(number)(i+1.0)); //gsl_ran_gaussian_ziggurat is another option
	  		(FCoeffs[i])[2*j+1] = gsl_ran_gaussian(Generator, sigma/(number)(i+1.0));
	  	}
	}

	//gsl_rng_free(Generator);
	Grown = true;
}

// getting spatial coord
template <uint Dim>
const Point<Dim> FLoop<Dim>::X(const number& t) const {
	Point<Dim> x;
	double wt = 2.0*PI*t, wti;
	for (uint i=0; i<Size; i++) {
		wti = wt*(number)(i+1.0);
		for (uint j=0; j<Dim; j++) {
			x[j] += (FCoeffs[i])[2*j]*cos(wti) + (FCoeffs[i])[2*j+1]*sin(wti);
		}
	}
	return x;
}

// getting spatial coord
template <uint Dim>
Point<Dim> FLoop<Dim>::X(const number& t) {
	Point<Dim> x;
	double wt = 2.0*PI*t, wti;
	for (uint i=0; i<Size; i++) {
		wti = wt*(number)(i+1.0);
		for (uint j=0; j<Dim; j++) {
			x[j] += (FCoeffs[i])[2*j]*cos(wti) + (FCoeffs[i])[2*j+1]*sin(wti);
		}
	}
	return x;
}

// getting d(spatial coord)
template <uint Dim>
const Point<Dim> FLoop<Dim>::dX(const number& t) const {
	Point<Dim> dx;
	double c = 2.0*PI;
	double wt = c*t, wti;
	for (uint i=0; i<Size; i++) {
		wti = wt*(number)(i+1.0);
		for (uint j=0; j<Dim; j++) {
			dx[j] += c*(i+1.0)*(-(FCoeffs[i])[2*j]*cos(wti) + (FCoeffs[i])[2*j+1]*sin(wti));
		}
	}
	return dx;
}

// getting d(spatial coord)
template <uint Dim>
Point<Dim> FLoop<Dim>::dX(const number& t) {
	Point<Dim> dx;
	double c = 2.0*PI;
	double wt = c*t, wti;
	for (uint i=0; i<Size; i++) {
		wti = wt*(number)(i+1.0);
		for (uint j=0; j<Dim; j++) {
			dx[j] += c*(i+1.0)*(-(FCoeffs[i])[2*j]*sin(wti) + (FCoeffs[i])[2*j+1]*cos(wti));
		}
	}
	return dx;
}

// getting dd(spatial coord)
template <uint Dim>
const Point<Dim> FLoop<Dim>::ddX(const number& t) const {
	Point<Dim> ddx;
	double c = 2.0*PI, ci;
	double wt = c*t, wti;
	c *= c;
	for (uint i=0; i<Size; i++) {
		wti = wt*(number)(i+1.0);
		ci = c*(i+1.0)*(i+1.0);
		for (uint j=0; j<Dim; j++) {
			ddx[j] += -ci*((FCoeffs[i])[2*j]*cos(wti) + (FCoeffs[i])[2*j+1]*sin(wti));
		}
	}
	return ddx;
}

// getting d(spatial coord)
template <uint Dim>
Point<Dim> FLoop<Dim>::ddX(const number& t) {
	Point<Dim> ddx;
	double c = 2.0*PI, ci;
	double wt = c*t, wti;
	c *= c;
	for (uint i=0; i<Size; i++) {
		wti = wt*(number)(i+1.0);
		ci = c*(i+1.0)*(i+1.0);
		for (uint j=0; j<Dim; j++) {
			ddx[j] += -ci*((FCoeffs[i])[2*j]*cos(wti) + (FCoeffs[i])[2*j+1]*sin(wti));
		}
	}
	return ddx;
}


// stream <<
template <uint Dim>
ostream& operator<< (ostream& os,const FLoop<Dim>& l) {
	for (uint j=0; j<l.Size; j++)
		os << l.FCoeffs[j] << endl;
	return os;
}

// stream >>
template <uint Dim>
istream& operator>> (istream& is, FLoop<Dim>& l) {
	for (uint j=0; j<l.Size; j++)
		is >> l.FCoeffs[j];
	l.Grown = true;
	return is;
}

/*----------------------------------------------------------------------------------------------------------------------------
	4 - FLoop functions
		
----------------------------------------------------------------------------------------------------------------------------*/

// argL
template <uint Dim>
static number argL(number t, void* void_fl) {
	FLoop<Dim>* fl = (FLoop<Dim> *) void_fl;
	return sqrt(Dot(fl->dX(t),fl->dX(t)));
}

// L
template <uint Dim> //agrees with evaluation on mathematica
number L (FLoop<Dim>& fl, \
			const number& tol, const uint& wsize, gsl_integration_workspace* ws, number& error) {
	gsl_function F;
	F.function = &argL<Dim>;
	F.params = &fl;
	number result;
	gsl_integration_qags(&F, 0.0, 1.0, 0.0, tol, wsize, ws, &result, &error); 
	return result;
}
			
// L
template <uint Dim> //agrees with evaluation on mathematica
number L (FLoop<Dim>& l, number& error) {
	number tol = 1.0e-7;
	uint wsize = 1e3;
	gsl_integration_workspace* ws = gsl_integration_workspace_alloc (wsize);;
	number result = L(l,tol,wsize,ws,error);
	gsl_integration_workspace_free(ws);
	return result;
}

// S0
template <uint Dim> //agrees with evaluation on mathematica
number S0 (const FLoop<Dim>& fl) {
	number result = 0.0;
	for (uint i=0; i<fl.size(); i++) {
		for (uint j=0; j<2*Dim; j++) {
			result += (i+1.0)*(i+1.0)*(fl[i])[j]*(fl[i])[j];
		}
	}
	return result*PI*PI/2.0;
}

// I0
template <uint Dim>
number I0 (const FLoop<Dim>& fl) {
	return 0.0;
}

// paramsV1r
template <uint Dim>
struct paramsV1r {
	FLoop<Dim>* fl;
	number a;
};

// argV1r
template <uint Dim>
static number argV1r(number* t, size_t dim, void* params) {
	(void)(dim); /* avoid unused parameter warnings */
	struct paramsV1r<Dim> *p = (struct paramsV1r<Dim> *) params;
	number a = p->a;
	FLoop<Dim>* fl = p->fl;
	number denom = a*a, result;
	denom += DistanceSquared(fl->X(t[0]),fl->X(t[1]));
	result = Dot(fl->dX(t[0]),fl->dX(t[1]));
	return result/denom;
}

// V1r
template <uint Dim>
number V1r (FLoop<Dim>& fl, const number& a, const number& tol, const uint& calls, number& error) {
	number result;

	number tl[2] = {0.0,0.0}, tu[2] = {1.0,1.0};
	gsl_rng *Generator = gsl_rng_alloc(gsl_rng_taus); // could also use gsl_rng_default

	paramsV1r<Dim> params;
	params.a = a;
	params.fl = &fl;	
	gsl_monte_function gmcf = { &argV1r<Dim>, 2, &params };

	gsl_rng_env_setup ();
	
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);

	gsl_monte_vegas_integrate (&gmcf, tl, tu, 2, calls/10, Generator, s, &result, &error);

	uint counter = 0;
	do {
		gsl_monte_vegas_integrate (&gmcf, tl, tu, 2, calls/10, Generator, s, &result, &error);
		counter++;
	}
	while (abs(gsl_monte_vegas_chisq(s) - 1.0) > 0.1 && counter<9);

	gsl_monte_vegas_free (s);

	gsl_rng_free(Generator);
	
	return result;			 
}
			 
// V1r
template <uint Dim>
number V1r (FLoop<Dim>& fl, const number& a, number& error) {
	number tol = 1.0e-2;
	uint calls = 1e5;
	return V1r(fl,a,tol,calls,error);
}

/*----------------------------------------------------------------------------------------------------------------------------
	5 - FMetropolis
----------------------------------------------------------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------------------------------------------------------
	6 - explicit template instantiation
----------------------------------------------------------------------------------------------------------------------------*/
// Dim=2
template class FCoeff<2>;
template ostream& operator<< <2>(ostream& os,const FCoeff<2>& p);
template FCoeff<2> operator+ <2>(const FCoeff<2>& lhs,const FCoeff<2>& rhs);
template FCoeff<2> operator- <2>(const FCoeff<2>& lhs,const FCoeff<2>& rhs);
template FCoeff<2> operator* <2>(const number& lhs,const FCoeff<2>& rhs);
template FCoeff<2> operator/ <2>(const FCoeff<2>& lhs,const number& rhs);
template bool operator== <2>(const FCoeff<2>& lhs, const FCoeff<2>& rhs);
template bool operator^= <2>(const FCoeff<2>& lhs, const FCoeff<2>& rhs);
template class FLoop<2>;
template ostream& operator<< <2>(ostream& os,const FLoop<2>& l);
template number S0<2> (const FLoop<2>& l);
template number argL<2>(number t, void* void_fl);
template number L<2>(FLoop<2>& l, number& error);
template number L<2>(FLoop<2>& l, const number& tol, const uint& wsize, gsl_integration_workspace* ws, number& error);
template struct paramsV1r<2>;
template number argV1r<2>(number* t, size_t dim, void* void_fl);
template number V1r<2>(FLoop<2>& fl, const number& a, const number& tol, const uint& calls, number& error);
template number V1r<2>(FLoop<2>& fl, const number& a, number& error);

// Dim=4
template class FCoeff<4>;
template ostream& operator<< <4>(ostream& os,const FCoeff<4>& p);
template FCoeff<4> operator+ <4>(const FCoeff<4>& lhs,const FCoeff<4>& rhs);
template FCoeff<4> operator- <4>(const FCoeff<4>& lhs,const FCoeff<4>& rhs);
template FCoeff<4> operator* <4>(const number& lhs,const FCoeff<4>& rhs);
template FCoeff<4> operator/ <4>(const FCoeff<4>& lhs,const number& rhs);
template bool operator== <4>(const FCoeff<4>& lhs, const FCoeff<4>& rhs);
template bool operator^= <4>(const FCoeff<4>& lhs, const FCoeff<4>& rhs);
template class FLoop<4>;
template ostream& operator<< <4>(ostream& os,const FLoop<4>& l);
template number S0<4> (const FLoop<4>& l);
template number argL<4>(number t, void* void_fl);
template number L<4>(FLoop<4>& l, number& error);
template number L<4>(FLoop<4>& l, const number& tol, const uint& wsize, gsl_integration_workspace* ws, number& error);
template struct paramsV1r<4>;
template number argV1r<4>(number* t, size_t dim, void* void_fl);
template number V1r<4>(FLoop<4>& fl, const number& a, const number& tol, const uint& calls, number& error);
template number V1r<4>(FLoop<4>& fl, const number& a, number& error);

// I0, for Dim=4
template <> number I0 <4>(const FLoop<4>& fl) {
	number result = 0.0;
	for (uint i=0; i<fl.size(); i++) {
		result += (i+1.0)*((fl[i])[2*3+1]*(fl[i])[2*2]-(fl[i])[2*2+1]*(fl[i])[2*3]);
	}
	return result*PI;
}
