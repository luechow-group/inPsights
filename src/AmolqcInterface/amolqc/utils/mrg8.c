#include <math.h>
#include <stdbool.h>
#include <stdio.h>

unsigned long coeff[8] = {
	1089656042, 1906537547, 1764115693, 1304127872,
  189748160, 1984088114, 626062218, 1927846343
};
long long v[8];
int y;

bool iset;
double gset;

static const long q = 4294967295;
static const int m = 2147483647;
static const int a = 1812433253;
static const int s = 30;

double init_mrgran_c_(int*);
long long mrg_intran_c_();
double mrg_ran_c_();
double mrg_gran_c_();

double init_mrgran_c_(int* seed_) {
	int seed = *seed_;
	v[0] = seed;
	for(int i = 1; i < 8; i++) {
		v[i] = (a * (v[i-1] ^ (v[i-1] >> s)) + i) & q;
	}
	y = 0;
	return mrg_ran_c_();
}

long long mrg_intran_c_() {
	long long c = 0;
	for(int i = 0; i < y; i++) {
		c += coeff[i] * v[y-i-1]; // coeff from 0 to y-1, v from y-1 to 0
	}
	for(int i = y; i < 8; i++) {
		c += coeff[i] * v[7-i+y]; // coeff from y to 7, v from 7 to y
	}
	c = c % m;
	v[y] = c;
	y = (y + 1) & 7;

	return c;
}

double mrg_ran_c_() {
	static const long long m = 9007199254740992LL;
	static const int k = 1073741823;
	double r = 0;
	while(r == 0) {
		long long p = mrg_intran_c_() & k;
		long long q = mrg_intran_c_() & k;
		r = ((p << 30) | q) & (m - 1);
		r = r/m;
	}
	return r;
}

double mrg_gran_c_() {
	if(iset) {
		iset = false;
		return gset;
	}

	double rsq = 0, v1, v2;
	while(rsq == 0.0 || rsq >= 1.0) {
		v1 = 2.0 * mrg_ran_c_() - 1.0;
		v2 = 2.0 * mrg_ran_c_() - 1.0;
		rsq = v1 * v1 + v2 * v2;
	}
	double fac = sqrt(-2.0*log(rsq)/rsq);
	iset = true;
	gset = v1 * fac;
	return v2 * fac;
}
