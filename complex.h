#ifndef _COMPLEX_H_
#define _COMPLEX_H_

#include <stdio.h>
#include <string>
#include <math.h>
#include "mtx.h"

struct complex 
{
	double re, im;

	complex() :re(0.0), im(0.0) {}
	complex(double r) :re(r), im(0.0) {}
	complex(double r, double i) :re(r), im(i) {}

	inline double real() { return re; }
	inline double image() { return im; }

	inline complex operator+(complex c) {
		return complex(re + c.re, im + c.im);
	}
	inline complex operator-(complex c) {
		return complex(re - c.re, im - c.im);
	}
	inline complex operator*(complex c) {
		return complex(re * c.re - im*c.im, re * c.im + im*c.re);
	}
	inline complex operator/(complex c) {
		double s = 1/(c.re*c.re + c.im*c.im);
		return complex(s*(re * c.re + im*c.im), s*(re * c.im - im*c.re));
	}

	inline void operator+=(complex c) {
		re += c.re, im += c.im;
	}
	inline void operator*=(complex c) {
		complex c1 = complex(re, im);
		c1 = c1*c;
		*this = complex(c1.re, c1.im);
	}	
	inline void operator/=(complex c) {
		complex c1 = complex(re, im);
		c1 = c1/c;
		*this = complex(c1.re, c1.im);
	}
	inline complex conj() {
		return complex(re, -im);
	}


	inline void print() {
		if (fabs(re) < 1E-8) *this = complex(0.0, im);
		if (fabs(im) < 1E-8) *this = complex(re, 0.0);
		if (im >= 0) {
			printf("%10f + %10f i", re, im);
		}
		else {
			printf("%10f - %10f i", re, -im);
		}
	}


};

complex polar2complex(double radius, double rad);


struct mtxc {
	int row, col;
	complex *v;
};

void mtxc_set_zero(mtxc *m);
void mtxc_new(mtxc *m, int row, int col);
void mtxc_delete(mtxc *m);

inline complex mtxc_get_value(mtxc *m, int i, int j);
void mtxc_set_value(mtxc *m, int i, int j, complex v);
complex mtxc_elements_sum(mtxc *m);

void mtx2mtxc(mtx *m1, mtxc *m2);
void mtxc2mtx(mtxc *m1, mtx *m2);

void mtxc_print(mtxc *m);



#endif //_COMPLEX_H_
