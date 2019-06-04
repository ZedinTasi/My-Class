#pragma once
#ifndef _GAUSSIAN_BLUR_H_
#define _GAUSSIAN_BLUR_H_

#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265359
#include "mtx.h"
#include "complex.h"
#include "fft.h"

struct gaussfilter {
	double se;
	int size;
	mtx kernel;
};

void gaussfilter_new(gaussfilter *G, int size);
void gaussfilter_delete(gaussfilter *G);

void gaussfilter_compute(gaussfilter *G);
void gaussfilter_compute(gaussfilter *G, double sigma);
void gaussfilter_inv(gaussfilter *G1, gaussfilter *G2, double SNR);
void gaussfilter_inv(gaussfilter *G1, gaussfilter *G2);

double gaussfilter_get_value(gaussfilter *G, int i, int j);



#endif //_GAUSSIAN_BLUR_H_