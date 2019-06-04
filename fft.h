#ifndef _FFT_H_
#define _FFT_H_

#include "complex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

complex* DFT_naive(complex* x, int N);
complex* FFT(complex* input, int N, int N1, int N2);
complex* IFFT(complex* input, int N, int N1, int N2);

void DFT_2D(mtxc * m, mtxc * M);

void IDFT_2D(mtxc * M, mtxc * m);




#endif //_FFT_H_
