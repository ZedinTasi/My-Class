
#include "fft.h"

#define PI 3.14159265359

#define FOR(i,n0,n) for(int i=n0;i<n;i++)
#define MAX(a,b) (a>b)?a:b

complex * DFT_naive(complex * x, int N)
{
	complex *X = (complex*)malloc(sizeof(complex)*N);
	
	FOR(k, 0, N) {
		X[k].re = 0.0;
		X[k].im = 0.0;
		FOR(n, 0, N) {
			X[k] += x[n] * polar2complex(1.0, -2 * k*n* PI / N);
		}
	}
	return X;
}

/** Implements the Cooley-Tukey FFT algorithm.
*
* @expects: N1*N2 = N
*/
complex* FFT(complex* input, int N, int N1, int N2) {
	int k1, k2;

	/* Allocate columnwise matrix */
	complex** columns = (complex**)malloc(sizeof(complex*) * N1);
	for (k1 = 0; k1 < N1; k1++) {
		columns[k1] = (complex*)malloc(sizeof(complex) * N2);
	}

	/* Allocate rowwise matrix */
	complex** rows = (complex**)malloc(sizeof(complex*) * N2);
	for (k2 = 0; k2 < N2; k2++) {
		rows[k2] = (complex*)malloc(sizeof(complex) * N1);
	}

	/* Reshape input into N1 columns */
	for (k1 = 0; k1 < N1; k1++) {
		for (k2 = 0; k2 < N2; k2++) {
			columns[k1][k2] = input[N1*k2 + k1];
		}
	}

	/* Compute N1 DFTs of length N2 using naive method */
	for (k1 = 0; k1 < N1; k1++) {
		columns[k1] = DFT_naive(columns[k1], N2);
	}

	/* Multiply by the twiddle factors  ( e^(-2*pi*j/N * k1*k2)) and transpose */
	for (k1 = 0; k1 < N1; k1++) {
		for (k2 = 0; k2 < N2; k2++) {
			rows[k2][k1] = polar2complex(1, -2.0*PI*k1*k2 / N)*columns[k1][k2];
		}
	}

	/* Compute N2 DFTs of length N1 using naive method */
	for (k2 = 0; k2 < N2; k2++) {
		rows[k2] = DFT_naive(rows[k2], N1);
	}

	/* Flatten into single output */
	complex* output = (complex*)malloc(sizeof(complex) * N);
	for (k1 = 0; k1 < N1; k1++) {
		for (k2 = 0; k2 < N2; k2++) {
			output[N2*k1 + k2] = rows[k2][k1];
		}
	}

	/* Free all alocated memory except output and input arrays */
	for (k1 = 0; k1 < N1; k1++) {
		free(columns[k1]);
	}
	for (k2 = 0; k2 < N2; k2++) {
		free(rows[k2]);
	}
	free(columns);
	free(rows);
	return output;
}

complex * IFFT(complex * input, int N, int N1, int N2)
{
	complex* output = (complex*)malloc(sizeof(complex) * N);
	complex* tmp = (complex*)malloc(sizeof(complex) * N);

	output = FFT(input, N, N1, N2);
	FOR(i, 0, N) {
		output[i] /= N;
		tmp[i] = output[i];
	}
	FOR(i, 1, N) {
		output[i] = tmp[N-i];
	}
	free(tmp);

	return output;
}

void DFT_2D(mtxc * m, mtxc * M)
{
	if (m->col != M->col || m->row != M->row)return;

	int w = m->col, h = m->row;
	int k = MAX(w, h);
	complex *x = (complex*)malloc(k * sizeof(complex));
	complex *X = (complex*)malloc(k * sizeof(complex));
	mtxc TMP;
	mtxc_new(&TMP, h, w);
	FOR(i, 0, h) {
		FOR(j, 0, w) {
			x[j] = mtxc_get_value(m, i, j);
		}
		X = FFT(x, w, 1, w);
		FOR(j, 0, w) {
			mtxc_set_value(&TMP, i, j,X[j]);
		}
	}

	FOR(j, 0, w) {
		FOR(i, 0, h) {
			x[i] = mtxc_get_value(&TMP, i, j);
		}
		X = FFT(x, h, 1, h);
		FOR(i, 0, h) {
			mtxc_set_value(M, i, j, X[i]);
		}
	}

	free(x);
	free(X);
	mtxc_delete(&TMP);
}

void IDFT_2D(mtxc * M, mtxc * m)
{
	if (m->col != M->col || m->row != M->row)return;

	int w = M->col, h = M->row;
	int k = MAX(w, h);
	complex *X = (complex*)malloc(k * sizeof(complex));
	complex *x = (complex*)malloc(k * sizeof(complex));
	mtxc TMP;
	mtxc_new(&TMP, h, w);

	FOR(j, 0, w) {
		FOR(i, 0, h) {
			X[i] = mtxc_get_value(M, i, j);
		}
		x = IFFT(X, h, 1, h);
		FOR(i, 0, h) {
			mtxc_set_value(&TMP, i, j, x[i]);
		}
	}

	FOR(i, 0, h) {
		FOR(j, 0, w) {
			X[j] = mtxc_get_value(&TMP, i, j);
		}
		x = IFFT(X, w, 1, w);
		mtxc_set_value(m, i, 0, x[0]);

		FOR(j, 1, w) {
			mtxc_set_value(m, i, j, x[w-j]);
		}
	}

	free(x);
	free(X);
	mtxc_delete(&TMP);
}
