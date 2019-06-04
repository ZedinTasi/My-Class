#include "gaussian_blur.h"


void gaussfilter_new(gaussfilter * G, int size)
{
	if (size < 1) return;
	G->se = 1;
	G->size = size;
	mtx_new(&G->kernel, size, size);
}

void gaussfilter_delete(gaussfilter *G)
{
	mtx_delete(&G->kernel);
}

void gaussfilter_compute(gaussfilter * G)
{
	if (G->se < 0 || !&G->kernel) return;
	
	double v = 2*G->se*G->se;
	int ori = (G->size - 1) / 2;
	double x, y;
	double Gij;

	for (int i = 0; i < G->size; i++) {
		for (int j = 0; j < G->size; j++) {
			y = i - ori;
			x = j - ori;
			Gij = exp((-x*x-y*y)/v)/v/PI;
			mtx_set_value(&G->kernel, i, j, Gij);
		}
	}
	mtx_normalize(&G->kernel);

}

void gaussfilter_compute(gaussfilter * G, double sigma)
{
	G->se = sigma;
	gaussfilter_compute(G);
}

void gaussfilter_inv(gaussfilter * G1, gaussfilter * G2, double SNR)
{
	int N = G1->size;

	mtxc m, M;
	mtxc_new(&m, N, N);
	mtxc_new(&M, N, N);
	mtx2mtxc(&G1->kernel, &m);

	DFT_2D(&m, &M);
	if (fabs(SNR) > 1e-6) {
		complex H, H_;
		double nsr = 1.0 / SNR;	
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				H = mtxc_get_value(&M, i, j);
				H_ = H.conj();
				H = H_ / (H*H_ + nsr);
				mtxc_set_value(&M, i, j, H);
			}
		}
	}
	else {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				mtxc_set_value(&M, i, j, complex());
			}
		}
	}

	IDFT_2D(&M, &m);
	mtxc2mtx(&m, &G2->kernel);
	mtxc_delete(&m);
	mtxc_delete(&M);
}

void gaussfilter_inv(gaussfilter * G1, gaussfilter * G2)
{
	gaussfilter_inv(G1, G2, 1E16);
	//int N = G1->size;

	//mtxc m, M;
	//mtxc_new(&m, N, N);
	//mtxc_new(&M, N, N);
	//mtx2mtxc(&G1->kernel, &m);

	//DFT_2D(&m, &M);
	//for (int i = 0; i < N; i++) {
	//	for (int j = 0; j < N; j++) {
	//		mtxc_set_value(&M, i, j, complex(1.0) / mtxc_get_value(&M, i, j));
	//	}
	//}
	//IDFT_2D(&M, &m);
	//mtxc2mtx(&m, &G2->kernel);
	//mtxc_delete(&m);
	//mtxc_delete(&M);
}


double gaussfilter_get_value(gaussfilter * G, int i, int j)
{
	return mtx_get_value(&G->kernel, i, j);
}
