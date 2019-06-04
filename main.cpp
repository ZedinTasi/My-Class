#include "complex.h"
#include "fft.h"
//#include "mtx.h"
#include "gaussian_blur.h"

#include <stdio.h>
#include <string>
#include <math.h>

#define PI 3.14159265359
#define FOR(i,n0,n) for(int i=n0;i<n;i++)

int main() {


	//int N = 9;
	//mtxc m, M, m_;
	//mtxc_new(&m, N, N);
	//mtxc_new(&M, N, N);
	//mtxc_new(&m_, N, N);
	////complex weight[3] = { 0.150342,0.094907,0.023792 };
	////complex s[3] = { 1.0 ,weight[1] / weight[0] , weight[2] / weight[0] };
	//complex weight[5] = { 0.146634,0.092566,0.023205,0.002289,0.000088 };
	//complex s[5] = { 1.0 ,weight[1] / weight[0] , weight[2] / weight[0],weight[3] / weight[0],weight[4] / weight[0] };
	//
	//int i = (N-1)/2;
	//mtxc_set_value(&m, i, i, weight[0]);
	//FOR(x, 1, i + 1) {
	//	mtxc_set_value(&m, i + x, i, weight[x]);
	//	mtxc_set_value(&m, i - x, i, weight[x]);
	//	mtxc_set_value(&m, i, i + x, weight[x]);
	//	mtxc_set_value(&m, i, i - x, weight[x]);
	//}
	//FOR(x, 1, i+1) {
	//	FOR(y, 1, i+1) {
	//		mtxc_set_value(&m, i + x, i + y, weight[x] * s[y]);
	//		mtxc_set_value(&m, i - x, i + y, weight[x] * s[y]);
	//		mtxc_set_value(&m, i + x, i - y, weight[x] * s[y]);
	//		mtxc_set_value(&m, i - x, i - y, weight[x] * s[y]);
	//	}
	//}
	//DFT_2D(&m,&M);
	//
	//FOR(i,0,N) {
	//	FOR(j,0,N) {
	//		mtxc_set_value(&M, i, j, complex(1.0) / mtxc_get_value(&M, i, j));
	//	}
	//}

	//IDFT_2D(&M, &m_);

	//mtxc_print(&m);
	//mtxc_print(&M);
	//mtxc_print(&m_);

	//mtxc_elements_sum(&m_).print();

	//mtxc_delete(&m);
	//mtxc_delete(&M);
	//mtxc_delete(&m_);

	//system("pause");

	//========================================================
	int N = 7;

	gaussfilter g,G;
	gaussfilter_new(&g, N);
	gaussfilter_new(&G, N);
	gaussfilter_compute(&g, 0.5);

	mtx_print(&g.kernel);

	gaussfilter_inv(&g, &G, 10.0);
	printf("SNR = 10:\n");
	mtx_print(&G.kernel);

	gaussfilter_inv(&g, &G);
	printf("SNR = inf:\n");
	mtx_print(&G.kernel);

	gaussfilter_inv(&g, &G, 0.0);
	printf("SNR = 0:\n");
	mtx_print(&G.kernel);

	gaussfilter_delete(&g);
	gaussfilter_delete(&G);

	system("pause");

	

}