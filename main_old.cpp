#include "mtx.h"

//#include "complex.h"
//#include "fft.h"
//#include "gaussian_blur.h"
//#include "camera_calibration.h"
//#include "jinc.h"
#include "kalman_filter.h"
#include "random.h"


#include <stdio.h>
#include <string>
#include <math.h>

#define PI 3.14159265359
#define FOR(i,n0,n) for(int i=n0;i<n;i++)

int main0() {

	// //test the complex matrix and FFT

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
	//
	//mtxc_print(&m);
	//mtxc_print(&M);
	//mtxc_print(&m_);

	//mtxc_elements_sum(&m_).print();

	//mtxc_delete(&m);
	//mtxc_delete(&M);
	//mtxc_delete(&m_);


	//========================================================

	// //test the gauss filter

	//int N = 7;

	//gaussfilter g,G;
	//gaussfilter_new(&g, N);
	//gaussfilter_new(&G, N);
	//gaussfilter_compute(&g, 0.5);

	//mtx_print(&g.kernel);

	//gaussfilter_inv(&g, &G, 10.0);
	//printf("SNR = 10:\n");
	//mtx_print(&G.kernel);

	//gaussfilter_inv(&g, &G);
	//printf("SNR = inf:\n");
	//mtx_print(&G.kernel);

	//gaussfilter_inv(&g, &G, 0.0);
	//printf("SNR = 0:\n");
	//mtx_print(&G.kernel);

	//gaussfilter_delete(&g);
	//gaussfilter_delete(&G);


	//========================================================

	// test the camera calibration  
	// use Levenberg-Marquardt algorithm to solve least square problem.

	//int N = 6;
	//clb C;
	//clb_new(&C, N);
	//double coo[24] = {
	//	-5., -5., -5.,      0.,      0.,      0.,
	//	1.,      2.,      3.,      1.,      2.,      3.,
	//	2331.,   2354.,   2343.,   3743.,   3782.,   3754.,
	//	1419.,   1974.,   2545.,   1295.,   2011.,   2742.
	//};	
	//double p0[12] = {
	//	0.7237861,
	//	0.0176256,
	//	0.6893513,
	//	-0.2285411,
	//	-0.0206042,
	//	1.0003858,
	//	-0.0058583,
	//	0.089679,
	//	-0.6889121,
	//	-0.0093001,
	//	0.7247852,
	//	0.1614007
	//};
	//FOR(i, 0, 12) {
	//	p0[i] += 0.1;
	//}

	//clb_initialize(&C, coo, p0);

	//printf("coo:");
	//mtx_print(&C.coo);
	//printf("\nW:");
	//mtx_print(&C.W);
	//printf("\nV:");
	//mtx_print(&C.V);
	//printf("\np0:");
	//mtx_print(&C.p0);

	//clb_LMlsqrsolve(&C);
	//clb_evaOrien(&C);

	//printf("\nsolp:");
	//mtx_print(&C.solp);	
	//printf("\nR:");
	//mtx_print(&C.R);
	//printf("\norien:");
	//mtx_print(&C.orien);


	//clb_delete(&C);
	//system("pause");

	//mtx ori0, p;
	//mtx_new(&ori0, 6, 1);
	//mtx_new(&p, 12, 1);
	//double a[6] = { 1.1,0.1,3.1,0.1,45.1,10.1 };
	//FOR(i, 0, 6) {
	//	mtx_set_value(&ori0, i, 0, a[i]);
	//}

	//clb_ori2p(&ori0, &p, 8.0);

	//printf("\np:");
	//mtx_print(&p);

	//clb_p2ori(&p, &ori0, 8.0);

	//printf("\nori:");
	//mtx_print(&ori0);



	//mtx_delete(&ori0);
	//mtx_delete(&p);


	//========================================================

	// //test the kalman filter

	//int n = 3, m = 1, r = 1;

	//double A[9] = {
	//	1.1269,-1.0144,0.1129,
	//	1,0,0,
	//	0,1,0
	//};
	//double B[3] = {
	//	-0.3832,
	//	0.5919,
	//	0.5191
	//};
	//double C[3] = {
	//	1,0,0
	//};

	//LIT_kalm kft; // kalman_filter_test
	//LIT_kalm_new(&kft, 3, 1, 1);

	//mtx xr, y0, u0; // real x
	//mtx_new(&xr, n, 1);
	//mtx_new(&y0, m, 1);
	//mtx_new(&u0, r, 1);


	//double x[3] = {
	//	10, 0, 0
	//};
	//double P[9] = {
	//	100, 0, 0,
	//	0,0,0,
	//	0,0,0
	//};

	//mtx_copy(&kft.lit.A, A);
	//mtx_copy(&kft.lit.B, B);
	//mtx_copy(&kft.lit.C, C);
	//mtx_copy(&xr, x);
	//mtx_copy(&kft.lit.x0, &xr);
	//mtx_multiply(&y0, &kft.lit.C, &xr);

	//mtx_copy(&kft.P, P);
	//mtx_set_value(&kft.Q, 0, 0, 0);
	//mtx_set_value(&kft.R, 0, 0, 1);



	//mtx_print_transpose(&kft.lit.x0);
	//FOR(i, 0, 300) {

	//	//predict the estimate x_ of x1 by x0 and input u0, and then x0 be x_{-1}.
	//	kalm_predict(&kft, &u0);

	//	// y0 is the measured output
	//	mtx_multiply(&xr, &kft.lit.A, &xr);
	//	mtx_multiply(&y0, &kft.lit.C, &xr);
	//	//mtx_add(&y0, &y0, (double)(i+1) / 10.0);
	//	mtx_add(&y0, &y0, random_normal());


	//	//update the estimate x0 by x_ and the measurement of output y0.
	//	kalm_update(&kft, &y0);

	//	mtx_print_transpose(&kft.lit.x0);
	//}


	//LIT_kalm_delete(&kft);
	//mtx_delete(&xr);
	//mtx_delete(&y0);
	//mtx_delete(&u0);


	//========================================================

	// //test the kalman filter

int n = 3, m = 1, r = 1;

double A[9] = {
	1.1269,-1.0144,0.1129,
	1,0,0,
	0,1,0
};
double B[3] = {
	-0.3832,
	0.5919,
	0.5191
};
double C[3] = {
	1,0,0
};

LTI_kalm kft; // kalman_filter_test
LTI_kalm_new(&kft, 3, 1, 1);

mtx xr, y0, u0; // real x
mtx_new(&xr, n, 1);
mtx_new(&y0, m, 1);
mtx_new(&u0, r, 1);


double x[3] = {
	10, 0, 0
};
double P[9] = {
	100, 0, 0,
	0,0,0,
	0,0,0
};

mtx_copy(&kft.lti.A, A);
mtx_copy(&kft.lti.B, B);
mtx_copy(&kft.lti.C, C);
mtx_copy(&xr, x);
mtx_copy(&kft.lti.x0, &xr);
mtx_multiply(&y0, &kft.lti.C, &xr);

mtx_copy(&kft.P, P);
mtx_set_value(&kft.Q, 0, 0, 0);
mtx_set_value(&kft.R, 0, 0, 1);

mtx_print_transpose(&kft.lti.x0);
FOR(i, 0, 300) {

	//predict the estimate x_ of x1 by x0 and input u0, and then x0 be x_{-1}.
	kalm_predict(&kft, &u0);

	// y0 is the measured output
	mtx_multiply(&xr, &kft.lti.A, &xr);
	mtx_multiply(&y0, &kft.lti.C, &xr);
	//mtx_add(&y0, &y0, (double)(i+1) / 10.0);
	mtx_add(&y0, &y0, random_normal());


	//update the estimate x0 by x_ and the measurement of output y0.
	kalm_update(&kft, &y0);

	mtx_print_transpose(&kft.lti.x0);
}


LTI_kalm_delete(&kft);
mtx_delete(&xr);
mtx_delete(&y0);
mtx_delete(&u0);


//========================================================

	system("pause");
	return 0;
}