#pragma once
#ifndef _KALMAN_FILTER_
#define _KALMAN_FILTER_

#include "mtx.h"

struct LTI {
	//	x(k+1) = A*xk + B*uk + wk
	//		yk = C*xk + D*uk + vk
	// xk	:nx1
	// yk	:mx1
	// uk	:rx1

	int n, m, r;
	mtx A;	//nxn
	mtx B;  //nxr
	mtx C;	//mxn
	mtx D;	//mxr

	mtx x0;

};

struct LTI_kalm {
	//	x(k+1) = A*xk + B*uk + wk
	//		yk = C*xk + D*uk + vk
	// xk	:nx1
	// yk	:mx1
	// uk	:rx1

	LTI lti;
	mtx P; // updated error covariance		at t=0 based on data up to t=0 

	mtx x_; // updated estimate			at t=1 based on data up to t=0 
	mtx P_; // updated error covariance	at t=1 based on data up to t=0 

	mtx K; // Kalman gain at t=1
	mtx i; // innovance at t=1
	mtx s; // innovance covariance at t=1
	
	mtx Q; // covariance of wk
	mtx R; // covariance of vk
};

void LTI_new(LTI *syst, int n, int m, int r);
void LTI_delete(LTI * syst);

void LTI_kalm_new(LTI_kalm *syst, int n, int m, int r);
void LTI_kalm_delete(LTI_kalm *syst);

//int kalm_update(mtx *x1, mtx *P1, LTI *syst, mtx *y1, mtx *x0, mtx *P0, mtx *Q, mtx *R);

int kalm_predict(LTI_kalm *syst, mtx *u0);
int kalm_update(LTI_kalm *syst, mtx *y0);

#endif //_KALMAN_FILTER_