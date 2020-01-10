#pragma once
#ifndef _CAMERA_CALIBRATION_H_
#define _CAMERA_CALIBRATION_H_

/*
We use the camera Canon EOS700D, the lens "EF-S10-22mm f/3.5-4.5 USM".
The degrees of (horizontal) visual angle is 97.1666667.

coo:	The first 2 elements in each columns of coo are the index of calibration pts in sketch,
		and the last 2 are the pixels of calibration pts in image.
		coo:4xn
W:		the real coordinate of calibration points. W:3xn
V:		the measured calibration points (from pixels) in the image with f = 1. V:2xn
P:		column vector with length 2n.

inv(R) = [P(1:3) P(5:7) P(9:11)]'
p0 = [P(4);P(8);P(12)]

tVi = inv(R)*Wi-inv(R)*p0
*/

#include <stdio.h>
#include <stdlib.h>

#include "mtx.h"

struct clb{
	int num_of_pts;
	mtx coo;
	mtx V;
	mtx W;

	mtx p0;
	mtx solp;

	mtx R;
	mtx orien;
};

void clb_new(clb *C, int n);
void clb_delete(clb *C);
void clb_initialize(clb *C, double *coo, double *p0);

void clb_set_cpt(clb*C, int sx, int sy, int px, int py, int count);
int clb_first_cpt_nonzero(clb *C);

void clb_evaf(clb *C, mtx* ans, mtx* p);
void clb_evaDf(clb *C, mtx* J, mtx* p);
void clb_LMlsqrsolve(clb *C);

void clb_evaOrien(clb *C);

void clb_evaCameraCoord(mtx *V, mtx *coo);
void clb_evaDomeCoord(mtx *W, mtx *coo);
void clb_evaVW(clb *C, double* coo);
void clb_evaVW(clb *C);


void clb_p2ori(mtx *p, mtx *ori, double dome_radius);
void clb_ori2p(mtx *ori, mtx *p, double dome_radius);



#endif //_CAMERA_CALIBRATION_H_