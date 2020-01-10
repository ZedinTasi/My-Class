#include "kalman_filter.h"

void LTI_new(LTI * syst, int n, int m, int r)
{
	syst->n = n;
	syst->m = m;
	syst->r = r;
	
	mtx_new(&syst->A, n, n);
	mtx_new(&syst->B, n, r);
	mtx_new(&syst->C, m, n);
	mtx_new(&syst->D, m, r);
	mtx_new(&syst->x0, n, 1);
}

void LTI_delete(LTI * syst)
{
	mtx_delete(&syst->A);
	mtx_delete(&syst->B);
	mtx_delete(&syst->C);
	mtx_delete(&syst->D);
	mtx_delete(&syst->x0);
}

void LTI_kalm_new(LTI_kalm * syst, int n, int m, int r)
{
	LTI_new(&syst->lti, n, m, r);
	mtx_new(&syst->x_, n, 1);
	mtx_new(&syst->P_, n, n);
	mtx_new(&syst->K, n, m);
	mtx_new(&syst->i, m, 1);
	mtx_new(&syst->s, m, m);

	mtx_new(&syst->P, n, n);
	mtx_new(&syst->Q, r, r);
	mtx_new(&syst->R, m, m);
}

void LTI_kalm_delete(LTI_kalm *syst)
{
	mtx_delete(&syst->x_);
	mtx_delete(&syst->P_);
	mtx_delete(&syst->K);
	mtx_delete(&syst->i);
	mtx_delete(&syst->s);
	mtx_delete(&syst->P);
	mtx_delete(&syst->Q);
	mtx_delete(&syst->R);

	LTI_delete(&syst->lti);
}

int kalm_predict(LTI_kalm * syst, mtx * u0)
{
	int n, m, r;
	n = syst->lti.n;
	m = syst->lti.m;
	r = syst->lti.r;

	mtx At, Bt;
	mtx tmpnn, tmpnr, tmpn1;

	mtx_new(&At, n, n);
	mtx_new(&Bt, r, n);	
	mtx_transpose(&At, &syst->lti.A);
	mtx_transpose(&Bt, &syst->lti.B);

	mtx_new(&tmpnn, n, n);
	mtx_new(&tmpnr, n, r);
	mtx_new(&tmpn1, n, 1);


	// x_ = A*x0 + B*u0
	mtx_multiply(&syst->x_, &syst->lti.A, &syst->lti.x0);
	mtx_multiply(&tmpn1, &syst->lti.B, u0);
	mtx_add(&syst->x_, &syst->x_, &tmpn1);

	// P_ = A*P0*A' + B*Q*B';
	mtx_multiply(&tmpnn, &syst->lti.A, &syst->P);
	mtx_multiply(&syst->P_, &tmpnn, &At);
	mtx_multiply(&tmpnr, &syst->lti.B, &syst->Q);
	mtx_multiply(&tmpnn, &tmpnr, &Bt);
	mtx_add(&syst->P_, &syst->P_, &tmpnn);



	mtx_delete(&At);
	mtx_delete(&Bt);
	mtx_delete(&tmpnn);
	mtx_delete(&tmpnr);
	mtx_delete(&tmpn1);

	return 0;
}

int kalm_update(LTI_kalm * syst, mtx *y0)
{
	int n, m, r;
	n = syst->lti.n;
	m = syst->lti.m;
	r = syst->lti.r;
	mtx Ct, Inn;
	mtx tmpnn, tmpnm, tmpn1, tmpmn, tmpmm, tmpm1;

	mtx_new(&Ct, n, m);
	mtx_new(&Inn, n, n);
	mtx_transpose(&Ct, &syst->lti.C);
	mtx_set_identity(&Inn);

	mtx_new(&tmpnn, n, n);
	mtx_new(&tmpnm, n, m);
	mtx_new(&tmpn1, n, 1);
	mtx_new(&tmpmn, m, n);
	mtx_new(&tmpmm, m, m);
	mtx_new(&tmpm1, m, 1);


	// i = y0 - C*x1_
	mtx_multiply(&tmpm1, &syst->lti.C, &syst->x_);
	mtx_subtract(&syst->i, y0, &tmpm1);

	// s = C*P_*C' + R (the covariance of i)
	mtx_multiply(&tmpmn, &syst->lti.C, &syst->P_);
	mtx_multiply(&tmpmm, &tmpmn, &Ct);
	mtx_add(&syst->s, &tmpmm, &syst->R);

	// K = P_*C'*s^(-1)
	mtx_solve_inverse(&syst->s, &syst->s);
	mtx_multiply(&tmpnm, &syst->P_, &Ct);
	mtx_multiply(&syst->K, &tmpnm, &syst->s);

	// x0 = x_ + K*i
	mtx_multiply(&tmpn1, &syst->K, &syst->i);
	mtx_add(&syst->lti.x0, &syst->x_, &tmpn1);

	// P1 = (I-K*C)*P_;
	mtx_multiply(&tmpnn, &syst->K, &syst->lti.C);
	mtx_subtract(&tmpnn, &Inn, &tmpnn);
	mtx_multiply(&syst->P, &tmpnn, &syst->P_);


	mtx_delete(&Ct);
	mtx_delete(&Inn);

	mtx_delete(&tmpnn);
	mtx_delete(&tmpnm);
	mtx_delete(&tmpn1);
	mtx_delete(&tmpmn);
	mtx_delete(&tmpmm);
	mtx_delete(&tmpm1);
	return 0;
}


//int kalm_update(mtx * x1, mtx * P1, LTI *syst, mtx *y1, mtx * x0, mtx * P0, mtx * Q, mtx * R)
//{
//	int n, m, r;
//	n = syst->n;
//	m = syst->m;
//	r = syst->r;
//
//	if (x1->row != n || x1->col != 1 || x0->row != n || x0->col != 1 ) return -1;
//	if (P1->row != n || P1->col != n || P0->row != n || P0->col != n ) return -1;
//	if (y1->row != m || y1->col != 1) return -1;
//	if (Q->row != r || Q->col != r || R->row != m || R->col != m ) return -1;
//	
//	mtx x1_, P1_, K, i, s;
//	mtx	At, Bt, Ct, Inn;
//	mtx tmpnn, tmpnm, tmpnr, tmpn1, tmpmn, tmpmm, tmpm1;
//	
//	mtx_new(&x1_, n, 1);
//	mtx_new(&P1_, n, n);
//	mtx_new(&K, n, m);
//	mtx_new(&i, m, 1);
//	mtx_new(&s, m, m);
//
//	mtx_new(&At, n, n);
//	mtx_new(&Bt, r, n);
//	mtx_new(&Ct, n, m);
//	mtx_new(&Inn, n, n);
//	mtx_transpose(&At, &syst->A);
//	mtx_transpose(&Bt, &syst->B);
//	mtx_transpose(&Ct, &syst->C);
//	mtx_set_identity(&Inn);
//
//	mtx_new(&tmpnn, n, n);
//	mtx_new(&tmpnm, n, m);
//	mtx_new(&tmpnr, n, r);
//	mtx_new(&tmpn1, n, 1);
//	mtx_new(&tmpmn, m, n);
//	mtx_new(&tmpmm, m, m);
//	mtx_new(&tmpm1, m, 1);
//
//
//	// x1_ = A*x0 (+ Bu0)
//	mtx_multiply(&x1_, &syst->A, x0); 
//
//	// P1_ = A*P0*A' + B*Q*B';
//	mtx_multiply(&tmpnn, &syst->A, P0);
//	mtx_multiply(P0, &tmpnn, &At);
//	mtx_multiply(&tmpnr, &syst->B, Q);
//	mtx_multiply(&tmpnn, &tmpnr, &Bt);
//	mtx_add(&P1_, P0, &tmpnn);
//
//	// i = y1 - C*x1_
//	mtx_multiply(&tmpm1, &syst->C, &x1_);
//	mtx_subtract(&i, y1, &tmpm1);
//
//	// s = C*P1_*C' + R (the covariance of i)
//	mtx_multiply(&tmpmn, &syst->C, &P1_);
//	mtx_multiply(&tmpmm, &tmpmn, &Ct);
//	mtx_add(&s, &tmpmm, R);
//	
//	// K = P1_*C'*s^(-1)
//	mtx_solve_inverse(&s,&s);
//	mtx_multiply(&tmpnm, &P1_, &Ct);
//	mtx_multiply(&K, &tmpnm, &s);
//
//	// x1 = x1_ + K*i
//	mtx_multiply(&tmpn1, &K, &i);
//	mtx_add(x1, &x1_, &tmpn1);
//
//	// P1 = (I-K*C)*P1_;
//	mtx_multiply(&tmpnn, &K, &syst->C);
//	mtx_subtract(&tmpnn, &Inn, &tmpnn);
//	mtx_multiply(P1, &tmpnn, &P1_);
//
//
//	mtx_delete(&x1_);
//	mtx_delete(&P1_);
//	mtx_delete(&K);
//	mtx_delete(&i);
//	mtx_delete(&s);
//
//	mtx_delete(&At);
//	mtx_delete(&Bt);
//	mtx_delete(&Ct);
//	mtx_delete(&Inn);
//
//	mtx_delete(&tmpnn);
//	mtx_delete(&tmpnm);
//	mtx_delete(&tmpnr);
//	mtx_delete(&tmpn1);
//	mtx_delete(&tmpmn);
//	mtx_delete(&tmpmm);
//	mtx_delete(&tmpm1);
//
//}