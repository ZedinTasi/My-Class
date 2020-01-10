#include "center_shift.h"

#define max_itr_times 100

void pts_shift(mtx * w, mtx * v, mtx * x)
{
	if (w->row != 3 || v->row != 3 || x->row != 6) return;
	if (w->col != 6 || v->col != 6 || x->col != 1) return;

	mtx Ra, Rb, Rc, R;
	mtx_new(&Ra, 3, 3);
	mtx_new(&Rb, 3, 3);
	mtx_new(&Rc, 3, 3);
	mtx_new(&R, 3, 3);
	mtx_rotation(&Ra, mtx_get_value(x, 3, 0), 0);
	mtx_rotation(&Rb, mtx_get_value(x, 4, 0), 1);
	mtx_rotation(&Rc, mtx_get_value(x, 5, 0), 2);
	mtx_multiply(&R, &Ra, &Rb);
	mtx_multiply(&R, &R, &Rc);
	mtx vi;
	mtx_new(&vi, 3, 1);

	for (int i = 0; i < 6; i++) {
		mtx_set_value(&vi, 0, 0, mtx_get_value(v, 0, i));
		mtx_set_value(&vi, 1, 0, mtx_get_value(v, 1, i));
		mtx_set_value(&vi, 2, 0, mtx_get_value(v, 2, i));
		mtx_multiply(&vi, &R, &vi);

		mtx_set_value(&vi, 0, 0, mtx_get_value(&vi, 0, 0) + mtx_get_value(x, 0, 0));
		mtx_set_value(&vi, 1, 0, mtx_get_value(&vi, 1, 0) + mtx_get_value(x, 1, 0));
		mtx_set_value(&vi, 2, 0, mtx_get_value(&vi, 2, 0) + mtx_get_value(x, 2, 0));
		for (int j = 0; j < 3; j++) {
			mtx_set_value(w, j, i, mtx_get_value(&vi, j, 0));
		}
	}
	mtx_delete(&vi);	
	mtx_delete(&Ra);
	mtx_delete(&Rb);
	mtx_delete(&Rc);
	mtx_delete(&R);
}

void eva_F(mtx * F, mtx * x, mtx * v, mtx * w)
{
	if (F->col != 1 || F->row != 6) return;
	if (x->col != 1 || x->row != 6) return;
	if (v->col != 6 || v->row != 3)return;
	if (w->col != 6 || w->row != 3)return;

	double a, b, c,tmpd;
	mtx Ra, Rb, Rc, R, tmp, wi;
	a = mtx_get_value(x, 3, 0);
	b = mtx_get_value(x, 4, 0);
	c = mtx_get_value(x, 5, 0);
	mtx_new(&Ra, 3, 3);
	mtx_new(&Rb, 3, 3);
	mtx_new(&Rc, 3, 3);
	mtx_new(&R, 3, 3);
	mtx_rotation(&Ra, a, 0);
	mtx_rotation(&Rb, b, 1);
	mtx_rotation(&Rc, c, 2);

	mtx_multiply(&R, &Ra, &Rb);
	mtx_multiply(&R, &R, &Rc);

	mtx_new(&tmp, 3, 1);
	mtx_new(&wi, 3, 1);
	for (int i = 0; i < 6; i++) 
	{
		mtx_set_value(&tmp, 0, 0, mtx_get_value(v, 0, i));
		mtx_set_value(&tmp, 1, 0, mtx_get_value(v, 1, i));
		mtx_set_value(&tmp, 2, 0, mtx_get_value(v, 2, i));
		mtx_multiply(&tmp, &R, &tmp);

		mtx_set_value(&tmp, 0, 0, mtx_get_value(&tmp, 0, 0)+ mtx_get_value(x, 0, 0));
		mtx_set_value(&tmp, 1, 0, mtx_get_value(&tmp, 1, 0) + mtx_get_value(x, 1, 0));
		mtx_set_value(&tmp, 2, 0, mtx_get_value(&tmp, 2, 0) + mtx_get_value(x, 2, 0));

		mtx_set_value(&wi, 0, 0, mtx_get_value(w, 0, i));
		mtx_set_value(&wi, 1, 0, mtx_get_value(w, 1, i));
		mtx_set_value(&wi, 2, 0, mtx_get_value(w, 2, i));

		tmpd = mtx_get_value(&wi, 0, 0)*mtx_get_value(&tmp, 0, 0) + mtx_get_value(&wi, 1, 0)*mtx_get_value(&tmp, 1, 0) + mtx_get_value(&wi, 2, 0)*mtx_get_value(&tmp, 2, 0);
		mtx_set_value(F, i, 0, tmpd*tmpd - mtx_Fnorm(&tmp)*mtx_Fnorm(&wi)*mtx_Fnorm(&tmp)*mtx_Fnorm(&wi));
	}


	mtx_delete(&tmp);
	mtx_delete(&wi);
	mtx_delete(&Ra);
	mtx_delete(&Rb);
	mtx_delete(&Rc);
	mtx_delete(&R);
}

void eva_DF(mtx * J, mtx * x, mtx * v, mtx * w)
{
	if (J->col != 6 || J->row != 6) return;
	if (x->col != 1 || x->row != 6) return;
	if (v->col != 6 || v->row != 3)return;
	if (w->col != 6 || w->row != 3)return;

	double a, b, c, tmpd1, tmpd2, tmpd3, witwi;
	mtx Ra, Rb, Rc, DRa, DRb, DRc, DRvi, R, tmp, wi, vi;
	a = mtx_get_value(x, 3, 0);
	b = mtx_get_value(x, 4, 0);
	c = mtx_get_value(x, 5, 0);
	mtx_new(&Ra, 3, 3);
	mtx_new(&Rb, 3, 3);
	mtx_new(&Rc, 3, 3);
	mtx_new(&R, 3, 3);
	mtx_new(&DRa, 3, 3);
	mtx_new(&DRb, 3, 3);
	mtx_new(&DRc, 3, 3);

	mtx_rotation(&Ra, a, 0);
	mtx_rotation(&Rb, b, 1);
	mtx_rotation(&Rc, c, 2);

	mtx_set_value(&DRa, 0, 0, -sin(a));
	mtx_set_value(&DRa, 0, 1, -cos(a));
	mtx_set_value(&DRa, 1, 0, cos(a));
	mtx_set_value(&DRa, 1, 1, -sin(a));

	mtx_set_value(&DRb, 2, 2, -sin(b));
	mtx_set_value(&DRb, 2, 0, -cos(b));
	mtx_set_value(&DRb, 0, 2, cos(b));
	mtx_set_value(&DRb, 0, 0, -sin(b));

	mtx_set_value(&DRc, 1, 1, -sin(c));
	mtx_set_value(&DRc, 1, 2, -cos(c));
	mtx_set_value(&DRc, 2, 1, cos(c));
	mtx_set_value(&DRc, 2, 2, -sin(c));

	mtx_multiply(&R, &Ra, &Rb);
	mtx_multiply(&DRa, &DRa, &Rb);
	mtx_multiply(&DRa, &DRa, &Rc);	
	mtx_multiply(&DRb, &Ra, &DRb);
	mtx_multiply(&DRb, &DRb, &Rc);
	mtx_multiply(&DRc, &R, &DRc);
	mtx_multiply(&R, &R, &Rc);

	mtx_new(&tmp, 3, 1);
	mtx_new(&wi, 3, 1);
	mtx_new(&vi, 3, 1);
	mtx_new(&DRvi, 3, 1);
	for (int i = 0; i < 6; i++) {
		mtx_set_value(&tmp, 0, 0, mtx_get_value(v, 0, i));
		mtx_set_value(&tmp, 1, 0, mtx_get_value(v, 1, i));
		mtx_set_value(&tmp, 2, 0, mtx_get_value(v, 2, i));
		mtx_multiply(&tmp, &R, &tmp);

		mtx_set_value(&tmp, 0, 0, mtx_get_value(&tmp, 0, 0) + mtx_get_value(x, 0, 0));
		mtx_set_value(&tmp, 1, 0, mtx_get_value(&tmp, 1, 0) + mtx_get_value(x, 1, 0));
		mtx_set_value(&tmp, 2, 0, mtx_get_value(&tmp, 2, 0) + mtx_get_value(x, 2, 0));

		mtx_set_value(&wi, 0, 0, mtx_get_value(w, 0, i));
		mtx_set_value(&wi, 1, 0, mtx_get_value(w, 1, i));
		mtx_set_value(&wi, 2, 0, mtx_get_value(w, 2, i));
		mtx_set_value(&vi, 0, 0, mtx_get_value(v, 0, i));
		mtx_set_value(&vi, 1, 0, mtx_get_value(v, 1, i));
		mtx_set_value(&vi, 2, 0, mtx_get_value(v, 2, i));

		tmpd1 = mtx_get_value(&wi, 0, 0)*mtx_get_value(&tmp, 0, 0) + mtx_get_value(&wi, 1, 0)*mtx_get_value(&tmp, 1, 0) + mtx_get_value(&wi, 2, 0)*mtx_get_value(&tmp, 2, 0);
		witwi = mtx_Fnorm(&wi)*mtx_Fnorm(&wi);

		mtx_set_value(J, i, 0, 2 * tmpd1*mtx_get_value(&wi, 0, 0) - 2 * mtx_get_value(&tmp, 0, 0)*witwi);
		mtx_set_value(J, i, 1, 2 * tmpd1*mtx_get_value(&wi, 1, 0) - 2 * mtx_get_value(&tmp, 1, 0)*witwi);
		mtx_set_value(J, i, 2, 2 * tmpd1*mtx_get_value(&wi, 2, 0) - 2 * mtx_get_value(&tmp, 2, 0)*witwi);
		
		mtx_multiply(&DRvi, &DRa, &vi);
		tmpd2 = mtx_get_value(&wi, 0, 0)*mtx_get_value(&DRvi, 0, 0) + mtx_get_value(&wi, 1, 0)*mtx_get_value(&DRvi, 1, 0) + mtx_get_value(&wi, 2, 0)*mtx_get_value(&DRvi, 2, 0);
		tmpd3 = mtx_get_value(&DRvi, 0, 0)*mtx_get_value(&tmp, 0, 0) + mtx_get_value(&DRvi, 1, 0)*mtx_get_value(&tmp, 1, 0) + mtx_get_value(&DRvi, 2, 0)*mtx_get_value(&tmp, 2, 0);
		mtx_set_value(J, i, 3, 2 * (tmpd1*tmpd2 - tmpd3*witwi));		
		mtx_multiply(&DRvi, &DRb, &vi);
		tmpd2 = mtx_get_value(&wi, 0, 0)*mtx_get_value(&DRvi, 0, 0) + mtx_get_value(&wi, 1, 0)*mtx_get_value(&DRvi, 1, 0) + mtx_get_value(&wi, 2, 0)*mtx_get_value(&DRvi, 2, 0);
		tmpd3 = mtx_get_value(&DRvi, 0, 0)*mtx_get_value(&tmp, 0, 0) + mtx_get_value(&DRvi, 1, 0)*mtx_get_value(&tmp, 1, 0) + mtx_get_value(&DRvi, 2, 0)*mtx_get_value(&tmp, 2, 0);
		mtx_set_value(J, i, 4, 2 * (tmpd1*tmpd2 - tmpd3*witwi));
		mtx_multiply(&DRvi, &DRc, &vi);
		tmpd2 = mtx_get_value(&wi, 0, 0)*mtx_get_value(&DRvi, 0, 0) + mtx_get_value(&wi, 1, 0)*mtx_get_value(&DRvi, 1, 0) + mtx_get_value(&wi, 2, 0)*mtx_get_value(&DRvi, 2, 0);
		tmpd3 = mtx_get_value(&DRvi, 0, 0)*mtx_get_value(&tmp, 0, 0) + mtx_get_value(&DRvi, 1, 0)*mtx_get_value(&tmp, 1, 0) + mtx_get_value(&DRvi, 2, 0)*mtx_get_value(&tmp, 2, 0);
		mtx_set_value(J, i, 5, 2 * (tmpd1*tmpd2 - tmpd3*witwi));
	}

}

void center_shift_solve(mtx * x, mtx * v, mtx * w)
{
	double diff;
	mtx xn, tmpx, F, J;
	mtx_new(&xn, 6, 1);
	mtx_new(&tmpx, 6, 1);
	mtx_new(&F, 6, 1);
	mtx_new(&J, 6, 6);

	for (int itr = 0; itr < max_itr_times; itr++) {
		for (int i = 0; i < 6; i++) {
			mtx_set_value(&tmpx, i, 0, mtx_get_value(&xn, i, 0));
		}
		eva_DF(&J, &xn, v, w);
		eva_F(&F, &xn, v, w);
		
		mtx_solve_inverse(&J, &J);
		mtx_multiply(&F, &J, &F);
		mtx_subtract(&xn, &xn, &F);

		mtx_subtract(&F, &xn, &tmpx);
		diff = mtx_Fnorm(&F);
		if (diff < 1e-5) {
			printf("itr times = %d\n", itr);
			break;
		}
		if (itr == max_itr_times) {
			printf("itr times = %d\n", itr);
		}
	}
	
	for (int i = 0; i < 6; i++) {
		mtx_set_value(x, i, 0, mtx_get_value(&xn, i, 0));
	}



	mtx_delete(&xn);
	mtx_delete(&tmpx);
	mtx_delete(&F);
	mtx_delete(&J);
}
