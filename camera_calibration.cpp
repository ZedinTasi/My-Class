#include "camera_calibration.h"

#define PI 3.14159265359
#define DEG2RAD(x) ((x)*0.01745329251)
#define RAD2DEG(x) ((x)*57.2957795131)
#define FOR(i,n0,n) for(int i=n0;i<n;i++)
#define MAX(a,b) (a>b)?a:b

void clb_new(clb *C, int n)
{
	C->num_of_pts = n;
	mtx_new(&C->coo, 4, n);
	mtx_new(&C->V, 2, n);
	mtx_new(&C->W, 3, n);

	mtx_new(&C->p0, 12, 1);
	mtx_new(&C->solp, 12, 1);

	mtx_new(&C->R, 3, 3);
	mtx_new(&C->orien, 6, 1);
}

void clb_delete(clb *C)
{
	mtx_delete(&C->V);
	mtx_delete(&C->W);
	mtx_delete(&C->coo);
	mtx_delete(&C->p0);
	mtx_delete(&C->solp);
	mtx_delete(&C->R);
	mtx_delete(&C->orien);
}

void clb_initialize(clb * C, double *coo, double * p0)
{
	clb_evaVW(C, coo);
	mtx_copy(&C->p0, p0);
}

void clb_set_cpt(clb * C, int sx, int sy, int px, int py, int count)
{
	if (count >= C->num_of_pts)return;
	mtx_set_value(&C->coo, 0, count, sx);
	mtx_set_value(&C->coo, 1, count, sy);
	mtx_set_value(&C->coo, 2, count, px);
	mtx_set_value(&C->coo, 3, count, py);
}

int clb_first_cpt_nonzero(clb * C)
{
	FOR(i, 0, C->num_of_pts) {
		if (mtx_get_value(&C->coo, 0, i) + mtx_get_value(&C->coo, 1, i) + mtx_get_value(&C->coo, 2, i) + mtx_get_value(&C->coo, 3, i) < 0.01)
			return i;
	}
	return C->num_of_pts;
}



void clb_evaf(clb * C, mtx* ans, mtx* p)
{
	int n = C->num_of_pts;
	double Xi, Yi, Zi;
	double nume_x, nume_y, deno;
	

	FOR(i, 0, n) {
		Xi = mtx_get_value(&C->W, 0, i); Yi = mtx_get_value(&C->W, 1, i); Zi = mtx_get_value(&C->W, 2, i);
		nume_x = mtx_get_value(p, 0, 0)*Xi
			+ mtx_get_value(p, 1, 0)*Yi
			+ mtx_get_value(p, 2, 0)*Zi
			+ mtx_get_value(p, 3, 0);
		nume_y = mtx_get_value(p, 4, 0)*Xi
			+ mtx_get_value(p, 5, 0)*Yi
			+ mtx_get_value(p, 6, 0)*Zi
			+ mtx_get_value(p, 7, 0);
		deno = mtx_get_value(p, 8, 0)*Xi
			+ mtx_get_value(p, 9, 0)*Yi
			+ mtx_get_value(p, 10, 0)*Zi
			+ mtx_get_value(p, 11, 0);
		mtx_set_value(ans, 2 * i, 0, nume_x / deno - mtx_get_value(&C->V, 0, i));
		mtx_set_value(ans, 2 * i+1, 0, nume_y / deno - mtx_get_value(&C->V, 1, i));
	}
	mtx r1, r2, r3;
	mtx_new(&r1, 3, 1);
	mtx_new(&r2, 3, 1);
	mtx_new(&r3, 3, 1);
	FOR(i, 0, 3) {
		mtx_set_value(&r1, i, 0, mtx_get_value(p, i, 0));
		mtx_set_value(&r2, i, 0, mtx_get_value(p, i + 4, 0));
		mtx_set_value(&r3, i, 0, mtx_get_value(p, i + 8, 0));
	}

	mtx_set_value(ans, 2 * n, 0, mtx_inner_product(&r1, &r1) - 1);
	mtx_set_value(ans, 2 * n + 1, 0, mtx_inner_product(&r2, &r2) - 1);
	mtx_set_value(ans, 2 * n + 2, 0, mtx_inner_product(&r3, &r3) - 1);
	mtx_set_value(ans, 2 * n + 3, 0, mtx_inner_product(&r1, &r2));
	mtx_set_value(ans, 2 * n + 4, 0, mtx_inner_product(&r1, &r3));
	mtx_set_value(ans, 2 * n + 5, 0, mtx_inner_product(&r2, &r3));

	double det = mtx_get_value(&r1, 0, 0)*mtx_get_value(&r2, 1, 0)*mtx_get_value(&r3, 2, 0)
		+ mtx_get_value(&r1, 1, 0)*mtx_get_value(&r2, 2, 0)*mtx_get_value(&r3, 0, 0)
		+ mtx_get_value(&r1, 2, 0)*mtx_get_value(&r2, 0, 0)*mtx_get_value(&r3, 1, 0)
		- mtx_get_value(&r1, 0, 0)*mtx_get_value(&r2, 2, 0)*mtx_get_value(&r3, 1, 0)
		- mtx_get_value(&r1, 1, 0)*mtx_get_value(&r2, 0, 0)*mtx_get_value(&r3, 2, 0)
		- mtx_get_value(&r1, 2, 0)*mtx_get_value(&r2, 1, 0)*mtx_get_value(&r3, 0, 0);
		
	mtx_set_value(ans, 2 * n + 6, 0, det - 1);



	mtx_delete(&r1);
	mtx_delete(&r2);
	mtx_delete(&r3);
}

void clb_evaDf(clb * C, mtx * J, mtx * p)
{
	int n = C->num_of_pts;
	double Xi, Yi, Zi;
	double nume_x, nume_y, deno;
	double r11, r12, r13, r21, r22, r23, r31, r32, r33;
	r11 = mtx_get_value(p, 0, 0);
	r12 = mtx_get_value(p, 1, 0);
	r13 = mtx_get_value(p, 2, 0);
	r21 = mtx_get_value(p, 4, 0);
	r22 = mtx_get_value(p, 5, 0);
	r23 = mtx_get_value(p, 6, 0);	
	r31 = mtx_get_value(p, 8, 0);
	r32 = mtx_get_value(p, 9, 0);
	r33 = mtx_get_value(p, 10, 0);

	FOR(i, 0, n) {
		Xi = mtx_get_value(&C->W, 0, i); Yi = mtx_get_value(&C->W, 1, i); Zi = mtx_get_value(&C->W, 2, i);
		nume_x = mtx_get_value(p, 0, 0)*Xi
			+ mtx_get_value(p, 1, 0)*Yi
			+ mtx_get_value(p, 2, 0)*Zi
			+ mtx_get_value(p, 3, 0);
		nume_y = mtx_get_value(p, 4, 0)*Xi
			+ mtx_get_value(p, 5, 0)*Yi
			+ mtx_get_value(p, 6, 0)*Zi
			+ mtx_get_value(p, 7, 0);
		deno = mtx_get_value(p, 8, 0)*Xi
			+ mtx_get_value(p, 9, 0)*Yi
			+ mtx_get_value(p, 10, 0)*Zi
			+ mtx_get_value(p, 11, 0);
		mtx_set_value(J, 2 * i, 0, Xi / deno);
		mtx_set_value(J, 2 * i, 1, Yi / deno);
		mtx_set_value(J, 2 * i, 2, Zi / deno);
		mtx_set_value(J, 2 * i, 3, 1 / deno);
		mtx_set_value(J, 2 * i, 8, -nume_x*Xi / deno / deno);
		mtx_set_value(J, 2 * i, 9, -nume_x*Yi / deno / deno);
		mtx_set_value(J, 2 * i, 10, -nume_x*Zi / deno / deno);
		mtx_set_value(J, 2 * i, 11, -nume_x / deno / deno);

		mtx_set_value(J, 2 * i + 1, 4, Xi / deno);
		mtx_set_value(J, 2 * i + 1, 5, Yi / deno);
		mtx_set_value(J, 2 * i + 1, 6, Zi / deno);
		mtx_set_value(J, 2 * i + 1, 7, 1 / deno);
		mtx_set_value(J, 2 * i + 1, 8, -nume_y*Xi / deno / deno);
		mtx_set_value(J, 2 * i + 1, 9, -nume_y*Yi / deno / deno);
		mtx_set_value(J, 2 * i + 1, 10, -nume_y*Zi / deno / deno);
		mtx_set_value(J, 2 * i + 1, 11, -nume_y / deno / deno);
	}

	FOR(j, 0, 3) {
		mtx_set_value(J, 2 * n, j, 2 * mtx_get_value(p, j, 0));
		mtx_set_value(J, 2 * n + 1, j + 4, 2 * mtx_get_value(p, j + 4, 0));
		mtx_set_value(J, 2 * n + 2, j + 8, 2 * mtx_get_value(p, j + 8, 0));

		mtx_set_value(J, 2 * n + 3, j, mtx_get_value(p, j + 4, 0));
		mtx_set_value(J, 2 * n + 3, j + 4, mtx_get_value(p, j, 0));
		mtx_set_value(J, 2 * n + 4, j, mtx_get_value(p, j + 8, 0));
		mtx_set_value(J, 2 * n + 4, j + 8, mtx_get_value(p, j, 0));
		mtx_set_value(J, 2 * n + 5, j + 4, mtx_get_value(p, j + 8, 0));
		mtx_set_value(J, 2 * n + 5, j + 8, mtx_get_value(p, j + 4, 0));
	}
	mtx_set_value(J, 2 * n + 6, 0, r22*r33 - r23*r32);
	mtx_set_value(J, 2 * n + 6, 1, r23*r31 - r21*r33);
	mtx_set_value(J, 2 * n + 6, 2, r21*r32 - r22*r31);
	mtx_set_value(J, 2 * n + 6, 4, r13*r32 - r12*r33);
	mtx_set_value(J, 2 * n + 6, 5, r11*r33 - r13*r31);
	mtx_set_value(J, 2 * n + 6, 6, r12*r31 - r11*r32);
	mtx_set_value(J, 2 * n + 6, 8, r12*r23 - r13*r22);
	mtx_set_value(J, 2 * n + 6, 9, r13*r21 - r11*r23);
	mtx_set_value(J, 2 * n + 6, 10, r11*r22 - r12*r21);
}

double evarho(double F, double F_new, double u, mtx *hlm, mtx *Jtf) {
	double ans;
	mtx tmp;
	mtx_new(&tmp, 12, 1);

	mtx_scalar(&tmp, hlm, u);
	mtx_subtract(&tmp, &tmp, Jtf);

	ans = (F - F_new) * 2 / mtx_inner_product(hlm, &tmp);
	mtx_delete(&tmp);
	return ans;
}


void clb_LMlsqrsolve(clb * C)
{
	int n = C->num_of_pts;

	int itr_max = 50;
	double e1 = 1e-4, e2 = 1e-4;
	int v = 2;
	double u, F, F_new, rho;
	mtx p, f, J, Jt, JtJ, Jtf, uI;
	mtx diagJtJ, hlm, p_new, f_new;
	mtx_new(&p, 12, 1);
	mtx_new(&f, 2 * n + 7, 1);
	mtx_new(&J, 2 * n + 7, 12);
	mtx_new(&Jt, 12, 2 * n + 7);
	mtx_new(&JtJ, 12, 12);
	mtx_new(&Jtf, 12, 1);
	mtx_new(&uI, 12, 12);
	mtx_new(&diagJtJ, 12, 1);
	mtx_new(&hlm, 12, 1);
	mtx_new(&p_new, 12, 1);
	mtx_new(&f_new, 2 * n + 7, 1);

	mtx_copy(&p, &C->p0);
	clb_evaf(C, &f, &p);
	F = mtx_inner_product(&f, &f) / 2;
	clb_evaDf(C, &J, &p);
	mtx_transpose(&J, &Jt);
	mtx_multiply(&JtJ, &Jt, &J);
	mtx_multiply(&Jtf, &Jt, &f);

	mtx_diag(&JtJ, &diagJtJ);
	u = mtx_elements_max(&diagJtJ);

	FOR(itr, 0, itr_max) {
		mtx_set_identity(&uI);
		mtx_scalar(&uI, &uI, u);
		mtx_add(&uI, &JtJ, &uI);
		mtx_solve_inverse(&uI, &uI);
		mtx_scalar(&uI, &uI, -1.0);
		mtx_multiply(&hlm, &uI, &Jtf);

		if (mtx_Fnorm(&hlm) <= e2*(mtx_Fnorm(&p) + e2)){
			break;
		}

		mtx_add(&p_new, &p, &hlm);
		clb_evaf(C, &f_new, &p_new);
		F_new = mtx_inner_product(&f_new, &f_new) / 2;
		rho = evarho(F, F_new, u, &hlm, &Jtf);

		if (rho > 0)
		{
			mtx_copy(&p, &p_new);
			mtx_copy(&f, &f_new);
			F = F_new;
			clb_evaDf(C, &J, &p);
			mtx_transpose(&J, &Jt);
			mtx_multiply(&JtJ, &Jt, &J);
			mtx_multiply(&Jtf, &Jt, &f);
			if (mtx_norm_inf(&Jtf) < e1) {
				break;
			}
			u = MAX(1 / 3, 1 - pow(2 * rho - 1, 3));
			v = 2;
		}
		else
		{
			u = u*v; 
			v = 2 * v;
		}
	}

	double norm = sqrt(mtx_get_value(&p, 0, 0)*mtx_get_value(&p, 0, 0)
		+ mtx_get_value(&p, 1, 0)*mtx_get_value(&p, 1, 0)
		+ mtx_get_value(&p, 2, 0)*mtx_get_value(&p, 2, 0));
	mtx_scalar(&p, &p, 1.0 / norm);
	mtx_copy(&C->solp, &p);

	mtx_delete(&p);
	mtx_delete(&f);
	mtx_delete(&J);
	mtx_delete(&JtJ);
	mtx_delete(&Jtf);
	mtx_delete(&uI);
	mtx_delete(&diagJtJ);
	mtx_delete(&hlm);
	mtx_delete(&p_new);
	mtx_delete(&f_new);
}

void clb_evaOrien(clb * C)
{
	mtx Rt;
	mtx_new(&Rt, 3, 3);
	FOR(i, 0, 3) {
		mtx_set_value(&Rt, 0, i, mtx_get_value(&C->solp, i, 0));
		mtx_set_value(&Rt, 1, i, mtx_get_value(&C->solp, i + 4, 0));
		mtx_set_value(&Rt, 2, i, mtx_get_value(&C->solp, i + 8, 0));
	}
	mtx_transpose(&Rt, &C->R);

	//double r11, r12, r13, r21, r22, r23, r31, r32, r33;
	//r11 = mtx_get_value(&C->R, 0, 0);
	//r12 = mtx_get_value(&C->R, 0, 1);
	//r13 = mtx_get_value(&C->R, 0, 2);
	//r21 = mtx_get_value(&C->R, 1, 0);
	//r22 = mtx_get_value(&C->R, 1, 1);
	//r23 = mtx_get_value(&C->R, 1, 2);
	//r31 = mtx_get_value(&C->R, 2, 0);
	//r32 = mtx_get_value(&C->R, 2, 1);
	//r33 = mtx_get_value(&C->R, 2, 2);

	//double a, b, c, x, y, z;
	//x = -(r11*mtx_get_value(&C->solp, 3, 0)
	//	+ r12*mtx_get_value(&C->solp, 7, 0)
	//	+ r13*mtx_get_value(&C->solp, 11, 0));
	//y = -(r21*mtx_get_value(&C->solp, 3, 0)
	//	+ r22*mtx_get_value(&C->solp, 7, 0)
	//	+ r23*mtx_get_value(&C->solp, 11, 0));
	//z = -(r31*mtx_get_value(&C->solp, 3, 0)
	//	+ r32*mtx_get_value(&C->solp, 7, 0)
	//	+ r33*mtx_get_value(&C->solp, 11, 0));

	//c = atan(r32 / r33);
	//b = atan(-r31 / sqrt(r32*r32+ r33*r33));
	//a = atan(r21 / r11);

	//mtx_set_value(&C->orien, 0, 0, x);
	//mtx_set_value(&C->orien, 1, 0, y);
	//mtx_set_value(&C->orien, 2, 0, z);
	//mtx_set_value(&C->orien, 3, 0, a);
	//mtx_set_value(&C->orien, 4, 0, b);
	//mtx_set_value(&C->orien, 5, 0, c);

	clb_p2ori(&C->solp, &C->orien, 8.0);

	mtx_delete(&Rt);
}


void clb_evaCameraCoord(mtx *V, mtx *coo)
{
	int w = 5184;
	int h = 3456;
	int n = V->col;
	double angle = 97.166666667 / 2 / 180 * PI;
	double f = 2592 / tan(angle);

	FOR(j, 0, n) {
		mtx_set_value(V, 0, j, (mtx_get_value(coo, 2, j) - w / 2) / f);
		mtx_set_value(V, 1, j, (mtx_get_value(coo, 3, j) - h / 2) / f);
	}
}

void clb_evaDomeCoord(mtx *W, mtx *coo)
{
	double theta = PI / 24;
	double x[8] = { -0.5, -0.24625, 0.025, 0.29425, 0.54175, 0.749125, 0.899625, 1. };
	int n = W->col;
	double tmp;
	FOR(j, 0, n) {
		tmp = sqrt(1 - x[(int)mtx_get_value(coo, 1, j)] * x[(int)mtx_get_value(coo, 1, j)]);

		mtx_set_value(W, 0, j, tmp*sin(theta*mtx_get_value(coo, 0, j)));
		mtx_set_value(W, 1, j, x[(int)mtx_get_value(coo, 1, j)]);
		mtx_set_value(W, 2, j, tmp*cos(theta*mtx_get_value(coo, 0, j)));
	}
}

void clb_evaVW(clb * C, double * coo)
{
	mtx_copy(&C->coo, coo);
	clb_evaCameraCoord(&C->V, &C->coo);
	clb_evaDomeCoord(&C->W, &C->coo);
}

void clb_evaVW(clb * C)
{
	clb_evaCameraCoord(&C->V, &C->coo);
	clb_evaDomeCoord(&C->W, &C->coo);
}

void clb_p2ori(mtx *p, mtx *ori, double dome_radius)
{
	mtx R, pos;
	mtx_new(&R, 3, 3);
	mtx_new(&pos, 3, 1);
	FOR(i, 0, 3) {
		mtx_set_value(&R, 0, i, mtx_get_value(p, i, 0));
		mtx_set_value(&R, 1, i, mtx_get_value(p, i+4, 0));
		mtx_set_value(&R, 2, i, mtx_get_value(p, i+8, 0));
	}
	mtx_transpose(&R);
	mtx_set_value(&pos, 0, 0, mtx_get_value(p, 3, 0));
	mtx_set_value(&pos, 1, 0, mtx_get_value(p, 7, 0));
	mtx_set_value(&pos, 2, 0, mtx_get_value(p, 11, 0));
	mtx_multiply(&pos, &R, &pos);

	FOR(i, 0, 3) {
		mtx_set_value(ori, i, 0, -dome_radius*mtx_get_value(&pos, i, 0));
	}
	mtx_set_value(ori, 5, 0, RAD2DEG(atan(mtx_get_value(&R, 2, 1) / mtx_get_value(&R, 2, 2))));
	mtx_set_value(ori, 4, 0, RAD2DEG(atan(-mtx_get_value(&R, 2, 0) / 
		sqrt(mtx_get_value(&R, 2, 1)*mtx_get_value(&R, 2, 1) + mtx_get_value(&R, 2, 2)*mtx_get_value(&R, 2, 2)))));
	mtx_set_value(ori, 3, 0, RAD2DEG(atan(mtx_get_value(&R, 1, 0) / mtx_get_value(&R, 0, 0))));

	mtx_delete(&R);
	mtx_delete(&pos);
}

void clb_ori2p(mtx *ori, mtx *p, double dome_radius)
{
	mtx tmp, R;
	mtx_new(&tmp, 3, 1);
	mtx_new(&R, 3, 3);

	FOR(i, 0, 3) {
		mtx_set_value(&tmp, i, 0, DEG2RAD(mtx_get_value(ori, i+3, 0)) );
	}
	mtx_rotation3(&R, &tmp);
	mtx_transpose(&R);

	FOR(i, 0, 3) {
		mtx_set_value(&tmp, i, 0, mtx_get_value(ori, i, 0) / dome_radius);
	}
	mtx_multiply(&tmp, &R, &tmp);

	mtx_set_value(p, 3, 0, -mtx_get_value(&tmp, 0, 0));
	mtx_set_value(p, 7, 0, -mtx_get_value(&tmp, 1, 0));
	mtx_set_value(p, 11, 0, -mtx_get_value(&tmp, 2, 0));

	FOR(i, 0, 3) {
		mtx_set_value(p, i, 0, mtx_get_value(&R, 0, i));
		mtx_set_value(p, i+4, 0, mtx_get_value(&R, 1, i));
		mtx_set_value(p, i+8, 0, mtx_get_value(&R, 2, i));
	}

	mtx_delete(&R);
	mtx_delete(&tmp);
}
