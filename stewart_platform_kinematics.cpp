#include "stewart_platform_kinematics.h"

#define PI (3.14159265359f)
#define DEG2RAD(x) ((x)*0.01745331111)
#define RAD2DEG(x) ((x)*57.2957795131)
void stewart_kinematics_inverse_transform(mtx *CL, mtx *CLvector, mtx *posture, mtx *P, mtx *B);
void stewart_kinematics_get_Jaccobian(mtx *J, mtx *posture, mtx *P, mtx *B);
void stewart_kinematics_get_diff(mtx *diff, mtx *posture, mtx *CL, mtx *P, mtx *B);

void SPK_update_parametersPB(SPK * spk)
{
	double thetaP = spk->angleP;
	double thetaB = spk->angleB;
	double rP = spk->radiusP;
	double rB = spk->radiusB;
	
	double lamPi, lamBi;
	for (int i = 0; i < 3; i++) {
		lamPi = (2 * i + 1)*PI / 3 - thetaP / 2;
		lamBi = (2 * i + 1)*PI / 3 - thetaB / 2;
		mtx_set_value(&spk->parametersP, 0, 2 * i, rP*cos(lamPi));
		mtx_set_value(&spk->parametersP, 1, 2 * i, rP*sin(lamPi));
		mtx_set_value(&spk->parametersP, 0, 2 * i + 1, rP*cos(lamPi + thetaP));
		mtx_set_value(&spk->parametersP, 1, 2 * i + 1, rP*sin(lamPi + thetaP));
		mtx_set_value(&spk->parametersB, 0, 2 * i, rB*cos(lamBi));
		mtx_set_value(&spk->parametersB, 1, 2 * i, rB*sin(lamBi));
		mtx_set_value(&spk->parametersB, 0, 2 * i + 1, rB*cos(lamBi + thetaB));
		mtx_set_value(&spk->parametersB, 1, 2 * i + 1, rB*sin(lamBi + thetaB));
	}
}

void SPK_new(SPK * spk, double aP, double aB, double rP, double rB, double Lmin, double str)
{
	spk->angleP = DEG2RAD(aP);
	spk->angleB = DEG2RAD(aB);
	spk->radiusP = rP*0.001;
	spk->radiusB = rB*0.001;
	spk->Lengthsmin = Lmin*0.001;
	spk->stroke = str*0.001;

	mtx_new(&spk->parametersB, 3, 6);
	mtx_new(&spk->parametersP, 3, 6);
	mtx_new(&spk->CylindersLengths, 6, 1);
	mtx_new(&spk->Posture, 6, 1);
	mtx_new(&spk->CylindersVectors, 3, 6);

	SPK_update_parametersPB(spk);
}

void SPK_delete(SPK * spk)
{
	mtx_delete(&spk->parametersB);
	mtx_delete(&spk->parametersP);
	mtx_delete(&spk->CylindersLengths);
	mtx_delete(&spk->Posture);
	mtx_delete(&spk->CylindersVectors);
}

void SPK_inverse_kinematics(SPK * spk, mtx * posture)
{
	if (!posture || posture->row != 6 || posture->col != 1) return;

	mtx_copy(&spk->Posture, posture);
	double tmp;
	mtx R, pt, Pi, Bi, Li;
	mtx_new(&R, 3, 3);
	mtx_new(&pt, 3, 1);
	mtx_new(&Pi, 3, 1);
	mtx_new(&Bi, 3, 1);
	mtx_new(&Li, 3, 1);

	mtx_rotation3(&R, posture);
	mtx_set_value(&pt, 0, 0, mtx_get_value(posture, 0, 0));
	mtx_set_value(&pt, 1, 0, mtx_get_value(posture, 1, 0));
	mtx_set_value(&pt, 2, 0, mtx_get_value(posture, 2, 0));

	for (int i = 0; i < 6; i++) {
		mtx_set_value(&Pi, 0, 0, mtx_get_value(&spk->parametersP, 0, i));
		mtx_set_value(&Pi, 1, 0, mtx_get_value(&spk->parametersP, 1, i));
		mtx_set_value(&Bi, 0, 0, mtx_get_value(&spk->parametersB, 0, i));
		mtx_set_value(&Bi, 1, 0, mtx_get_value(&spk->parametersB, 1, i));

		//  calculate the length of legs
		mtx_multiply(&Li, &R, &Pi);
		mtx_add(&Li, &pt, &Li);
		for (int j = 0; j < 3; j++) {
			mtx_set_value(&spk->CylindersVectors, j, i, mtx_get_value(&Li, j, 0));
		}
		mtx_subtract(&Li, &Li, &Bi);

		tmp = mtx_Fnorm(&Li);
		mtx_set_value(&spk->CylindersLengths, i, 0, tmp);
	}

	mtx_delete(&R);
	mtx_delete(&pt);
	mtx_delete(&Pi);
	mtx_delete(&Bi);
	mtx_delete(&Li);
}

void SPK_forward_kinematics(SPK * spk, mtx * CL)
{
	if (!CL || CL->row != 6 || CL->col != 1) return;

	mtx_copy(&spk->CylindersLengths, CL);

	mtx x, tmp, J, diff;
	mtx_new(&x, 6, 1);
	mtx_new(&tmp, 6, 1);
	mtx_new(&J, 6, 6);
	mtx_new(&diff, 6, 1);

	//initial value x = [0,0,1,0,0,0]^t
	mtx_set_value(&x, 2, 0, 1.0f);

	int max_itr_step = 1000;
	int itr_step = 0;
	double eps = 1e-8;

	for (itr_step = 0; itr_step < max_itr_step; itr_step++) {
		// tmp = x
		for (int j = 0; j < 6; j++) {
			mtx_set_value(&tmp, j, 0, mtx_get_value(&x, j, 0));
		}

		//  x = x - J \ eva_diff(x,LL)
		stewart_kinematics_get_diff(&diff, &x, CL, &spk->parametersP, &spk->parametersB);
		stewart_kinematics_get_Jaccobian(&J, &x, &spk->parametersP, &spk->parametersB);
		mtx_solve_inverse(&J, &J);
		mtx_multiply(&diff, &J, &diff);
		mtx_subtract(&x, &x, &diff);

		// if ||tmp-x||^2<eps, then break;
		mtx_subtract(&tmp, &tmp, &x);
		if (mtx_Fnorm(&tmp) < eps) {
			break;
		}
	}

	for (int j = 0; j < 6; j++) {
		mtx_set_value(&spk->Posture, j, 0, mtx_get_value(&x, j, 0));
	}

	mtx_delete(&x);
	mtx_delete(&tmp);
	mtx_delete(&J);
	mtx_delete(&diff);
}

void stewart_kinematics_inverse_transform(mtx *CL, mtx *CLvector, mtx *posture, mtx *P, mtx *B) {
	if (!posture || posture->row != 6 || posture->col != 1) return;
	if (!CL || CL->row != 6 || CL->col != 1) return;
	if (!CLvector || CLvector->row != 3 || CLvector->col != 6) return;

	double tmp;
	mtx R, pt, Pi, Bi, Li;
	mtx_new(&R, 3, 3);
	mtx_new(&pt, 3, 1);
	mtx_new(&Pi, 3, 1);
	mtx_new(&Bi, 3, 1);
	mtx_new(&Li, 3, 1);

	mtx_rotation3(&R, posture);
	mtx_set_value(&pt, 0, 0, mtx_get_value(posture, 0, 0));
	mtx_set_value(&pt, 1, 0, mtx_get_value(posture, 1, 0));
	mtx_set_value(&pt, 2, 0, mtx_get_value(posture, 2, 0));

	for (int i = 0; i < 6; i++) {
		mtx_set_value(&Pi, 0, 0, mtx_get_value(P, 0, i));
		mtx_set_value(&Pi, 1, 0, mtx_get_value(P, 1, i));
		mtx_set_value(&Bi, 0, 0, mtx_get_value(B, 0, i));
		mtx_set_value(&Bi, 1, 0, mtx_get_value(B, 1, i));

		//  calculate the length of legs
		mtx_multiply(&Li, &R, &Pi);
		mtx_add(&Li, &pt, &Li);
		for (int j = 0; j < 3; j++) {
			mtx_set_value(CLvector, j, i, mtx_get_value(&Li, j, 0));
		}
		mtx_subtract(&Li, &Li, &Bi);

		tmp = mtx_Fnorm(&Li);
		mtx_set_value(CL, i, 0, tmp);
	}

	mtx_delete(&R);
	mtx_delete(&pt);
	mtx_delete(&Pi);
	mtx_delete(&Bi);
	mtx_delete(&Li);
}

void stewart_kinematics_get_Jaccobian(mtx * J, mtx * posture, mtx * P, mtx * B)
{
	if (!J || J->row != 6 || J->col != 6) return;
	if (!posture || posture->row != 6 || posture->col != 1) return;

	double a, b, c;
	mtx R, Ra, Rb, Rc, DRa, DRb, DRc, pt, Pi, Bi, Li, tmp33, tmp31;
	mtx_new(&R, 3, 3);
	mtx_new(&Ra, 3, 3);
	mtx_new(&Rb, 3, 3);
	mtx_new(&Rc, 3, 3);
	mtx_new(&DRa, 3, 3);
	mtx_new(&DRb, 3, 3);
	mtx_new(&DRc, 3, 3);
	mtx_new(&pt, 3, 1);
	mtx_new(&Pi, 3, 1);
	mtx_new(&Bi, 3, 1);
	mtx_new(&Li, 3, 1);
	mtx_new(&tmp33, 3, 3);
	mtx_new(&tmp31, 3, 1);

	a = mtx_get_value(posture, 3, 0);
	b = mtx_get_value(posture, 4, 0);
	c = mtx_get_value(posture, 5, 0);

	mtx_rotation3(&R, posture);
	mtx_set_value(&pt, 0, 0, mtx_get_value(posture, 0, 0));
	mtx_set_value(&pt, 1, 0, mtx_get_value(posture, 1, 0));
	mtx_set_value(&pt, 2, 0, mtx_get_value(posture, 2, 0));
	mtx_rotation(&Ra, a, 0);
	mtx_rotation(&Rb, b, 1);
	mtx_rotation(&Rc, c, 2);

	mtx_set_value(&DRa, 0, 0, -sin(a));
	mtx_set_value(&DRa, 0, 1, -cos(a));
	mtx_set_value(&DRa, 1, 0, cos(a));
	mtx_set_value(&DRa, 1, 1, -sin(a));
	mtx_set_value(&DRb, 0, 0, -sin(b));
	mtx_set_value(&DRb, 2, 0, -cos(b));
	mtx_set_value(&DRb, 0, 2, cos(b));
	mtx_set_value(&DRb, 2, 2, -sin(b));
	mtx_set_value(&DRc, 1, 1, -sin(c));
	mtx_set_value(&DRc, 1, 2, -cos(c));
	mtx_set_value(&DRc, 2, 1, cos(c));
	mtx_set_value(&DRc, 2, 2, -sin(c));

	for (int i = 0; i < 6; i++) {
		mtx_set_value(&Pi, 0, 0, mtx_get_value(P, 0, i));
		mtx_set_value(&Pi, 1, 0, mtx_get_value(P, 1, i));
		mtx_set_value(&Bi, 0, 0, mtx_get_value(B, 0, i));
		mtx_set_value(&Bi, 1, 0, mtx_get_value(B, 1, i));

		mtx_multiply(&Li, &R, &Pi);
		mtx_add(&Li, &pt, &Li);
		mtx_subtract(&Li, &Li, &Bi);
		mtx_set_value(J, i, 0, 2 * mtx_get_value(&Li, 0, 0));
		mtx_set_value(J, i, 1, 2 * mtx_get_value(&Li, 1, 0));
		mtx_set_value(J, i, 2, 2 * mtx_get_value(&Li, 2, 0));

		mtx_multiply(&tmp33, &DRa, &Rb);
		mtx_multiply(&tmp33, &tmp33, &Rc);
		mtx_multiply(&tmp31, &tmp33, &Pi);
		mtx_set_value(J, i, 3, 2 * mtx_inner_product(&Li, &tmp31));
		mtx_multiply(&tmp33, &Ra, &DRb);
		mtx_multiply(&tmp33, &tmp33, &Rc);
		mtx_multiply(&tmp31, &tmp33, &Pi);
		mtx_set_value(J, i, 4, 2 * mtx_inner_product(&Li, &tmp31));
		mtx_multiply(&tmp33, &Ra, &Rb);
		mtx_multiply(&tmp33, &tmp33, &DRc);
		mtx_multiply(&tmp31, &tmp33, &Pi);
		mtx_set_value(J, i, 5, 2 * mtx_inner_product(&Li, &tmp31));
	}

	mtx_delete(&R);
	mtx_delete(&Ra);
	mtx_delete(&Rb);
	mtx_delete(&Rc);
	mtx_delete(&DRa);
	mtx_delete(&DRb);
	mtx_delete(&DRc);
	mtx_delete(&pt);
	mtx_delete(&Pi);
	mtx_delete(&Bi);
	mtx_delete(&Li);
	mtx_delete(&tmp31);
	mtx_delete(&tmp33);
}

void stewart_kinematics_get_diff(mtx * diff, mtx * posture, mtx * CL, mtx * P, mtx * B)
{
	if (!posture || posture->row != 6 || posture->col != 1) return;
	if (!CL || CL->row != 6 || CL->col != 1) return;
	if (!diff || diff->row != 6 || diff->col != 1) return;

	double tmp;
	mtx L;
	mtx_new(&L, 3, 6);
	stewart_kinematics_inverse_transform(diff, &L, posture, P, B);
	for (int i = 0; i < 6; i++) {
		tmp = mtx_get_value(diff, i, 0)*mtx_get_value(diff, i, 0) - mtx_get_value(CL, i, 0)*mtx_get_value(CL, i, 0);
		mtx_set_value(diff, i, 0, tmp);
	}

	mtx_delete(&L);
}
