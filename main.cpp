#include "mtx.h"
#include <stdio.h>
#include <string>
#include <math.h>
#include "stewart_platform_kinematics.h"

#define PI 3.14159265359
#define FOR(i,n0,n) for(int i=n0;i<n;i++)

int main() {

	// test stewart platform kinematics
	double thetaP = 109.007266910599f;
	double thetaB = 8.53112053424833f;
	double radiusP = 939.60275858147f;
	double radiusB = 1210.01377400165f;
	double Lengthsmin = 1220;
	double stroke = 323;

	SPK spk;
	SPK_new(&spk, thetaP, thetaB, radiusP, radiusB, Lengthsmin, stroke);

	mtx_print(&spk.parametersP);

	SPK_delete(&spk);



	system("pause");
}