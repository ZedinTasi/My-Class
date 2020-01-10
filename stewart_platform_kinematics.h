#pragma once
#ifndef _STEWART_PLATFORM_KINEMATICS_H_
#define _STEWART_PLATFORM_KINEMATICS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mtx.h"


// stewart platform kinematics
struct SPK {
	double angleP, angleB, radiusP, radiusB, Lengthsmin, stroke;
	mtx parametersP;
	mtx parametersB;

	mtx CylindersLengths;
	mtx Posture;

	// the vectors from 0 to the 6 points at the motion platform
	mtx CylindersVectors;
};

void SPK_update_parametersPB(SPK *spk);

// angle: degrees, radius: mm
void SPK_new(SPK *spk, double aP, double aB, double rP, double rB, 
	double Lmin, double str);
void SPK_delete(SPK *spk);

void SPK_inverse_kinematics(SPK *spk, mtx *posture);
void SPK_forward_kinematics(SPK *spk, mtx *CL);




#endif //_STEWART_PLATFORM_KINEMATICS_H_