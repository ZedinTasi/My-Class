#include "jinc.h"

#include <math.h>
#define _USE_MATH_DEFINES // for C
#include <stdio.h>

#define r1 72362614232.0f
#define r2 -7895059235.0f
#define r3 242396853.1f
#define r4 -2972611.439f
#define r5 15704.48260f
#define r6 -30.16036606f

#define s1 144725228442.0f
#define s2 2300535178.0f
#define s3 18583304.74f
#define s4 99447.43394f
#define s5 376.9991397f
#define s6 1.0f

#define p1 1.0f
#define p2 0.183105e-2f
#define p3 -0.3516396496e-4f
#define p4 .2457520174e-5f
#define p5 -.240337019e-6f

#define q1 0.04687499995f
#define q2 -.2002690873e-3f
#define q3 .8449199096e-5f
#define q4 -.88228987e-6f
#define q5 .105787412e-6f


double bessj1(double x) {
	double ans, y, ax, xx, z, sgn;
	ax = fabs(x);
	if (ax < 8.0f)
	{
		y = x*x;
		ans = x*(r1 + y*(r2 + y*(r3 + y*(r4 + y*(r5 + y*r6)))))
			/(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))));
	}
	else
	{
		z = 8 / ax;
		y = z*z;
		xx = ax - 2.356194491f; //ax - pi*3/4
		sgn = (ax>1e-8)? 1 / ax : 0.0f;
		ans = sqrt(0.636619772f / ax)*(cos(xx)*(p1 + y*(p2 + y*(p3 + y*(p4 + y*p5)))) - z*sin(xx)*(q1 + y*(q2 + y*(q3 + y*(q4 + y*q5)))))*sgn;
	}
	return ans;

}


double jinc(double x, double v)
{
	//return bessj1(x) / x;
	double ans, y, ax, xx, z, sgn;
	ax = fabs(x);
	if (ax < v)
	{
		y = x*x;
		ans = (r1 + y*(r2 + y*(r3 + y*(r4 + y*(r5 + y*r6)))))
			/ (s1 + y*(s2 + y*(s3 + y*(s4 + y*(s5 + y*s6)))));
	}
	else
	{
		z = v / ax;
		y = z*z;
		xx = ax - 2.356194491f; //ax - pi*3/4
		sgn = (ax>1e-8) ? 1 / ax : 0.0f;
		ans = sqrt(0.636619772f / ax)*(cos(xx)*(p1 + y*(p2 + y*(p3 + y*(p4 + y*p5)))) - z*sin(xx)*(q1 + y*(q2 + y*(q3 + y*(q4 + y*q5)))))*sgn;
	}
	return ans;

}

//int factorial(int n)
//{
//	if (n < 0) return 0;
//	else if (n < 2) return 1;
//	
//	int ans = 1;
//	for (int i = 2; i <= n; i++) {
//		ans *= i;
//	}
//	return ans;
//}
//double Bessel1(double x)
//{
//	double ans = 0, coef;
//	for (int n = 0; n <= 17; n++) {
//		coef = 1 / (pow(2, 2*n + 1)*(double)(factorial(n)*factorial(1 + n)));
//		if (n % 2 == 1) coef *= -1;
//		ans += coef*pow(x, 2*n + 1);
//	}
//	return ans;
//}
//double jinc2(double x)
//{
//	return Bessel1(x) / x;
//}