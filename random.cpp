#include "random.h"
#define PI 3.14159265359

double random_uniform()
{
	return (double)rand() / RAND_MAX;
}

double random_uniform(double a, double b)
{
	double ans;
	ans = (double)rand() / RAND_MAX*(b - a);
	ans = ans + a;
	return ans;
}

double random_normal()
{
	srand(time(0));
	double u1 = (double)rand() / RAND_MAX;
	double u2 = (double)rand() / RAND_MAX;
	double ans;
	ans = sqrt(-2 * log(u1))*cos(2 * PI*u2);
	return ans;
}

double random_normal(double mean, double std)
{
	if (std <= 0) return 0;

	double ans = random_normal();
	ans = ans*std + mean;

	return ans;
}

