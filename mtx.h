#ifndef _MTX_H_
#define _MTX_H_

#include <stdio.h>
#include <string>
#include <math.h>

class vec3
{
public:
	vec3() { x = y = z = 0; };
	vec3(float a, float b, float c) { x = a; y = b; z = c; };

	vec3 operator+(vec3 v) { vec3 ans; ans.x = x + v.x; ans.y = y + v.y; ans.z = z + v.z; return ans; };
	vec3 operator-(vec3 v) { vec3 ans; ans.x = x - v.x; ans.y = y - v.y; ans.z = z - v.z; return ans; };
	void operator*=(float v) { x *= v; y *= v; z *= v; };
	void operator/=(float v) { x /= v; y /= v; z /= v; };
	void normalize()
	{
		float tmp = sqrt(x*x + y*y + z*z);
		x /= tmp; y /= tmp; z /= tmp;
	};
	float x, y, z;
};



struct Line
{
	vec3 v[2];
};

struct Plane
{
	vec3 p[4];
};

struct mtx
{
	double *v;
	int row, col;
};

inline void mtx_set_value(mtx *m, int i, int j, double v);
inline double mtx_get_value(mtx *m, int i, int j);

void mtx_set_zero(mtx *m);
void mtx_set_identity(mtx *m);
void mtx_multiply(mtx *m, mtx *m1, mtx *m2);
void mtx_add(mtx *m, mtx *m1, mtx *m2);
void mtx_subtract(mtx *m, mtx *m1, mtx *m2);

void mtx_new(mtx *m, int row, int col);
void mtx_delete(mtx *m);

void mtx_print(mtx *m);

void mtx_transpose(mtx *m);
void mtx_swap_row(mtx *m, int row1, int row2);
void mtx_equal(mtx *m1, mtx *m2);
void mtx_solve_inverse(mtx *m, mtx *im);

double mtx_elements_sum(mtx *m);
void mtx_normalize(mtx *m);


double mtx_Fnorm(mtx *m);
void mtx_rotation(mtx *m, double theta, int n);
void mtx_rotation3(mtx *m, mtx *shift);


#endif //_MTX_H_