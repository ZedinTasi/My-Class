#include "mtx.h"


#define PI 3.14159265359
#define DEG2RAD(x) ((x)*0.01745329251)
#define RAD2DEG(x) ((x)*57.2957795131)

#define SWAP(type, x, y) {type swaptemp; swaptemp = x; x = y; y = swaptemp;}

inline void mtx_set_value(mtx *m, int i, int j, double v)
{
	m->v[i*m->col + j] = v;
}

inline double mtx_get_value(mtx *m, int i, int j)
{
	return m->v[i*m->col + j];
}

void mtx_print(mtx *m)
{
	printf("\n");
	for (int i = 0; i < m->row; i++)
	{
		for (int j = 0; j < m->col; j++)
		{
			printf("%14.8f ", mtx_get_value(m, i, j));
		}
		printf("\n");
	}
}

void mtx_set_zero(mtx *m)
{
	memset(m->v, 0, m->row*m->col * sizeof(double));
}

void mtx_transpose(mtx * m)
{
	double *v1 = (double*)malloc(m->row*m->col * sizeof(double));
	for (int i = 0; i < m->col; i++) {
		for (int j = 0; j < m->row; j++) {
			v1[i*m->row + j] = mtx_get_value(m, j, i);
		}
	}
	int tmp = m->row;
	m->row = m->col;
	m->col = tmp;
	free(m->v);
	m->v = v1;
}

void mtx_swap_row(mtx *m, int row1, int row2)
{
	double tmp;
	for (int i = 0; i < m->col; i++)
	{
		tmp = mtx_get_value(m, i, row1);
		mtx_set_value(m, i, row1, mtx_get_value(m, i, row1));
		mtx_set_value(m, i, row2, tmp);
	}
}

void mtx_new(mtx *m, int row, int col)
{
	m->row = row;
	m->col = col;
	m->v = (double*)malloc(m->row*m->col * sizeof(double));
	mtx_set_zero(m);
}

void mtx_delete(mtx *m)
{
	free(m->v);
}

void mtx_set_identity(mtx *m)
{
	int n;
	mtx_set_zero(m);
	n = (m->col < m->row) ? m->col : m->row;
	for (int i = 0; i < n; i++) mtx_set_value(m, i, i, 1.0);
}

void mtx_equal(mtx *m1, mtx *m2)
{
	memcpy(m1->v, m2->v, m1->col*m1->row * sizeof(double));
}

#define F_ABS(x) ((x > 0) ? (x) : -(x))


void mtx_multiply(mtx *m, mtx *m1, mtx *m2)
{
	// m1*m2 = m
	if (m->row != m1->row || m->col != m2->col || m1->col != m2->row) return;

	double v;
	mtx mm;
	mtx_new(&mm, m->row, m->col);

	for (int i = 0; i < m->row; i++)
	{
		for (int j = 0; j < m->col; j++)
		{
			v = 0.0;
			for (int k = 0; k < m1->col; k++)
			{
				v += mtx_get_value(m1, i, k)*mtx_get_value(m2, k, j);
			}
			mtx_set_value(&mm, i, j, v);
		}
	}
	for (int i = 0; i < m->row; i++)
	{
		for (int j = 0; j < m->col; j++)
		{
			mtx_set_value(m, i, j, mtx_get_value(&mm, i, j));
		}
	}

	mtx_delete(&mm);
}

void mtx_add(mtx * m, mtx * m1, mtx * m2)
{
	if (m->row != m1->row || m->row != m2->row || m->col != m1->col || m->col != m2->col) return;
	for (int i = 0; i < m->row; i++) {
		for (int j = 0; j < m->col; j++) {
			mtx_set_value(m, i, j, mtx_get_value(m1, i, j) + mtx_get_value(m2, i, j));
		}
	}
}

void mtx_subtract(mtx * m, mtx * m1, mtx * m2)
{
	if (m->row != m1->row || m->row != m2->row || m->col != m1->col || m->col != m2->col) return;
	for (int i = 0; i < m->row; i++) {
		for (int j = 0; j < m->col; j++) {
			mtx_set_value(m, i, j, mtx_get_value(m1, i, j) - mtx_get_value(m2, i, j));
		}
	}
}

void mtx_solve_inverse(mtx *m, mtx *im)
{
	mtx A, B, C;
	int size, *order;
	size = m->col;
	A = *m;
	mtx_new(&B, m->row, m->col);
	//mtx_new(&C, m->row, m->col);
	C = *im;
	mtx_equal(&B, &A);
	mtx_set_identity(&C);

	order = (int*)malloc(size * sizeof(double));
	for (int i = 0; i < size; i++) order[i] = i;

	int ii;
	for (int i = 0; i < size; i++)
	{
		ii = i + 1;
		while ((fabs(mtx_get_value(&B, i, i)) < 1e-8) && (ii < size))
		{
			//SWAP(double*, B.v[i], B.v[ii]);
			//SWAP(double*, C.v[i], C.v[ii]);
			mtx_swap_row(&B, i, ii);
			mtx_swap_row(&C, i, ii);
			SWAP(int, order[i], order[ii]);
			ii++;
		}

		for (int j = i + 1; j < size; j++)
		{
			double tmp = mtx_get_value(&B, j, i) / mtx_get_value(&B, i, i);
			for (int k = 0; k < size; k++)
			{
				mtx_set_value(&B, j, k,
					mtx_get_value(&B, j, k) - mtx_get_value(&B, i, k)*tmp);
				mtx_set_value(&C, j, k,
					mtx_get_value(&C, j, k) - mtx_get_value(&C, i, k)*tmp);
			}
		}
	}

	for (int i = size - 1; i >= 0; i--)
	{
		for (int j = i - 1; j >= 0; j--)
		{
			//double temp = B.v[j][i] / B.v[i][i];
			double tmp = mtx_get_value(&B, j, i) / mtx_get_value(&B, i, i);
			for (int k = size - 1; k >= 0; k--)
			{
				mtx_set_value(&B, j, k,
					mtx_get_value(&B, j, k) - mtx_get_value(&B, i, k)*tmp);
				mtx_set_value(&C, j, k,
					mtx_get_value(&C, j, k) - mtx_get_value(&C, i, k)*tmp);
			}
		}
	}

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			//C.v[i][j] /= B.v[i][i];
			mtx_set_value(&C, i, j, mtx_get_value(&C, i, j) / mtx_get_value(&B, i, i));
		}
	}

	for (int i = 0; i < size; i++)
	{
		while (order[i] != i)
		{
			//SWAP(double*, C.v[order[i]], C.v[i]);
			mtx_swap_row(&C, order[i], i);
			SWAP(int, order[order[i]], order[i]);
		}
	}

	mtx_delete(&B);
	//mtx_delete(&C);
	free(order);
	//return C;
}

double mtx_elements_sum(mtx * m)
{
	double sum = 0;
	for (int i = 0; i < m->col; i++) {
		for (int j = 0; j < m->row; j++) {
			sum += mtx_get_value(m, i, j);
		}
	}
	return sum;
}

void mtx_normalize(mtx * m)
{
	double sum = mtx_elements_sum(m);
	for (int i = 0; i < m->col; i++) {
		for (int j = 0; j < m->row; j++) {
			mtx_set_value(m, i, j, mtx_get_value(m, i, j) / sum);
		}
	}

}

double mtx_Fnorm(mtx * m)
{
	double ans = 0.0f;
	for (int i = 0; i < m->row; i++) {
		for (int j = 0; j < m->col; j++) {
			ans += mtx_get_value(m, i, j)*mtx_get_value(m, i, j);
		}
	}
	ans = sqrt(ans);
	return ans;
}

void mtx_rotation(mtx * m, double theta, int n)
{
	if (m->col != 3 || m->row != 3) return;

	switch (n)
	{
	case 0:
		mtx_set_value(m, 0, 0, cos(theta));
		mtx_set_value(m, 0, 1, -sin(theta));
		mtx_set_value(m, 1, 0, sin(theta));
		mtx_set_value(m, 1, 1, cos(theta));
		mtx_set_value(m, 2, 2, 1);
		break;
	case 1:
		mtx_set_value(m, 2, 2, cos(theta));
		mtx_set_value(m, 2, 0, -sin(theta));
		mtx_set_value(m, 0, 2, sin(theta));
		mtx_set_value(m, 0, 0, cos(theta));
		mtx_set_value(m, 1, 1, 1);
		break;
	case 2:
		mtx_set_value(m, 1, 1, cos(theta));
		mtx_set_value(m, 1, 2, -sin(theta));
		mtx_set_value(m, 2, 1, sin(theta));
		mtx_set_value(m, 2, 2, cos(theta));
		mtx_set_value(m, 0, 0, 1);
	}
}

void mtx_rotation3(mtx * m, mtx * shift)
{
	if (m->row != 3 || m->col != 3)return;
	if (shift->col != 1)return;

	int n = shift->row;
	mtx Ra, Rb, Rc;
	mtx_new(&Ra, 3, 3);
	mtx_new(&Rb, 3, 3);
	mtx_new(&Rc, 3, 3);
	mtx_rotation(&Ra, mtx_get_value(shift, n - 3, 0), 0);
	mtx_rotation(&Rb, mtx_get_value(shift, n - 2, 0), 1);
	mtx_rotation(&Rc, mtx_get_value(shift, n - 1, 0), 2);
	mtx_multiply(m, &Ra, &Rb);
	mtx_multiply(m, m, &Rc);
	mtx_delete(&Ra);
	mtx_delete(&Rb);
	mtx_delete(&Rc);
}
