#include "complex.h"
#define FOR(i,n0,n) for(int i=n0;i<n;i++)

complex polar2complex(double radius, double rad)
{
	return complex(radius*cos(rad), radius*sin(rad));
}

void mtxc_set_zero(mtxc * m)
{
	memset(m->v, 0, m->row * m->col * sizeof(complex));
}

void mtxc_new(mtxc * m, int row, int col)
{
	m->col = col;
	m->row = row;
	m->v = (complex*)malloc(row*col*sizeof(complex));
	mtxc_set_zero(m);
}

void mtxc_delete(mtxc * m)
{
	free(m->v);
}

inline complex mtxc_get_value(mtxc * m, int i, int j)
{
	return m->v[i*m->col + j];
}

void mtxc_set_value(mtxc * m, int i, int j, complex v)
{
	m->v[i*m->col + j] = v;
}

complex mtxc_elements_sum(mtxc * m)
{
	complex sum = 0.0;
	for (int i = 0; i < m->col; i++) {
		for (int j = 0; j < m->row; j++) {
			sum += mtxc_get_value(m, i, j);
		}
	}
	return sum;
}

void mtx2mtxc(mtx * m1, mtxc * m2)
{
	if (m1->row > m2->row || m1->col > m2->col)return;


	FOR(i, 0, m1->row) {
		FOR(j, 0, m1->col) {
			m2->v[i*m2->col + j] = mtx_get_value(m1, i, j);
		}
	}

}

void mtxc2mtx(mtxc * m1, mtx * m2)
{
	if (m1->row > m2->row || m1->col > m2->col)return;


	FOR(i, 0, m1->row) {
		FOR(j, 0, m1->col) {
			m2->v[i*m2->col + j] = mtxc_get_value(m1, i, j).real();
		}
	}

}

void mtxc_print(mtxc * m)
{
	printf("\n");
	for (int i = 0; i < m->row; i++)
	{
		for (int j = 0; j < m->col; j++)
		{
			mtxc_get_value(m, i, j).print();
			printf("\t");
		}
		printf("\n");
	}
}
