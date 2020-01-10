#ifndef _CENTER_SHIFT_H_
#define _CENTER_SHIFT_H_

#include "mtx.h"

void pts_shift(mtx *w, mtx *v, mtx *x);

void eva_F(mtx * F, mtx *x, mtx *v, mtx *w);
void eva_DF(mtx * J, mtx *x, mtx *v, mtx *w);

void center_shift_solve(mtx *x, mtx *v, mtx *w);

#endif
