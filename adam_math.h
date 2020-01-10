#ifndef _ADAM_MATH_H_
#define _ADAM_MATH_H_

#include <math.h>

template <typename T> int sgn(T val) {
	return (T(0) <val) - (val <T(0));
}

template <typename T> T MAX(T a, T b) {
	return (a > b) ? a : b;
}
template <typename T> T MIN(T a, T b) {
	return (a > b) ? b : a;
}




#endif //_ADAM_MATH_H_