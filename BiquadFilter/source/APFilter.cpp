//APFilter.cpp
#define _USE_MATH_DEFINES

#include "BiquadFilter.h"
#include <cmath>

using std::vector;

namespace BiquadFilter{
	APFilter::APFilter(double cutoff, double Q){
		this->cutoff = cutoff;
		this->Q = Q;

		alloc();

		//init filter coefficient
		double omega = 2.0 * M_PI* cutoff;
		double alpha = sin(omega) / 2.0 * Q;

		a[0] = 1.0 + alpha;
		a[1] = -2.0 * cos(omega);
		a[2] = 1.0 - alpha;
		b[0] = 1.0 - alpha;
		b[1] = -2.0 * cos(omega);
		b[2] = 1.0 + alpha;
	}
}