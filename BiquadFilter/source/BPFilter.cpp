//BPFilter.cpp
#define _USE_MATH_DEFINES

#include "BiquadFilter.h"
#include <cmath>

using std::vector;

namespace BiquadFilter{
	BPFilter::BPFilter(double low_edge, double high_edge){
		this->low_edge = low_edge;
		this->high_edge = high_edge;

		alloc();

		//init filter coefficient
		double bw = log2(high_edge / low_edge);
		double cutoff = low_edge * pow(2, bw/2);
		double omega = 2.0 * M_PI* cutoff;
		double alpha = sin(omega) * sinh(log(2.0)) / 2.0 * bw * omega / sin(omega);

		a[0] = 1.0 + alpha;
		a[1] = -2.0 * cos(omega);
		a[2] = 1.0 - alpha;
		b[0] = alpha;
		b[1] = 0.0;
		b[2] = -alpha;
	}
}