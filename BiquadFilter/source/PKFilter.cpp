//PKFilter.cpp
#define _USE_MATH_DEFINES

#include "BiquadFilter.h"
#include <cmath>

using std::vector;

namespace BiquadFilter{
	PKFilter::PKFilter(double low_edge, double high_edge, double gain){
		this->low_edge = low_edge;
		this->high_edge = high_edge;
		this->gain = gain;

		alloc();

		//init filter coefficient
		double bw = log2(high_edge / low_edge);
		double cutoff = low_edge * pow(2, bw / 2);
		double omega = 2.0 * M_PI* cutoff;
		double alpha = sin(omega) * sinh(log(2.0)) / 2.0 * bw * omega / sin(omega);
		double A = pow(10.0, (gain/40.0));

		a[0] = 1.0 + alpha / A;
		a[1] = -2.0 * cos(omega);
		a[2] = 1.0 - alpha / A;
		b[0] = 1.0 + alpha * A;
		b[1] = -2.0 * cos(omega);
		b[2] = 1.0 - alpha * A;
	}
}