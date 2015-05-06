//BiquadFilter.cpp

#include "BiquadFilter.h"
#include <cmath>

//using
using std::vector;

namespace BiquadFilter{
	Filter::~Filter(){}

	void Filter::filtering(vector<double> *x){
		vector<double> &in = *x;
		vector<double> out(in.size());
		double bin1 = 0.0, bin2 = 0.0;
		double bout1 = 0.0, bout2 = 0.0;

		for (unsigned int i = 0; i < out.size(); i++){
			out[i] = (b[0] / a[0]) * in[i] + (b[1] / a[0]) * bin1 + (b[2] / a[0]) * bin2 - (a[1] / a[0]) * bout1 - (a[2] / a[0]) * bout2;
			//update input buf
			bin2 = bin1;
			bin1 = in[i];
			//update output buf
			bout2 = bout1;
			bout1 = out[i];
		}

		//copy
		for (unsigned int i = 0; i < in.size(); i++){
			in[i] = out[i];
		}

	}

	void Filter::alloc(){
		a.resize(3);
		b.resize(3);
	}

}