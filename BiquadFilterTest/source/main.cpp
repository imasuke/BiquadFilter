//main.cpp
#include "BiquadFilter.h"
#include "audio.hpp"

#include <iostream>
#include <vector>

using std::vector;
using namespace BiquadFilter;

int main(void){
	Audio audio;
	audio.open("white_noise.wav");

	vector<double> data;
	data.reserve(audio.length());
	for (unsigned int i = 0; i < audio.length(); i++){
		data.push_back(audio[i]);
	}

	LSFilter lsf(0.2, 1.0, 12.0);
	lsf.filtering(&data);

	for (unsigned int i = 0; i < audio.length(); i++){
		audio[i] = data[i];
	}
	audio.save("out.wav");

	return 0;
}