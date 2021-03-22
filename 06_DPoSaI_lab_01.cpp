#include <cmath>
#include <iostream>

#include "Fourier.h"

typedef double T;
constexpr double pi = 3.141592653589793238;

constexpr unsigned int sampleNum = 32;
constexpr T stride = pi * 2 / sampleNum;
constexpr unsigned int deepth = sampleNum;

int main()
{
	std::vector<T> functionSamples(sampleNum, 0);

	T argument = 0;
	for (unsigned int i = 0; i < sampleNum; i++) {
		functionSamples[i] = cos(argument * 3) + sin(argument * 2);
		argument += stride;
	}

	std::cout << "Input function: " << std::endl;
	for (unsigned int i = 0; i < functionSamples.size(); i++) {
		std::cout << i + 1 << ". " << functionSamples[i] << std::endl;
	}

	{
		using namespace Fourier::trivial;
		std::vector<std::complex<T>> frequencies =
			directTransform(functionSamples, deepth);
		//std::vector<T> restored = inverseTransform(frequencies, sampleNum);

		std::cout << std::endl << "Direct Fourier transform: " << std::endl;
		for (unsigned int i = 0; i < frequencies.size(); i++) {
			std::cout << i + 1 << ". " << frequencies[i].real() << " + " << frequencies[i].imag() << "*i" << std::endl;
		}
		//std::cout << std::endl << "Function restoration with inverse Fourier transform: " << std::endl;
		//for (unsigned int i = 0; i < restored.size(); i++) {
		//	std::cout << i + 1 << ". " << restored[i] << std::endl;
		//}
	}
	{
		using namespace Fourier::fast;
		std::vector<std::complex<T>> frequencies = directTransform(functionSamples, sampleNum);
		std::cout << std::endl << "Direct fast Fourier transform: " << std::endl;
		for (unsigned int i = 0; i < frequencies.size(); i++) {
			std::cout << i + 1 << ". " << frequencies[i].real() << " + " << frequencies[i].imag() << "*i" << std::endl;
		}
	}
}