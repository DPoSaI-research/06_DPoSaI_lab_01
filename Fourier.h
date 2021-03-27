#pragma once

#include <complex>
#include <vector>
#include <valarray>

#include <cmath>

namespace Fourier {
	template <typename T>
	std::complex<T> exponentialFactor(int m, unsigned int n);

	namespace trivial {
		template <typename T>
		std::vector<std::complex<T>> directTransform(const std::vector<T>& x, unsigned int k);

		template <typename T>
		std::vector<T> inverseTransform(const std::vector<std::complex<T>>& c, unsigned int n);
	}

	namespace fast {
		template <typename T>
		std::vector<std::complex<T>> directTransform(const std::vector<T>& x, unsigned int k);

		// Recursive optimization overload
		template <typename T>
		std::vector<std::complex<T>> directTransform(const std::vector<T>& x);

		// Should be the fastest version
		template <typename T>
		std::vector<std::complex<T>> directTransformNonrecursive(const std::vector<T>& x);

		template <typename T>
		std::vector<T> inverseTransform(std::vector<std::complex<T>> c, unsigned int n);

		namespace {
			// Base version
			template <typename T>
			std::valarray<std::complex<T>> recursiveDirect(std::valarray<T> source) {
				std::valarray<std::complex<T>> result(source.size());
				if (result.size() > 1) {
					size_t halfLen = source.size() / 2;
			
					std::valarray<std::complex<T>> factor(1, halfLen);
					std::complex<T> exp = Fourier::exponentialFactor<T>(1, source.size());
					for (unsigned int i = 1; i < factor.size(); i++) {
						factor[i] = factor[i - 1] * exp;
					}
			
					std::valarray<std::complex<T>>&& tempOdd = recursiveFourier<T>(std::valarray<T>(source[std::slice(1, halfLen, 2)])) * factor;
					std::valarray<std::complex<T>>&& tempEven = recursiveFourier<T>(std::valarray<T>(source[std::slice(0, halfLen, 2)]));
			
					result[std::slice(0, halfLen, 1)] = tempEven + tempOdd;
					result[std::slice(halfLen, halfLen, 1)] = tempEven - tempOdd;
				}
				else {
					result[0] = std::complex<T>(source[0], 0);
				}
				return result;
			}

			// Optimized version, pass dist initialized with funciton values c[i] = x[i];
			template <typename T>
			std::gslice recursiveDirect(std::valarray<std::complex<T>>& dist, std::slice source) {
				if (source.size() > 1) {
					size_t halfLen = source.size() / 2;

					std::slice inputEven(source.start(), halfLen, source.stride() * 2);
					std::slice inputOdd(source.start() + source.stride(), halfLen, source.stride() * 2);

					std::gslice leftPart = tempMatrix<T>(dist, inputEven);
					std::gslice rightPart = tempMatrix<T>(dist, inputOdd);


					std::valarray<std::complex<T>> factor(-1, halfLen);
					std::complex<T> exp = Fourier::exponentialFactor<T>(1, source.size());
					for (unsigned int i = 1; i < factor.size(); i++) {
						factor[i] = factor[i - 1] * exp;
					}

					dist[rightPart] *= factor;
					dist[leftPart] -= dist[rightPart];
					dist[rightPart] += dist[rightPart];
					dist[rightPart] += dist[leftPart];

					std::valarray<size_t> outputSizes(leftPart.size().size() + 1);
					outputSizes[std::slice(1, leftPart.size().size(), 1)] = leftPart.size();
					outputSizes[0] = 2;

					std::valarray<size_t> outputStrides(leftPart.stride().size() + 1);
					outputStrides[std::slice(1, leftPart.stride().size(), 1)] = leftPart.stride();
					outputStrides[0] = source.stride();

					return std::gslice(leftPart.start(), outputSizes, outputStrides);
				}
				else {
					return std::gslice(source.start(), std::valarray<size_t>({ source.size() }), std::valarray<size_t>({ source.stride() }));
				}
			}

			template <typename T>
			std::valarray<T> recursiveInverse(std::valarray<std::complex<T>> source);
		}
	}

	namespace {
		constexpr double pi = 3.141592653589793238;
	}
}

template <typename T>
std::complex<T> Fourier::exponentialFactor(int m, unsigned int n) {
	return std::exp(std::complex<T>(0, -1) * (pi * 2 * m / n));
}

template <typename T>
std::vector<std::complex<T>> Fourier::trivial::directTransform(const std::vector<T>& x, unsigned int k) {
	std::vector<std::complex<T>> result(k, 0);
	for (unsigned int i = 0; i < k; i++) {
		for (int j = 0; j < x.size(); j++) {
			result[i] += x[j] * exponentialFactor<T>(j * i, x.size());
		}
		result[i] /= x.size();
	}
	return result;
}

template <typename T>
std::vector<T> Fourier::trivial::inverseTransform(const std::vector<std::complex<T>>& c, unsigned int n) {
	std::vector<T> result(n, 0);
	for (unsigned int i = 0; i < n; i++) {
		for (int j = 0; j < c.size(); j++) {
			std::complex<T>&& temp = c[j] * exponentialFactor<T>(-j * i, n);
			result[i] += temp.real();
		}
	}
	return result;
}

template <typename T>
std::vector<std::complex<T>> Fourier::fast::directTransform(const std::vector<T>& x, unsigned int k) {
	std::valarray<std::complex<T>> result = Fourier::fast::recursiveDirect<T>(std::valarray<T>(x.data(), x.size()));
	result /= x.size(); // something extremely strange
	return std::vector<std::complex<T>>(std::begin(result), std::end(result));
}

template <typename T>
std::vector<T> Fourier::fast::inverseTransform(std::vector<std::complex<T>> c, unsigned int n) {
	std::valarray<T> result = Fourier::fast::recursiveInverse<T>(std::valarray<std::complex<T>>(c.data(), c.size()));
	return std::vector<T>(std::begin(result), std::end(result));
}

template <typename T>
std::vector<std::complex<T>> Fourier::fast::directTransform(const std::vector<T>& x) {
	std::valarray<std::complex<T>> result(x.size());
	for (unsigned int i = 0; i < x.size(); i++) {
		result[i] = std::complex<T>(x[i], 0);
	}
	result /= x.size();
	std::slice inputValues(0, result.size(), 1);
	std::valarray<std::complex<T>> ordered = result[Fourier::fast::recursiveDirect(result, inputValues)];
	return std::vector<std::complex<T>>(std::begin(ordered), std::end(ordered));
}


template <typename T>
std::vector<std::complex<T>> Fourier::fast::directTransformNonrecursive(const std::vector<T>& x) {
	size_t const inputSize = x.size();
	if ((inputSize & (inputSize - 1)) != 0) {
		throw std::logic_error("Invalid input vector size. ");
	}
	std::valarray<std::complex<T>> values(inputSize);
	for (unsigned int i = 0; i < inputSize; i++) {
		values[i] = std::complex<T>(x[i], 0);
	}
	values /= inputSize;

	// get binary logarithm for approximately O(1)
	int dimensionNum;
	frexp(inputSize, &dimensionNum);
	dimensionNum--;

	std::valarray<size_t> sizes(2, dimensionNum + 1);
	sizes[dimensionNum] = 1;

	std::valarray<size_t> strides(dimensionNum + 1);
	strides[dimensionNum] = inputSize;
	for (int i = dimensionNum - 1; i >= 0; --i) {
		strides[i] = strides[i + 1] / 2;
	}

	std::gslice generalSlice(0, sizes, strides);

	std::valarray<std::complex<T>> factor(-1, inputSize / 2);
	std::complex<T> exp = Fourier::exponentialFactor<T>(1, inputSize);
	for (unsigned int i = 1; i < factor.size(); i++) {
		factor[i] = factor[i - 1] * exp;
	}
	
	for (int i = dimensionNum - 1; i >= 0; --i) {
		std::slice subarraySelection(i + 1, dimensionNum - i, 1);
		std::slice factorSlice(0, inputSize / generalSlice.stride()[i] / 2, generalSlice.stride()[i]);
		for (size_t j = 0; j < generalSlice.stride()[i]; j++) {
			std::gslice leftSlice(j, generalSlice.size()[subarraySelection], generalSlice.stride()[subarraySelection]);
			std::gslice rightSlice(j + generalSlice.stride()[i], generalSlice.size()[subarraySelection], generalSlice.stride()[subarraySelection]);
	
			values[rightSlice] *= factor[factorSlice];
			values[leftSlice] -= values[rightSlice];
			// tricky way to multiply by 2 without unnecessary memory transactions, the same as values[rightSlice] *= 2;
 			// unfortunately there is no way to multiply slice by scalar
			values[rightSlice] += values[rightSlice]; 
			values[rightSlice] += values[leftSlice];
		}
	}
	std::valarray<std::complex<T>> ordered = values[generalSlice];
	return std::vector<std::complex<T>>(std::begin(ordered), std::end(ordered));
}



// size_t totalRangeSize = source.size()[source.size().size() - 1] << source.size().size() - 1;
// size_t totalRangeSize = 1;
// for (auto val: source.size()) {
// 	totalRangeSize *= val;
// }

//std::slice outputLow = std::slice(output.start(), halfLen, output.stride());
//std::slice outputHigh = std::slice(output.start() + halfLen * output.stride(), halfLen, output.stride());

//std::valarray<size_t> outputSizes(leftPart.size().size() + 1);
//outputSizes[std::slice(0, leftPart.size().size(), 1)] = leftPart.size();
//outputSizes[leftPart.size().size() - 1] = 2;
//outputSizes[leftPart.size().size()] = halfLen;
//
//std::valarray<size_t> outputStrides(leftPart.stride().size() + 1);
//outputStrides[std::slice(0, leftPart.stride().size(), 1)] = leftPart.stride();
//outputStrides[leftPart.stride().size()] = outputStrides[leftPart.stride().size() - 1] * 2;

//template <typename T>
//std::valarray<std::complex<T>> recursiveFourier(std::valarray<T> source) {
//	std::valarray<std::complex<T>> result(source.size());
//	if (result.size() > 1) {
//		size_t halfLen = source.size() / 2;
//
//		std::valarray<std::complex<T>> factor(1, halfLen);
//		std::complex<T> exp = Fourier::exponentialFactor<T>(1, source.size());
//		for (unsigned int i = 1; i < factor.size(); i++) {
//			factor[i] = factor[i - 1] * exp;
//		}
//
//		std::valarray<std::complex<T>>&& tempOdd = recursiveFourier<T>(std::valarray<T>(source[std::slice(1, halfLen, 2)])) * factor;
//		std::valarray<std::complex<T>>&& tempEven = recursiveFourier<T>(std::valarray<T>(source[std::slice(0, halfLen, 2)]));
//
//		result[std::slice(0, halfLen, 1)] = tempEven + tempOdd;
//		result[std::slice(halfLen, halfLen, 1)] = tempEven - tempOdd;
//	}
//	else {
//		result[0] = std::complex<T>(source[0], 0);
//	}
//	return result;
//}

//template <typename T>
//void temp(const std::valarray<T>& source, std::vector<std::complex<T>> dist, size_t start, size_t size, size_t stride) {
//	if (size > 1) {
//		size_t nextStride = stride * 2;
//		size_t nextSize = size / 2;
//
//		temp(source, dist, start, nextSize, nextStride);
//		temp(source, dist, start + stride, nextSize, nextStride);
//
//		std::complex<T> exp = Fourier::exponentialFactor<T>(1, size);
//		std::complex<T> factor = 1;
//
//		size_t halfLinear = stride * (size - 1) / 2;
//		for (size_t i = start; i < start + halfLinear; i += stride) {
//			dist[i + halfLinear] *= factor;
//			dist[i] -= dist[i + halfLinear];
//			dist[i + halfLinear] = dist[i] + dist[i + halfLinear] * 2;
//			factor *= exp;
//		}
//	}
//	else {
//		dist[start] = std::complex<T>(source[start], 0);
//	}
//}