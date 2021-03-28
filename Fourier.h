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
		std::vector<std::complex<T>> directTransform(const std::vector<T>& x);

		template <typename T>
		std::vector<T> inverseTransform(const std::vector<std::complex<T>>& c);

		namespace {
			template <typename T>
			std::valarray<std::complex<T>> recursiveDirect(const std::valarray<T>& source) {
				std::valarray<std::complex<T>> result(source.size());
				if (result.size() > 1) {
					size_t halfLen = source.size() / 2;

					std::valarray<std::complex<T>> factor(1, halfLen);
					std::complex<T> exp = Fourier::exponentialFactor<T>(1, source.size());
					for (unsigned int i = 1; i < factor.size(); i++) {
						factor[i] = factor[i - 1] * exp;
					}

					std::valarray<std::complex<T>>&& tempOdd = recursiveDirect<T>(source[std::slice(1, halfLen, 2)]) * factor;
					std::valarray<std::complex<T>>&& tempEven = recursiveDirect<T>(source[std::slice(0, halfLen, 2)]);

					result[std::slice(0, halfLen, 1)] = tempEven + tempOdd;
					result[std::slice(halfLen, halfLen, 1)] = tempEven - tempOdd;
				}
				else {
					result[0] = std::complex<T>(source[0], 0);
				}
				return result;
			}

			template <typename T>
			std::valarray<std::complex<T>> recursiveInverse(const std::valarray<std::complex<T>>& source) {
				std::valarray<std::complex<T>> result(source.size());
				if (result.size() > 1) {
					size_t halfLen = source.size() / 2;

					std::valarray<std::complex<T>> factor(1, halfLen);
					std::complex<T> exp = Fourier::exponentialFactor<T>(-1, source.size());
					for (unsigned int i = 1; i < factor.size(); i++) {
						factor[i] = factor[i - 1] * exp;
					}

					std::valarray<std::complex<T>>&& tempOdd = recursiveInverse<T>(source[std::slice(1, halfLen, 2)]) * factor;
					std::valarray<std::complex<T>>&& tempEven = recursiveInverse<T>(source[std::slice(0, halfLen, 2)]);

					result[std::slice(0, halfLen, 1)] = tempEven + tempOdd;
					result[std::slice(halfLen, halfLen, 1)] = tempEven - tempOdd;
				}
				else {
					result[0] = source[0];
				}
				return result;
			}
		}
	}
	namespace utils {
		template <typename T>
		std::vector<std::complex<T>> realToComplex(const std::vector<T>& input) {
			size_t inputSize = input.size();
			std::vector<std::complex<T>> result(inputSize);
			for (size_t i = 0; i < inputSize; ++i) {
				result[i] = std::complex<T>(input[i], 0);
			}
			return result;
		}

		template <typename T>
		std::vector<T> complexToReal(const std::vector<std::complex<T>>& input) {
			size_t inputSize = input.size();
			std::vector<T> result(inputSize);
			for (size_t i = 0; i < inputSize; ++i) {
				result[i] = input[i].real();
			}
			return result;
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
std::vector<std::complex<T>> Fourier::fast::directTransform(const std::vector<T>& x) {
	std::valarray<std::complex<T>> result = Fourier::fast::recursiveDirect<T>(std::valarray<T>(x.data(), x.size()));
	result /= x.size();
	return std::vector<std::complex<T>>(std::begin(result), std::end(result));
}

template <typename T>
std::vector<T> Fourier::fast::inverseTransform(const std::vector<std::complex<T>>& c) {
	std::valarray<std::complex<T>> result = Fourier::fast::recursiveInverse<T>(std::valarray<std::complex<T>>(c.data(), c.size()));
	return Fourier::utils::complexToReal(std::vector<std::complex<T>>(std::begin(result), std::end(result)));
}