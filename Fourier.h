#pragma once

#include <complex>
#include <vector>

template <typename T>
std::complex<T> exponentialFactor(int m, unsigned int k, unsigned int n);

namespace trivialFourierTransform {
	template <typename T>
	std::vector<std::complex<T>> directTransform(const std::vector<T>& x, unsigned int k);

	template <typename T>
	std::vector<T> inverseTransform(const std::vector<std::complex<T>>& c, unsigned int n);
}

namespace fastFourierTransform {
	template <typename T>
	std::vector<std::complex<T>> directTransform(std::vector<T> x, unsigned int k);

	template <typename T>
	std::vector<T> inverseTransform(std::vector<std::complex<T>> c, unsigned int n);
}

namespace {
	constexpr double pi = 3.141592653589793238;
}

template <typename T>
std::complex<T> exponentialFactor(int m, unsigned int k, unsigned int n) {
	return std::exp(std::complex<T>(0, -1) * (pi * 2 * m * k / n));
}

template <typename T>
std::vector<std::complex<T>> trivialFourierTransform::directTransform(const std::vector<T>& x, unsigned int k) {
	std::vector<std::complex<T>> result(k, 0);
	for (unsigned int i = 0; i < k; i++) {
		for (int j = 0; j < x.size(); j++) {
			result[i] += x[j] * exponentialFactor<T>(j, i, x.size());
		}
		result[i] /= x.size();
	}
	return result;
}

template <typename T>
std::vector<T> trivialFourierTransform::inverseTransform(const std::vector<std::complex<T>>& c, unsigned int n) {
	std::vector<T> result(n, 0);
	for (unsigned int i = 0; i < n; i++) {
		for (int j = 0; j < c.size(); j++) {
			std::complex<T>&& temp = c[j] * exponentialFactor<T>(-j, i, n);
			result[i] += temp.real();
		}
	}
	return result;
}
