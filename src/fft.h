#pragma once

#include "common.h"
#include "divisors.h"

template<typename T>
auto FFT(auto& data, T const& root) {
	std::vector<T> out(data.size());
	std::vector<T> rpow(data.size());
	rpow[0] = 1;
	for (size_t i = 1; i < rpow.size(); ++i) {
		rpow[i] = rpow[i-1];
		rpow[i] *= root;
	}
	auto const RPow = [&](auto const p) { return rpow[p%data.size()]; };
	auto const InternalFFT = [&](auto const& self, auto const& data, auto& out, auto const& order, auto dbeg, auto dend, auto in, auto ib, auto on, auto ob, auto pn, auto pa) -> void {
		auto dradix = dbeg++;
		if (dbeg == dend) {
			for (size_t i = 0; i < order; ++i) {
				for (size_t o = 0; o < order; ++o) {
					auto add = RPow(pn * o * i);
					add *= data[i * in + ib];
					out[on * o + ob] += add;
				}
			}
		} else {
			auto ib2 = ib;
			auto ob2 = ob;
			auto pa2 = 0;
			for (size_t sh = 0; sh < *dradix; ++sh) {
				self(self, data, out, order / *dradix, dbeg, dend, in * *dradix, ib2, on, ob2, pn * *dradix, pa2);
				ib2 += in;
				ob2 += on * order / *dradix;
				pa2 += pn;
			}
			std::vector<T> mid(*dradix);
			ob2 = ob;
			for (size_t but = 0; but < order / *dradix; ++but) {
				self(self, out, mid, *dradix, dradix, dbeg, on * order / *dradix, ob2, 1, 0, pn * order / *dradix, 0);
				for (size_t i = 0; i < *dradix; ++i) {
					out[on * order / *dradix * i + ob2] = mid[i];
					mid[i] = 0;
				}
				ob2 += on;
			}
		}
		if (pa != 0) {
			for (size_t o = 0; o < order; ++o) {
				out[on * o + ob] *= RPow(o * pa);
			}
		}
	};
	std::vector<size_t> divs;
	for (auto const& [d, c] : Factorization(data.size())) {
		for (int i = 0; i < c; ++i) {
			divs.push_back(d);
		}
	}
	InternalFFT(InternalFFT, data, out, data.size(), divs.begin(), divs.end(), 1, 0, 1, 0, 1, 0);
	return out;
}

#if 0
void TestFFT() {
	auto const NaiveDFT = []<typename T>(auto const& data, T root) {
		std::vector<T> out(data.size());
		for (u64 i = 0; i < data.size(); ++i) {
			auto rp = Power(root, i);
			T xp = Power(root, data.size());
			for (u64 j = 0; j < data.size(); ++j) {
				auto add = xp;
				add *= data[j];
				out[i] += add;
				xp *= rp;
			}
		}
		return out;
	};
	auto const CheckRes = [&]<typename... Args>(Args&&... args) {
		auto r1 = NaiveDFT(std::forward<Args>(args)...);
		auto r2 = FFT(std::forward<Args>(args)...);
		assert_(r1.size() == r2.size());
		for (size_t i = 0; i < r1.size(); ++i) {
			assert_(r1[i] == r2[i]);
		}
	};
	{
		m64::SetGlobalMod set(11);
		m64 r(PrimitiveRoot(11));
		r *= r;
		std::vector<u64> test{ 1,2,3,4,5 };
		CheckRes(test, r);
	}
	{
		m64::SetGlobalMod set(5);
		m64 r(PrimitiveRoot(5));
		std::vector<u64> test{ 1,2,3,4 };
		CheckRes(test, r);
	}
	{
		m64::SetGlobalMod set(17);
		m64 r(PrimitiveRoot(17));
		r *= r;
		std::vector<u64> test{ 1,2,3,4,5,6,7,8 };
		CheckRes(test, r);
	}
	{
		m64::SetGlobalMod set(19);
		m64 r(PrimitiveRoot(19));
		r *= r;
		std::vector<u64> test{ 1,2,3,4,5,6,7,8,9 };
		CheckRes(test, r);
	}
	{
		m64::SetGlobalMod set(1951);
		m64 r(PrimitiveRoot(1951));
		r *= r;
		std::vector<u64> test;
		for (int i = 1; i <= 975; ++i) {
			test.push_back(i);
		}
		CheckRes(test, r);
	}
	{
		m64::SetGlobalMod set(65537);
		std::vector<m64> v1(64), v2(64);
		for (int i = 0; i < 32; ++i) {
			v1[i] = 1 + i * i;
			v2[i] = 1 + 2 * i;
		}
		std::vector<m64> p1(64);
		for (int i = 0; i < 32; ++i) {
			for (int j = 0; j < 32; ++j) {
				auto add = v1[i];
				add *= v2[j];
				p1[i + j] += add;
			}
		}
		m64 r(Power<m64>(PrimitiveRoot(65537ull), 1024));
		auto f1 = FFT(v1, r);
		auto f2 = FFT(v2, r);
		std::vector<m64> f(64);
		for (int i = 0; i < 64; ++i) {
			f[i] = f1[i];
			f[i] *= f2[i];
		}
		auto p2 = FFT(f, Power(r, 63));
		for (int i = 0; i < 64; ++i) {
			p2[i] /= 64;
			assert_(p1[i] == p2[i]);
		}
	}
}

auto PerformanceTestFFT() {
	std::vector<u64> vals(1ull << 20);
	for (u64 i = 0; i < vals.size(); ++i) {
		vals[i] = i + 1 + ((i * i) >> 20);
	}
	u64 o = 1;
	u64 p;
	for (;;) {
		p = o * vals.size() + 1;
		if (IsPrime(p)) break;
		++o;
	}
	auto start = std::chrono::steady_clock::now();
	{
		m64::SetGlobalMod set(p);
		auto res = FFT(vals, Power<m64>(PrimitiveRoot(p), o));
	}
	return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
}
#endif