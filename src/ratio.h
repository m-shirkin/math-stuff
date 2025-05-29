#pragma once

#include "common.h"
#include "divisors.h"

template<typename T>
struct Ratio;

template<typename T>
concept ratio = template_of<T, Ratio>;

template<typename T>
struct Ratio {
	T num;
	to_unsigned<T> denom;

	constexpr void Normalize() {
		if (num == 0) {
			denom = 1;
		} else {
			auto gcd = GCD(Abs(num), denom);
			num /= T(gcd);
			denom /= gcd;
		}
	}

	constexpr Ratio(T const& p = 0, T const& q = 1)
		: num(p), denom(q) 
	{
		assert_(q > 0);
		Normalize();
	}

	template<typename TO>
	constexpr Ratio& operator=(TO&& other) {
		if constexpr(ratio<TO>) {
			num = std::forward<TO>(other).num;
			denom = std::forward<TO>(other).denom;
		} else {
			num = std::forward<TO>(other);
			denom = 1;
		}
		return *this;
	}

	constexpr Ratio& operator+=(Ratio const& other) {
		auto gcd = GCD(other.denom, denom);
		num = T(other.denom / gcd) * num + other.num * T(denom / gcd);
		denom *= other.denom / gcd;
		Normalize();
		return *this;
	}

	constexpr Ratio& operator-=(Ratio const& other) {
		auto gcd = GCD(other.denom, denom);
		num = T(other.denom / gcd) * num - other.num * T(denom / gcd);
		denom *= other.denom / gcd;
		Normalize();
		return *this;
	}

	constexpr Ratio& operator*=(Ratio const& other) {
		if (num == 0) return *this;
		if (other.num == 0) {
			return *this = 0;
		}
		auto gcd1 = GCD(Abs(other.num), denom);
		auto gcd2 = GCD(Abs(num), other.denom);
		num /= T(gcd2);
		num *= other.num / T(gcd1);
		denom /= gcd1;
		denom *= other.denom / gcd2;
		return *this;
	}

	constexpr Ratio& operator/=(Ratio const& other) {
		assert_(other.num != 0);
		if (other.num < 0) num = -num;
		return *this *= Ratio(other.denom, Abs(other.num));
	}

	constexpr Ratio operator-() const& {
		return Ratio(-num, denom);
	}
};

template<typename OStrm, typename T>
auto& operator<<(OStrm& strm, Ratio<T> const& r) {
	return strm << r.num << "/" << r.denom;
}

template<typename T1, typename T2>
constexpr bool operator==(Ratio<T1> const& lhs, Ratio<T2> const& rhs) {
	return lhs.num == rhs.num && lhs.denom == rhs.denom;
}

#define DEFINE_RATIO_OP(op) \
    template<typename T1, typename T2> requires (ratio<T1> || ratio<T2>) \
    constexpr auto operator op(T1 const& t1, T2 const& t2) { \
        if constexpr(ratio<T1>) { \
            T1 t1c = t1; \
            t1c op##= T1(t2); \
            return t1c; \
        } else { \
            T2 t1c(t1); \
            t1c op##= t2; \
            return t1c; \
        } \
    }

DEFINE_RATIO_OP(+)
DEFINE_RATIO_OP(-)
DEFINE_RATIO_OP(*)
DEFINE_RATIO_OP(/)

#undef DEFINE_RATIO_OP
