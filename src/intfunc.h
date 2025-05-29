#pragma once

#include"common.h"
#include"bits.h"

template<typename TRet = void, typename T>
constexpr SubVoid<TRet, T> Factorial(T const& x) {
    assert_(x >= 0);
    if (x == 0) return SubVoid<TRet, T>(1);
    SubVoid<TRet, T> res(x);
    for (T i = 2; i < x; ++i) res *= i;
    return res;
}

template<typename TRet = void, typename TBase, typename TExp>
constexpr SubVoid<TRet, TBase> Power(TBase const& base, TExp const& exp) {
    assert_(exp >= 0);
    SubVoid<TRet, TBase> res(base);
    if (exp == 0) {
        res = 1;
    } else {
        assert_(BitVal(exp, BitLength(exp) - 1));
        for (size_t i = BitLength(exp) - 1; i > 0;) {
            --i;
            res *= res;
            if (BitVal(exp, i)) res *= base;
        }
    }
    return res;
}

template<typename TRet>
constexpr TRet BinomialC(size_t k, size_t n) {
    assert_(k <= n);
    if (2 * k > n) return BinomialC<TRet>(n - k, n);
    if (k == 0) return TRet(1);
    if constexpr(false/*TODO*/) {
        std::vector<TRet> vec;
        vec.reserve(k);
        for (size_t i = 0; i < k; ++i) {
            vec.push_back(n - k + 1 + i);
        }
        auto o = vec[0] - 1;
        for (size_t d = 2; d <= k; ++d) {
            vec[d - 1 - o % d] /= d;
        }
        TRet res(vec[0]);
        for (size_t i = 1; i < k; ++i) {
            res *= vec[i];
        }
        return res;
    } else {
        TRet res(1);
		for (size_t i = 0; i < k; ++i) {
			res *= n - k + 1 + i;
		}
		for (size_t i = 2; i <= k; ++i) {
			res /= i;
		}
        return res;
    }
}

double BinomialCLog(size_t k, size_t n) {
    return std::lgamma(double(n + 1)) - std::lgamma(double(k + 1)) - std::lgamma(double(n - k + 1));
}

bool IsPow2(auto x) {
    return (x & (x - 1)) == 0;
}

template<typename T>
constexpr auto Digits(T num, int base = 10) {
    assert_(base > 0);
    std::vector<int> res;
    res.reserve(std::ceil(BitLength(num) / std::log2(base) + 1e-3));
    if constexpr(signed_int<T>) {
        if (num < 0) num = -num;
    }
    if (IsPow2(to_unsigned<int>(base))) {
        to_unsigned<int> d = 0;
        size_t dsz = BitLength(to_unsigned<int>(base)) - 1;
        size_t numsz = BitLength(num);
        int c = 0;
        for (size_t i = 0; i < numsz; ++i) {
            if (BitVal(num, i)) d += to_unsigned<int>(1) << c;
            if (c == dsz - 1) {
                res.emplace_back(d);
                d = 0;
                c = 0;
            } else ++c;
        }
        if (d > 0) res.emplace_back(d);
    } else {
        while (num != 0) {
            auto num_ = num;
            T rem;
            DivMod<T>(num_, base, num, rem);
            res.emplace_back(rem);
        }
    }
    return res;
}

template<typename TRet = unsigned long long, typename It>
constexpr auto FromDigits(It begin, It end, int base = 10) {
    TRet res(0);
    TRet exp(1);
    for (;begin != end; ++begin) {
        TRet add = exp;
        add *= *begin;
        res += add;
        exp *= base;
    }
    return res;
}

template<typename TRet = void, typename T>
constexpr auto IntSqrt(T const& x) {
    if (x == 0) return SubVoid<TRet, T>(0);
    assert_(x > 0);
    T ans = 1;
    ans <<= (BitLength(x) + 1) / 2;
    T d;
    for(;;) {
        d = x / ans;
        if (d > ans) {
            d -= ans;
            if (d == 1) break;
            d >>= 1;
            ans += d;
        } else if (d < ans) {
            d = ans - d;
            if (d == 1) {
                --ans;
                break;
            }
            d >>= 1;
            ans -= d;
        } else break;
    }
    return SubVoid<TRet, T>(ans);
}

template<typename TRet = void, typename T>
constexpr auto IntNrt(T const& x, size_t n) {
    if (x == 0) return SubVoid<TRet, T>(0);
    if (n%2 == 0) {
        assert_(x > 0);
    } else if constexpr (signed_int<T>) {
        if (x < 0) return -IntNrt<TRet>(n, -x);
    }
    T ans = 1;
    ans <<= (BitLength(x) + n - 1) / n;
    T d;
    for(;;) {
        d = x / Power(ans, n-1);
        if (d > ans) {
            d -= ans;
            if (d < n) break;
            d /= n;
            ans += d;
        } else if (d < ans) {
            d = ans - d;
            if (d < n) {
                --ans;
                break;
            }
            d /= n;
            ans -= d;
        } else break;
    }
    return SubVoid<TRet, T>(ans);
}

template<typename TRet = void, typename T>
constexpr std::optional<SubVoid<TRet, T>> IfIntPower(T const& x, size_t p) {
    auto rt = IntNrt(x, p);
    if (Power(rt, p) == x) {
        return SubVoid<TRet, T>(rt);
    } else {
        return {};
    }
}

template<typename TRet>
TRet BinomialCCached(size_t k, size_t n) {
    assert_(k <= n);
    if (k == 0 || k == n) return TRet(1);
    static std::unordered_map<size_t, TRet> cache;
    size_t index = n*(n - 3)/2 + k;
    if (auto it = cache.find(index); it != cache.end()) {
        return it->second;
    } else {
        TRet res = BinomialCCached<TRet>(k, n - 1);
        res += BinomialCCached<TRet>(k - 1, n - 1);
        cache[index] = res;
        return res;
    }
}

template<typename T = unsigned long long>
constexpr void ForAllPalindromes(auto const& func, int const base = 10) {
    auto const pbase = [base](auto const i) { return Power<T>(base, i); };
	for (T i(1); i < 10; ++i) {
		func(i);
	}
	for (T i(11); i < 100; i += 11) {
		func(i);
	}
	for (unsigned int len = 3;; ++len) {
		auto p10c = pbase(len / 2);
        auto add0 = p10c + (len & 1 ? 0 : p10c / base);
		for (auto num = pbase(len - 1) + 1;;) {
            if (func(num)) return;
            auto numc = num / p10c;
            if (numc % base + 1 != base) {
                num += add0;
            } else {
				T p10i(base * base);
				while (numc % p10i + 1 == p10i) p10i *= 10;
				if (p10i > numc * base) break;
				if (len & 1) p10i /= base;
				num += p10c * (base + 1) / p10i;
            }
		}
	}
}

template<typename T = unsigned long long, typename TN>
constexpr T Fib(TN const& n, T a = 0, T b = 1) {
    struct FibPair {
        T a, b;

        constexpr FibPair (T vala, T valb = 0) : a(vala), b(valb) {}
        constexpr FibPair& operator=(T const& other) {
            a = other;
            b = 0;
            return *this;
        }

        constexpr FibPair& operator*=(FibPair const& other) {
			std::tie(a, b) = std::make_tuple(a * other.a + b * other.b, b*(other.a+other.b)+a*other.b);
            return *this;
        }
    };

    if (n == 0) return a;
    if (n == 1) return b;
    return Power(FibPair{a, b}, n).b;
}