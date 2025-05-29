#pragma once

#include"common.h"
#include"bits.h"
#include"intfunc.h"

namespace GCD_detail {
    template<typename T, typename... Args>
    constexpr auto GCD_Impl(T a, T b, Args... args) {
        if (a < b) {
            std::swap(a, b);
            if (sizeof...(Args) > 0) {
                [&]<typename T_, T_... i>(std::integer_sequence<T_, i...>) {
                    (std::swap(select_nth<2 * i>(args...), select_nth<2 * i + 1>(args...)), ...);
                }(std::make_index_sequence<sizeof...(Args) / 2>());
            }
        }
        for(;;) {
            if (b > (a >> 1)) {
                T t = std::move(a);
                a = std::move(b);
                b = t - b;
                if (b == 0) break;
                if (sizeof...(Args) > 0) {
                    [&] <typename T_, T_... i>(std::integer_sequence<T_, i...>) {
                        ([&](){
                            auto t = select_nth<2 * i>(std::move(args)...);
                            select_nth<2 * i>(args...) = select_nth<2 * i + 1>(std::move(args)...);
                            select_nth<2 * i + 1>(args...) = t - select_nth<2 * i + 1>(std::move(args)...);
                        }(),...);
                    }(std::make_index_sequence<sizeof...(Args) / 2>());
                }
            } else {
                T q = a / b;
                {
                    T t = std::move(a);
                    a = std::move(b);
                    b = t - q * b;
                }
                if (b == 0) break;
                if (sizeof...(Args) > 0) {
                    [&] <typename T_, T_... i>(std::integer_sequence<T_, i...>) {
                        ([&](){
                            auto t = select_nth<2 * i>(std::move(args)...);
                            select_nth<2 * i>(args...) = select_nth<2 * i + 1>(std::move(args)...);
                            select_nth<2 * i + 1>(args...) = t - q*select_nth<2 * i + 1>(std::move(args)...);
                        }(),...);
                    }(std::make_index_sequence<sizeof...(Args) / 2>());
                }
            }
        }
        if constexpr(sizeof...(Args) == 0) {
            return std::move(a);
        } else {
            return [&]<typename T_, T_... i>(std::integer_sequence<T_, i...>) {
                return std::make_tuple(std::move(a), select_nth<2 * i + 1>(std::move(args)...)...);
            }(std::make_index_sequence<sizeof...(Args) / 2>());
        }
    }
}

struct Extended_tag {} extended_tag;
struct Inverse_tag {} inverse_tag;

template<typename T>
constexpr auto GCD(T a, T b) {
	assert_(a > 0);
	assert_(b > 0);
	IF_CHECKS(T a0 = a; T b0 = b;)
		return checkval(GCD_detail::GCD_Impl(std::move(a), std::move(b)),
			assert_(a0 % _ == 0);
	assert_(b0 % _ == 0);
		);
}

template<typename T>
constexpr auto LCM(T a, T b) {
    if (a == 0 || b == 0) {
        return T(0);
    }
    return a / GCD(a, b) * b;
}

template<typename TTag, signed_int T>
constexpr auto GCD(TTag, T a, T b) {
    assert_(a > 0);
    assert_(b > 0);
    IF_CHECKS(T a0 = a; T b0 = b;)
    if constexpr (std::is_same_v<TTag, Extended_tag>) {
        return checkval(GCD_detail::GCD_Impl(std::move(a), std::move(b), T(1), T(0), T(0), T(1)),
            assert_(a0 % std::get<0>(_) == 0);
            assert_(b0 % std::get<0>(_) == 0);
            assert_(a0 * std::get<1>(_) + b0 * std::get<2>(_) == std::get<0>(_));
        );
    } else if constexpr (std::is_same_v<TTag, Inverse_tag>) {
        return checkval(GCD_detail::GCD_Impl(std::move(a), std::move(b), T(0), T(1)),
            assert_(a0 % std::get<0>(_) == 0);
            assert_(b0 % std::get<0>(_) == 0);
            assert_((b0 * std::get<1>(_) - std::get<0>(_)) % a0 == 0);
        );
    }
}

template<typename T>
constexpr auto GeneralizedModularInverse(T a, T mod) {
    assert_(a > 0);
    assert_(mod > 0);
    auto [gcd, res] = GCD(inverse_tag, to_signed<T>(mod), to_signed<T>(a));
    if (res < 0) res += mod;
    assert_(res > 0);
    return std::make_pair(gcd, to_unsigned<T>(res));
}

template<typename T>
constexpr auto ModularInverse(T a, T mod) {
    auto [gcd, res] = GeneralizedModularInverse(a, mod);
    assert_(gcd == 1);
    return res;
}

template<bool CheckOverflow = true, typename T>
constexpr std::pair<T, std::optional<T>> ChineseRemainder(T const& p, T const& a, T const& q, T const& b) {
    if constexpr (CheckOverflow) {
        if (max(BitLength(p), BitLength(q)) * 2 > (bit_size<T> -signed_int<T>)) {
            using TExt = std::conditional_t<signed_int<T>, intN<bit_size<T> * 2>, uintN<bit_size<T> * 2>>;
            auto res = ChineseRemainder<false>(TExt(p), TExt(a), TExt(q), TExt(b));
            assert_(BitLength(res.first) <= bit_size<T>);
            return { T(res.first), res.second ? std::make_optional(T(*(res.second))) : std::optional<T>() };
        }
    }
    auto res = GCD(extended_tag, to_signed<T>(p), to_signed<T>(q));
    T gcd(std::get<0>(res));
    T lcm = p * q / gcd;
    if (a % gcd == b % gcd) {
		if (std::get<1>(res) < 0) std::get<1>(res) += q;
		if (std::get<2>(res) < 0) std::get<2>(res) += p;
        T ret = (b * T(std::get<1>(res)) % q * (p / gcd) + a * T(std::get<2>(res)) % p * (q / gcd)) % lcm;
        assert_(ret % p == a && ret % q == b);
        return {lcm, std::make_optional(ret)};
    } else {
        return {lcm, std::nullopt};
    }
}

template<typename T = unsigned long long, size_t wheelsz = 4>
struct PrimeGenerator {
private:
    static auto constexpr wheelp = [](){
        std::array<int, wheelsz> res;
        res[0] = 2;
        res[1] = 3;
        res[2] = 5;
        res[3] = 7;
        int cand = 11;
        for (size_t nxt = 4; nxt < wheelsz; ++nxt) {
            for (int i = 0; i < nxt; ++i) {
                if (cand % res[i] == 0) {
                    i = -1;
                    ++cand;
                }
            }
            res[nxt] = cand;
        }
        return res;
    }();

    static auto constexpr wheellen = [](){
        size_t res = 1;
        for (auto p : wheelp) res *= p;
        return res;
    }();

    static auto constexpr offsets = ArrayFromVector([](){
        std::array<bool, wheellen> isp;
        std::fill(isp.begin(), isp.end(), true);
        isp[0] = false;
        for (auto p : wheelp) {
            for (auto i = p;;i+=p) {
                if (i >= isp.size()) break;
                isp[i] = false;
            }
        }
        std::vector<int> res;
        for (size_t i = 0; i < isp.size(); ++i) {
            if (isp[i]) res.emplace_back(i);
        }
        return res;
    });

    static auto constexpr wheelcount = offsets.size();

    static auto constexpr wheelmap = []() {
        std::array<size_t, wheellen> res;
        for (size_t i = 0; i < wheellen; ++i) {
            res[i] = NPOS;
        }
        for (size_t i = 0; i < wheelcount; ++i) {
            res[offsets[i]] = i;
        }
        return res;
    }();

    std::vector<T> primebuf;

    std::vector<char> wheelmask;
    size_t wheelHi = 0;
    size_t wheelLo = 0;

    T startval = 0;
 
    constexpr char* MaskAtVal(T const& pos) {
        T quo, rem;
        DivMod<T>(pos-startval, wheellen, quo, rem);
        if (wheelmap[size_t(rem)] == NPOS) return nullptr;
        return std::addressof(wheelmask[size_t(quo) * wheelcount + wheelmap[size_t(rem)]]);
    }

    constexpr void NewWheels(size_t count) {
        startval += wheelHi * wheellen;
        auto endval = startval + count * wheellen;
        wheelHi = 0;
        wheelLo = 0;
        for (T np = (primebuf.empty() ? wheelp[wheelsz - 1] : primebuf.back()) + 1; np * np < endval; ++np) {
            if (wheelmap[size_t(np % wheellen)] != NPOS) {
                bool success = true;
                for (T const& p : primebuf) {
                    if (p * p > np) break;
                    if (np % p == 0) {
                        success = false;
                        break;
                    }
                }
                if (success) primebuf.push_back(np);
            }
        }
        wheelmask.resize(count * wheelcount);
        std::fill(wheelmask.begin(), wheelmask.end(), 0);
        for (T const& p : primebuf) {
            for (T i = max((startval + p - 1) / p * p, p * p); i < endval; i += p) {
                char* loc = MaskAtVal(i);
                if (loc) *loc = 1;
            }
        }
    }

public:
    constexpr PrimeGenerator() {
        NewWheels(1);
        wheelLo = wheelcount - 1;
        wheelHi = NPOS; // Magic value
    }

    constexpr T operator*() const {
        return wheelHi > NPOS - wheelsz
            ? wheelp[NPOS - wheelHi]
            : startval + wheelHi * wheellen + offsets[wheelLo];
    }

    constexpr void operator++() {
        do {
            ++wheelLo;
            if (wheelLo == wheelcount) {
                wheelLo = 0;
                if (wheelHi > NPOS - wheelsz) { // First wheelsz primes
                    wheelLo = wheelcount - 1;
                    --wheelHi;
                    if (wheelHi == NPOS - wheelsz) {
                        wheelHi = 0;
                        wheelLo = 1; // Skip 1 and wheel primes
                    }
                    break;
                } else if (++wheelHi * wheelcount == wheelmask.size()) {
                    NewWheels(size_t((2 * primebuf.back() + wheellen) / wheellen));
                }
            }
        } while (wheelmask[wheelHi * wheelcount + wheelLo]);
    }
};

template<typename T = unsigned long long>
T CachedPrime(size_t i) {
    static PrimeGenerator<T> pg;
    static std::vector<T> data;
    while (i >= data.size()) {
        data.push_back(*pg);
        ++pg;
    }
    return data[i];
}

#include"modint.h"

size_t constexpr TRIVIAL_PRIME_COUNT = 100;
auto constexpr FIRST_PRIMES = [](){
    std::array<size_t, TRIVIAL_PRIME_COUNT> res;
    PrimeGenerator<size_t> pg;
    for (int i = 0; i < TRIVIAL_PRIME_COUNT; ++i) {
        res[i] = *pg; ++pg;
    }
    return res;
}();

namespace detail {
    template<typename T>
    constexpr bool MillerRabin(T const& x, int k) {
        auto xm1 = x - 1;
        size_t s = 0;
        while (!BitVal(xm1, s)) ++s;
        auto d = xm1 >> s;
        auto a = 2ull;
        for (;k > 0; --k) {
            auto n = Power(ModInt(x, a), d);
            if (n.val != 1 && n.val != xm1) {
                for (size_t i = 0; i < s; ++i) {
                    auto nn = n;
                    nn *= nn;
                    n = nn;
                    if (n.val == 1) return false;
                    if (n.val == xm1) break;
                }
                if (n.val != xm1) return false;
            }
            do {
                a = a * 0x12345ull % 0x7fffffffull;
            } while ([amod = a%x, &xm1]() {return amod <= 1 || amod == xm1; }());
        }
        return true;
    }
}

template<typename T>
constexpr bool IsPrime(T const& x) {
    if (x < 2) return false;
    for (auto const& p : FIRST_PRIMES) {
        if (p*p > x) {
            return true;
        }
        if (x % p == 0) {
            return false;
        }
    }
	if (BitLength(x) * 2 <= bit_size<T>
		? detail::MillerRabin(x, 20)
		: detail::MillerRabin(uintN<bit_size<T>*2>(x), 20)
    ){
        return true;
    } else {
        return false;
    }
}

template<typename T1, typename T2>
using FactorizationT = std::vector<std::pair<T1, T2>>;

template<typename T>
constexpr auto Factorization(T x) {
    assert_(x >= 1);
    FactorizationT<T, int> res;
    if (x == 1) return res;
    auto CheckDiv = [&](auto p){
        if (p * p > x) return 1;
        if (x % p == 0) {
            res.emplace_back(p, 1);
            x /= p;
            while (x % p == 0) {
                ++res.back().second;
                x /= p;
            }
            if (p * p > x) return 1;
        }
        return 0;
    };
    for (auto const& p : FIRST_PRIMES) {
        if (CheckDiv(p)) {
            if (x > 1) {
                res.emplace_back(x, 1);
            }
            return res;
        }
    }
    if (x <= FIRST_PRIMES.back() * FIRST_PRIMES.back() || (
        BitLength(x) * 2 <= bit_size<T>
            ? detail::MillerRabin(x, 20)
            : detail::MillerRabin(uintN<bit_size<T>*2>(x), 20)
    )) {
        res.emplace_back(x, 1);
        return res;
    }
    auto dp = FIRST_PRIMES.back() % 6 == 1 ? 4 : 2;
    for (auto i = FIRST_PRIMES.back() + dp;; i += dp) {
        if (CheckDiv(i)) {
            if (x > 1) {
                res.emplace_back(x, 1);
            }
            return res;
        }
        dp = 6 - dp;
    }
}

template<typename T1, typename T2>
auto Divisors(FactorizationT<T1, T2> const& fact) {
    std::vector<T1> res;
    auto const iter = [&](auto const& self, T1 cur, int i) -> void {
        if (i < fact.size()) {
            self(self, cur, i + 1);
            for (int j = 0; j < fact[i].second; ++j) {
                cur *= fact[i].first;
                self(self, cur, i + 1);
            }
        } else {
            res.push_back(cur);
        }
    };
    iter(iter, T1(1), 0);
    if constexpr (requires (T1 x, T1 y){ x < y; }) std::sort(res.begin(), res.end());
    return res;
}

template<typename T>
auto Divisors(T const& x) {
    return Divisors(Factorization(x));
}

template<typename T1, typename T2>
auto DivisorSum(FactorizationT<T1, T2> const& fact) {
    T1 res(1);
    for (auto d : fact) {
        if (d.second == 1) {
            ++d.first;
            res *= d.first;
        } else {
            auto dp = Power<T1>(d.first, d.second + 1);
            --dp;
            --d.first;
            res *= dp;
            res /= d.first;
        }
    }
    return res;
}

template<typename T>
auto DivisorSum(T const& x) {
    return DivisorSum(Factorization(x));
}

namespace detail {
    // mode
    // 0: addition
    // 1: subtraction
    template<auto mode, typename T1, typename T2>
    auto MergeFactors_Impl(FactorizationT<T1, T2> const& x, FactorizationT<T1, T2> const& y) {
        FactorizationT<T1, T2> res;
        auto itx = x.begin();
        auto ity = y.begin();
        for (;;) {
            if (itx == x.end()) {
                for (; ity != y.end(); ++ity) {
                    if constexpr(mode == 0) {
                        res.push_back(*ity);
                    } else {
                        assert_(0);
                    }
                }
                break;
            }
            if (ity == y.end()) {
                for (; itx != x.end(); ++itx) {
                    res.push_back(*itx);
                }
                break;
            }
            if (itx->first < ity->first) {
                res.push_back(*itx);
                ++itx;
            } else if (itx->first == ity->first) {
                if constexpr (mode == 0) {
                    res.emplace_back(itx->first, itx->second + ity->second);
                } else {
                    assert_(itx->second >= ity->second);
                    if (itx->second > ity->second) {
                        res.emplace_back(itx->first, itx->second - ity->second);
                    }
                }
                ++itx;
                ++ity;
            } else {
                if constexpr (mode == 0) {
                    res.push_back(*ity);
                } else {
                    assert_(0);
                }
                ++ity;
            }
        }
        return res;
    }
}

template<typename... Args>
auto MergeFactors(Args&&... args) {
    return detail::MergeFactors_Impl<0>(std::forward<Args>(args)...);
}

template<typename... Args>
auto SubtractFactors(Args&&... args) {
    return detail::MergeFactors_Impl<1>(std::forward<Args>(args)...);
}

template<typename T>
constexpr T PrimitiveRoot(T mod, T searchFrom = T(2)) {
    assert_(mod > 2);
    auto factors = Factorization(mod - 1);
    for (;;) {
        if (!(searchFrom < mod)) searchFrom = T(2);
        if ([&]() {
            for (auto const& f : factors) {
                if (Power(ModInt(mod, searchFrom), (mod - 1) / f.first).val == 1) {
                    return false;
                }
            }
            return true;
        }()) return searchFrom;
        ++searchFrom;
    }
}

template<typename T1, typename T2 = int>
struct BulkFactorizer {
    BulkFactorizer(T1 val): sm(val) {
        for (T1 i = 2; i * i <= val; ++i) {
            if (sm[i-1] != 0) continue;
            for (T1 j = i * i; j <= val; j += i) {
                if (sm[j - 1] != 0) continue;
                sm[j - 1] = i;
            }
        }
    }

    auto Factorize(T1 x) const& {
        assert_(x <= sm.size());
        FactorizationT<T1, T2> res;
        while(x > 1) {
            auto d = sm[x-1] == 0 ? x : sm[x-1];
            if (res.size() > 0 && res.back().first == d) {
                ++res.back().second;
            } else {
                res.emplace_back(d, 1);                
            }
            x /= d;
        }
        return res;
    }

    std::vector<T1> sm;
};

namespace detail {
    template<typename T>
    struct PrimeCounter {
        static bool constexpr sh3[6] = { 0, 1, 1, 0, 0, 1 };
        static constexpr auto Ind(auto const& x) {
            return x / 3 + sh3[x % 6];
        }

        static constexpr auto ComputeSh(auto itp, auto itpe) {
            T cprod = 1;
            T mult = 1;
            std::vector<int> rems(1, T(0));
            for (; itp != itpe; ++itp) {
                rems.resize(cprod * *itp);
                for (size_t i = cprod * *itp; i != 0;) {
                    --i;
                    rems[i] = rems[i % cprod] + (i / cprod) * mult - rems[i / *itp];
                }
                mult *= *itp - 1;
                cprod *= *itp;
            }
            for (size_t i = 0; i < cprod; ++i) {
                rems[i] -= i * mult / cprod;
            }
            return std::make_pair(mult, rems);
        }

        std::vector<T> pi;
        std::vector<T> prs;
        constexpr auto const& Pi(auto const& x) const& {
            return pi[Ind(x)];
        };

        constexpr PrimeCounter(T const nlim) :
            pi(Ind(nlim) + 1)
        {
            pi[0] = 1;
            pi[1] = 2;
            prs = { 2,3 };
            for (size_t pb = 2; pb < pi.size(); ++pb) {
                T p = pb * 3 - 1 - (pb & 1);
                if (pi[pb] != T(-1)) {
                    pi[pb] = pi[pb - 1] + 1;
                    prs.push_back(p);
                    auto i0 = Ind((pb * 3 + 1 + (pb & 1)) * p);
                    auto d = i0 - Ind((pb * 3 - 1 - (pb & 1)) * p);
                    for (; i0 < pi.size(); i0 += 2 * p) {
                        pi[i0] = -1;
                        pi[i0 - d] = -1;
                    }
                    if (i0 - d < pi.size()) pi[i0 - d] = -1;
                } else {
                    pi[pb] = pi[pb - 1];
                }
            }
        }

        constexpr T Count(T const& n) const& {
            if (n == 2) return 1;
            if (Ind(n) < pi.size()) return Pi(n);
            T sqrtn = IntSqrt(n);
            T a = std::upper_bound(prs.begin(), prs.end(), IntNrt<3>(n) + 1) - prs.begin();
            T ans((Pi(sqrtn) - a + 1) * (Pi(sqrtn) + a - 2) / 2);
            {
                //std::vector<decltype(ComputeSh(prs.begin(), prs.begin()))> sh;
                //for (size_t i = 0; /*(sh.empty() ? 6 : sh.back().second.size()) * prs[i + 2] < prs[a]*/i < 1; ++i) {
                //    sh.push_back(ComputeSh(prs.begin(), prs.begin() + i + 3));
                //}
                auto const Phi = [&](auto const& self, auto const& x, auto o) -> to_signed<T> {
                    if (o == 1) {
                        return (x + 1) / 2;
                    } else if (o == 2) {
                        return Ind(x);
                    //} else if (o + 3 < sh.size()) {
                    //    auto const& sho = sh[o - 3];
                    //    return sho.first * x / sho.second.size() + sho.second[x % sho.second.size()];
                    } else if (x < prs[o] * prs[o] && Ind(x) < pi.size()) {
                        return 1 + (x < prs[o] ? 0 : Pi(x) - o);
                    } else {
                        to_signed<T> res = Ind(x);
                        for (--o; o > 1; --o) {
                            res -= self(self, x / prs[o], o);
                        }
                        return res;
                    }
                };
		        ans += Phi(Phi, n, a);
            }
            for (size_t b = a; b < Pi(sqrtn); ++b) {
                auto bdiv = Ind(n / prs[b]);
                if (bdiv < pi.size()) {
                    ans -= pi[bdiv];
                } else {
                    ans -= Count(n / prs[b]);
                }
            }
            return ans;
        }
    };
}

template<typename T>
auto constexpr PrimeCounter(T const& n) {
    return detail::PrimeCounter<T>(max(Power(IntNrt<3>(n), 2), 1000));
}

template<typename T>
auto constexpr PrimePi(T const& n) {
    return PrimeCounter(n).Count(n);
}

// Count over region G_lim defined by lim and sum function (size of G)
// "kx in G_kn iff x in G_n" should always hold
// dist_ function is to count how many values in [l, r) could be GCDs

template<typename T, typename TSum = std::type_identity<T>, typename TDist = std::type_identity<T>>
struct CoprimePairs {
private:
    template<typename TFn, auto... args>
    struct TFnRes {
        using type = decltype(std::declval<TFn>()(args...));
    };

	template<typename TFn, auto... args>
	struct TFnRes<std::type_identity<TFn>, args...> {
		using type = TFn;
	};

    using TSumRet = TFnRes<TSum, 0>::type;
    using TDistRet = TFnRes<TDist, 0, 0>::type;
    using TRet = decltype(std::declval<TSumRet>() * std::declval<TDistRet>());

	TSum sum_;
    TDist dist_;
    std::unordered_map<T, TRet> cache;

    auto sum(auto const n) {
		if constexpr (template_of<TSum, std::type_identity>) {
			TSumRet res(n / 2);
			res *= n - (n % 2 == 0);
			return res;
		} else {
			return sum_(n);
		}
    }

	auto dist(auto const l, auto const r) {
		if constexpr (template_of<TDist, std::type_identity>) {
			return TDistRet(r - l);
		} else {
			return dist_(l, r);
		}
	}

public:
    CoprimePairs(TSum sum_ = {}, TDist dist_ = {}):
        sum_(std::move(sum_)),
        dist_(std::move(dist_))
    {}

    TRet operator()(T lim) {
		if (lim == 1) return TRet(0);
		if (cache.count(lim) == 0) {
			TRet res(sum(lim));
			auto dlim = IntSqrt(lim);
			for (T i = 2; lim / i >= dlim; ++i) {
				auto mul = dist(i, i + 1);
				if (mul == 0) continue;
				res -= (*this)(lim / i) * mul;
			}
			for (T d = 2; d < dlim; ++d) {
				auto mul = dist(lim / (d + 1) + 1, lim / d + 1);
				if (mul == 0) continue;
				auto red = (*this)(d);
				red *= mul;
				res -= red;
			}
			cache[lim] = res;
		}
		return cache[lim];
    }
};