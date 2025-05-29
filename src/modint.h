#pragma once

#include"common.h"
#include"divisors.h"
#include"intfunc.h"

template<typename T>
struct ModInt;
    
template<typename T>
concept mod_int = template_of<T, ModInt>;

template<unsigned_int T>
struct ModInt<T> {
    T const mod;
    T val;

    static auto& GlobalMod() {
        static std::optional<T> g_mod = std::nullopt;
        return g_mod;
    }

    struct SetGlobalMod {
        SetGlobalMod(T mod) {GlobalMod() = std::move(mod);}
        ~SetGlobalMod() {GlobalMod() = std::nullopt;}
    };

    template<typename TVal = T>
    constexpr ModInt(T mod_, TVal&& val_):
        mod(std::move(mod_)),
        val([&]() -> T {
            if constexpr (signed_int<TVal>) {
                if (val_ < 0) {
                    return mod - 1 - (T(-val_) - 1) % mod;
                }
            }
            if (val_ < mod) {
                return std::forward<TVal>(val_);
            } else {
                return std::forward<TVal>(val_) % mod;
            }
        }())
    {
    }

    template<mod_int T_>
	constexpr ModInt(T_&& x) :
        ModInt(std::forward<T_>(x).mod, std::forward<T_>(x).val)
    {}

    template<typename T_>
    explicit ModInt(T_&& x):
        ModInt(GlobalMod() ? *GlobalMod() : std::forward<T_>(x), GlobalMod() ? std::forward<T_>(x) : 0)
    {}

    ModInt():
        ModInt(*GlobalMod(), 0)
    {}

    constexpr ModInt& operator=(ModInt const& other) {
		assert_(mod == other.mod);
		val = other.val;
        return *this;
    }

    template<typename T_>
    constexpr ModInt& operator=(T_&& other) {
        if constexpr(mod_int<T_>) {
			assert_(mod == other.mod);
            val = std::forward<T_>(other).val;
        } else {
			val = other % mod;
        }
        return *this;
    }

    template<typename TO>
    constexpr ModInt& operator+=(TO const& in) {
        if constexpr(mod_int<TO>) {
			assert_(mod == in.mod);
			val += in.val;
			if (val >= mod) val -= mod;
			return *this;
        } else {
            if constexpr(signed_int<TO>) if (in < 0) return *this -= -in;
			val = T((val + in) % mod);
			return *this;
        }
    }

    template<typename TO>
    constexpr ModInt& operator-=(TO const& in) {
		if constexpr (mod_int<TO>) {
			assert_(mod == in.mod);
			if (in.val > val) {
				val += mod;
				val -= in.val;
				if (val >= mod) val -= mod;
			} else {
				val -= in.val;
			}
			return *this;
		} else {
            if constexpr (signed_int<TO>) if (in < 0) return *this += -in;
			val = T((mod - in % mod + val) % mod);
			return *this;
        }
    }

    template<typename TO>
    constexpr ModInt& operator*=(TO const& in) {
        assert_(BitLength(mod) * 2 <= bit_size<T>);
		if constexpr (mod_int<TO>) {
			assert_(mod == in.mod);
			val *= in.val;
			val %= mod;
			return *this;
        } else {
            if constexpr(signed_int<TO>) if (in < 0) {
                val = T(-in % mod) * (mod - val) % mod;
                return *this;
            }
			val *= in % mod;
            val %= mod;
			return *this;
        }
    }

    template<typename TO>
    constexpr ModInt& operator/=(TO const& in) {
        assert_(BitLength(mod) * 2 <= bit_size<T>);
		if constexpr (mod_int<TO>) {
			assert_(mod == in.mod);
			assert_(in.val != 0);
			val *= ModularInverse<T>(T(in.val), mod);
			val %= mod;
			return *this;
        } else {
            assert_(in != 0);
            if constexpr (signed_int<TO>) if (in < 0) {
				val = ModularInverse<T>(T(-in % mod), mod) * (mod - val) % mod;
				return *this;
			}
			val *= ModularInverse<T>(T(in % mod), mod);
			val %= mod;
			return *this;
        }
    }

    constexpr ModInt& operator++() {
        return *this += 1;
    }

    constexpr ModInt operator++(int) const& {
        ModInt res = *this;
        return res += 1;
    }

    constexpr ModInt& operator--() {
        return *this -= 1;
    }

    constexpr ModInt operator--(int) const& {
        ModInt res = *this;
        return res -= 1;
    }

    constexpr ModInt operator-() const& {
        return ModInt(mod, val == 0 ? 0 : mod - val);
    }
};

template<typename T, typename T2>
ModInt(T, T2) -> ModInt<to_unsigned<T>>;

#define DEFINE_MODINT_OP(op) \
    template<typename T1, typename T2> requires (mod_int<T1> || mod_int<T2>) \
    constexpr auto operator op(T1 const& t1, T2 const& t2) { \
        if constexpr(mod_int<T1>) { \
            T1 t1c = t1; \
            t1c op##= t2; \
            return t1c; \
        } else { \
            T2 t1c(t2.mod, t1); \
            t1c op##= t2; \
            return t1c; \
        } \
    }

DEFINE_MODINT_OP(+)
DEFINE_MODINT_OP(-)
DEFINE_MODINT_OP(*)
DEFINE_MODINT_OP(/)

#undef DEFINE_MODINT_OP

template<typename OStrm, typename T>
auto& operator<<(OStrm& strm, ModInt<T> const& m) {
	return strm << m.val << " (mod " << m.mod << ")";
}

template<mod_int T1, mod_int T2>
constexpr bool operator==(T1 const& lhs, T2 const& rhs) {
    return lhs.mod == rhs.mod && lhs.val == rhs.val;
}

// Mod should be prime
template<typename T, typename TV>
std::optional<T> constexpr ModSqrt(T const& mod, TV const& val_) {
    assert_(BitLength(mod) * 2 <= bit_size<T>);
    if (auto val = val_ % mod; val == 1 || val == 0) return T(val);
    if (mod <= 3) return std::nullopt;
    ModInt<T> z(mod, 2);
    for(;;) {
        auto res = Power(z, (mod - 1) >> 1).val;
        if (res != 1) {
            assert_(res == mod-1);
            break;
        }
        ++z;
    } 
    size_t s = 0;
    T q = mod - 1;
    while (!(q & 1)) {
        ++s;
        q >>= 1;
    }
    z = Power(z, q);
    T val(val_ > 0 ? val_ % mod : mod - 1 - (-val_ - 1) % mod);
    if (Power(ModInt<T>(mod, val), (mod - 1) >> 1).val != 1) return std::nullopt;
    auto R = Power(ModInt<T>(mod, val), (q - 1) >> 1);
    auto t = R;
    t *= R;
    t *= val;
    R *= val;
    while (t.val != 1) {
        assert_(s > 1);
        --s;
        auto tsq = t;
        for (size_t m = 0; m < s - 1; ++m) {
            tsq *= tsq;
        }
        if (tsq.val == mod - 1) {
            t *= z;
            t *= z;
            R *= z;
        } else {
            assert_(tsq.val == 1);
        }
        z *= z;
    }
    assert_(R.val * R.val % mod == val);
    return R.val;
}