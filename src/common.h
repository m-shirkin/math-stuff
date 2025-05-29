#pragma once

#include<utility>
#include<iostream>
#include<array>
#include<vector>
#include<tuple>
#include<cmath>
#include<optional>
#include<random>
#include<unordered_map>
#include<tuple>
#include<bitset>

#ifdef assert_
    #undef assert_
#endif
#define FWD(...) __VA_ARGS__
#if CHECKS
    bool __assert_triggered = false;
    #define assert_(...) if (!(__VA_ARGS__)) { \
        if (std::is_constant_evaluated()) { \
             __assert_triggered = true; \
        } else { \
            std::cerr << "assert_(" << #__VA_ARGS__ << ")\nFile " << __FILE__ << "(" << __LINE__ << ")\n"; \
            std::terminate(); \
        } \
    }
    #define checkval(val, ...) ([&]<typename __TCheckVal>(__TCheckVal&& _) -> decltype(auto) {__VA_ARGS__; return static_cast<__TCheckVal>(_);}(val))

    #define IF_CHECKS(...) __VA_ARGS__
    #define IF_NO_CHECKS(...)
#else
    #define assert_(...)
    #define checkval(val, ...) val
    #define IF_CHECKS(...)
    #define IF_NO_CHECKS(...) __VA_ARGS__
#endif

template<typename TA, typename TB>
constexpr auto min(TA&& a, TB&& b) {
    return a < b ? std::forward<TA>(a) : std::forward<TB>(b);
}

template<typename TA, typename TB>
constexpr auto max(TA&& a, TB&& b) {
	return a < b ? std::forward<TA>(b) : std::forward<TB>(a);
}

template<auto N, typename Arg0, typename... Args>
constexpr decltype(auto) select_nth(Arg0&& arg0, Args&&... args) noexcept {
    if constexpr(N == 0) {
        return std::forward<Arg0>(arg0);
    } else {
        return select_nth<N-1>(std::forward<Args>(args)...);
    }
}

template<typename T, typename TSub>
using SubVoid = std::conditional_t<std::is_void_v<T>, TSub, T>;

namespace detail {
	template<size_t N>
	struct uintN_Impl {};

	template<size_t N>
    struct intN_Impl {};

    template<size_t N> requires (N <= 64)
	struct uintN_Impl<N> {
		using type =
			std::conditional_t<N <= 8, uint8_t,
			std::conditional_t<N <= 16, uint16_t,
			std::conditional_t<N <= 32, uint32_t,
			uint64_t
		>>>;
	};

    template<size_t N> requires (N <= 64)
    struct intN_Impl<N> {
        using type = 
            std::conditional_t<N <= 8, int8_t,
            std::conditional_t<N <= 16, int16_t,
            std::conditional_t<N <= 32, int32_t,
            int64_t
        >>>;
    };
}

template<size_t N> using uintN = typename detail::uintN_Impl<N>::type;
template<size_t N> using intN = typename detail::intN_Impl<N>::type;

namespace detail {
    template<typename T> struct signed_unsigned_impl {};

    template<std::integral T>
    struct signed_unsigned_impl<T> {
        using unsigned_t = uintN<sizeof(T)*8>;
        using signed_t = intN<sizeof(T)*8>;
    };

    template<typename T>
    struct bit_size_impl {
        static constexpr auto value = sizeof(T)*8;
    };
}

template<typename T>
constexpr auto bit_size = detail::bit_size_impl<T>::value;

template<typename T> using to_unsigned = typename detail::signed_unsigned_impl<T>::unsigned_t;
template<typename T> using to_signed = typename detail::signed_unsigned_impl<T>::signed_t;

template<typename T> concept signed_int = std::is_same_v<T, to_signed<T>>;
template<typename T> concept unsigned_int = std::is_same_v<T, to_unsigned<T>>;

namespace overload_detail {
	constexpr int MAX_OVERLOAD_COUNT = 2;
}

#define OVERLOAD_DECL(tmpl, extmpl, Name,...) \
namespace overload_detail { \
    template<int, typename> \
    struct Name##_impl {}; \
\
    template<tmpl, int __prio, typename __TToken, typename... __Args> \
    static constexpr auto Name##_base(__Args&&... args) { \
        if constexpr(__prio < 0) { \
            static_assert(false, "No suitable overload found"); \
        } else if constexpr(requires {Name##_impl<__prio, __TToken>::template Call<extmpl>(std::forward<__Args>(args)...);}) { \
            return Name##_impl<__prio, __TToken>::template Call<extmpl>(std::forward<__Args>(args)...); \
        } else { \
            return Name##_base<extmpl, __prio - 1, __TToken>(std::forward<__Args>(args)...); \
        } \
    } \
\
    template<typename __TToken> \
    struct Name##_call { \
        template<tmpl, typename... __Args> \
        static constexpr auto Call(__Args&&... args) { \
            return Name##_base<extmpl, MAX_OVERLOAD_COUNT - 1, __TToken>(std::forward<__Args>(args)...); \
        } \
    }; \
}

#define OVERLOAD_CALL(Name) overload_detail::Name##_call<decltype([](){})>::template Call

#define OVERLOAD_IMPL(tmpl, Name, prio,...) \
namespace overload_detail { \
    static_assert(prio < MAX_OVERLOAD_COUNT, "Too many overloads!"); \
    template<typename __TToken> \
    struct Name##_impl<prio, __TToken> { \
        template<tmpl> \
        static constexpr auto Call __VA_ARGS__ \
    };\
};

OVERLOAD_DECL(typename T, T, DivMod)
#define DivMod OVERLOAD_CALL(DivMod)

OVERLOAD_IMPL(typename T, DivMod, 0, (T const& a, T const& b, T& quo, T& rem) {
    quo = a / b;
    rem = a % b;
})

#define to_constexpr(...) []() consteval { return __VA_ARGS__; }()

auto constexpr ArrayFromVector(auto fn) {
    auto vec = fn();
    std::array<std::remove_cvref_t<decltype(vec[0])>, fn().size()> res;
    for (size_t i = 0; i < vec.size(); ++i) {
        res[i] = vec[i];
    }
    return res;
}

constexpr auto NPOS = size_t(-1);

template<typename T>
auto& to_lvalue(T&& val) {
    return static_cast<T&>(val);
}

template<typename T>
constexpr decltype(auto) Abs(T const& x) {
    if constexpr(unsigned_int<T>) return x; else {
		return to_unsigned<T>(x < 0 ? -x : x);
    }
}

namespace detail {
    template<template<typename> typename TMPL, typename T>
    struct IsTemplateOf;

    template<template<typename> typename TMPL, typename T>
    struct IsTemplateOf<TMPL, TMPL<T>> {};
}

template<typename T, template<typename> typename TMPL>
concept template_of = requires { sizeof(detail::IsTemplateOf<TMPL, typename std::remove_cvref_t<T>>); };
