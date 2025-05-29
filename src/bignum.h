#pragma once

#include"common.h"
#include<boost/multiprecision/cpp_int.hpp>

using iBig = boost::multiprecision::cpp_int;

namespace detail {
	template<size_t N> requires (N > 64)
	struct intN_Impl<N> {
        using type = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<0, ((N + 63) >> 6 << 6), boost::multiprecision::signed_magnitude, boost::multiprecision::IF_CHECKS(checked)IF_NO_CHECKS(unchecked), void>>;
    };

    template<size_t MinBits, size_t MaxBits, boost::multiprecision::cpp_integer_type SignType, boost::multiprecision::cpp_int_check_type Checked, typename Allocator>
    struct signed_unsigned_impl<boost::multiprecision::number<boost::multiprecision::cpp_int_backend<MinBits, MaxBits, SignType, Checked, Allocator>>> {
        using signed_t = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<MinBits, MaxBits, boost::multiprecision::signed_magnitude, Checked, Allocator>>;
        using unsigned_t = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<MinBits, MaxBits, boost::multiprecision::unsigned_magnitude, Checked, Allocator>>;
    };

    template<>
    struct bit_size_impl<iBig> {
        static constexpr auto value = size_t(-1);
    };

    template<size_t MinBits, size_t MaxBits, boost::multiprecision::cpp_integer_type SignType, boost::multiprecision::cpp_int_check_type Checked, typename Allocator>
    struct bit_size_impl<boost::multiprecision::number<boost::multiprecision::cpp_int_backend<MinBits, MaxBits, SignType, Checked, Allocator>>> {
        static constexpr auto value = MaxBits;
    };

	template<size_t N> requires (N > 64)
	struct uintN_Impl<N> {
		using type = to_unsigned<intN<N>>;
	};
}

OVERLOAD_IMPL(typename T, DivMod, 1, (T const& a, T const& b, T& quo, T& rem) {
	boost::multiprecision::divide_qr(a, b, quo, rem);
})

namespace detail {
    template<char... Ch>
    struct literal_impl {
        static consteval auto convert() {
            intN<sizeof...(Ch) * 392 / 59> res = 0;
            ([&](){res *= 10; res += Ch - '0';}(),...);
            return res;
        }
    };

    template<char... Ch>
    struct literal_impl<'0', Ch...> {
        static consteval auto convert() {
            intN<sizeof...(Ch) * 6> res = 0;
            ([&]() {res *= 8; res += Ch - '0';}(), ...);
            return res;
        }
    };

    template<char C0, char... Ch> requires (C0 == 'b' || C0 == 'B')
        struct literal_impl<'0', C0, Ch...> {
        static consteval auto convert() {
            intN<sizeof...(Ch) * 2> res = 0;
            ([&]() {res *= 2; res += Ch - '0';}(), ...);
            return res;
        }
    };

    template<char C0, char... Ch> requires (C0 == 'x' || C0 == 'X')
        struct literal_impl<'0', C0, Ch...> {
        static consteval auto convert() {
            intN<sizeof...(Ch) * 8> res = 0;
            ([&]() {res *= 16; res += Ch >= 'a' ? Ch - 'a' + 10 : Ch - '0';}(), ...);
            return res;
        }
    };
}

template<char... Ch>
consteval auto operator""_n() {
    return detail::literal_impl<Ch...>::convert();
}