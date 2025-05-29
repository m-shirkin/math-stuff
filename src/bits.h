#pragma once

#include "common.h"

#include<bit>
#include<cstdint>

template<typename T>
constexpr size_t BitLength(T const& val) {
    if constexpr(requires {msb(val);}) {
        return msb(val) + 1;
    } else {
        return sizeof(val)*8 - std::countl_zero(static_cast<to_unsigned<T>>(val));
    }
}

template<typename T>
constexpr bool BitVal(T const& val, size_t i) {
    if constexpr(requires {bit_test(val, i);}) {
        return bit_test(val, i);
    } else {
        using UT = to_unsigned<T>;
        return static_cast<bool>(static_cast<UT>(val) & (static_cast<UT>(1) << i));
    }
}

template<typename T>
constexpr size_t EndZeroes(T const& val) {
    assert_(val != 0);
	if constexpr (requires {lsb(val); }) {
        return lsb(val);
    } else {
        return std::countr_zero(val);
    }
}