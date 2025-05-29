#pragma once

#include "common.h"

template<auto r, typename T = long long>
struct IntPlusRoot {
	T a, b;

	template<typename TA, typename TB>
	constexpr IntPlusRoot(TA&& a_, TB&& b_):
		a(std::forward<TA>(a_)), b(std::forward<TB>(b_))
	{}

	template<typename TA>
	constexpr IntPlusRoot(TA&& a_) :
		IntPlusRoot(std::forward<TA>(a_), 0)
	{}
	
	template<typename TO>
	constexpr IntPlusRoot& operator=(TO&& other) {
		if constexpr(std::is_same_v<typename std::remove_cvref_t<TO>, IntPlusRoot>) {
			a = std::forward<TO>(other).a;
			b = std::forward<TO>(other).b;
		} else {
			a = std::forward<TO>(other);
			b = 0;
		}
		return *this;
	}

#define DEFINE_BROADCAST_OP(op) \
	constexpr IntPlusRoot& operator op(IntPlusRoot const& other) { \
		a op other.a; \
		b op other.b; \
		return *this; \
	}

	DEFINE_BROADCAST_OP(+=)
	DEFINE_BROADCAST_OP(-=)
#undef DEFINE_BROADCAST_OP

	constexpr IntPlusRoot& operator*=(IntPlusRoot const& other) {
		std::tie(a, b) = std::make_tuple(a * other.a + b * other.b * r, a * other.b + b * other.a);
		return *this;
	}

	constexpr IntPlusRoot& operator/=(IntPlusRoot const& other) {
		assert_(!unsigned_int<T>);
		*this *= IntPlusRoot(other.a, -other.b);
		return *this /= other.a * other.a - other.b * other.b * r;
	}

#define DEFINE_INT_OP(op, ...) \
	template<typename TO> \
	constexpr IntPlusRoot& operator op(TO const& other) { \
		a op other; \
		__VA_ARGS__; \
		return *this; \
	}

	DEFINE_INT_OP(+=)
	DEFINE_INT_OP(-=)
	DEFINE_INT_OP(*=, b *= other)
	DEFINE_INT_OP(/=, b /= other)
#undef DEFINE_INT_OP
};

template<auto n, auto m, typename T>
struct StaticMatrix {
	std::array<T, n*m> vals;

	constexpr T& el(auto x, auto y) { return vals[x * m + y]; }
	constexpr T const& el(auto x, auto y) const& { return vals[x * m + y]; }

	template<typename T2>
	constexpr StaticMatrix& operator=(T2&& other) {
		if constexpr(std::is_same_v<typename std::remove_cvref_t<T2>, StaticMatrix>) {
			vals = std::forward<T2>(other).vals;
		} else {
			static_assert(n == m);
			for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < m; ++j) {
				el(i, j) = (i == j) * other;
			}
		}
		return *this;
	}
	
	template<typename T2>
	explicit constexpr StaticMatrix(T2&& other) {
		*this = std::forward<T2>(other);
	}

	constexpr StaticMatrix():
		vals{}
	{}

#define DEFINE_BROADCAST_OP(op) \
	constexpr StaticMatrix& operator op(StaticMatrix const& other) { \
		for (size_t i = 0; i < vals.size(); ++i) { \
			el(i) op other.el(i); \
		} \
		return *this; \
	}

	DEFINE_BROADCAST_OP(+=)
	DEFINE_BROADCAST_OP(-=)
#undef DEFINE_BROADCAST_OP

	constexpr StaticMatrix& operator*=(StaticMatrix const& other) {
		static_assert(n == m);
		StaticMatrix res(0);
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < n; ++j) {
				for (size_t k = 0; k < n; ++k) {
					res.el(i, j) += el(i, k) * other.el(k, j);
				}
			}
		}
		return *this = res;
	}

	constexpr auto Apply(auto const& vec) const& {
		std::array<T, m> res = {};
		for (size_t i = 0; i < m; ++i) {
			for (size_t j = 0; j < n; ++j) {
				res[j] += vec[i] * el(i, j);
			}
		}
		return res;
	}
};

template<typename TSeq>
constexpr size_t SequenceLength = []<template <typename T, T... Ns> typename TSeqT, typename T, T... Ns>(std::type_identity<TSeqT<T, Ns...>>){
	return sizeof...(Ns);
}(std::type_identity<TSeq>());

template<typename TSeq>
using InverseSequence = decltype(
	[]<template <typename T, T... Ns> typename TSeqT, typename T, T... Ns, size_t... is>
	(std::index_sequence<is...>, std::type_identity<TSeqT<T, Ns...>>) {
		return std::type_identity<TSeqT<T, select_nth<sizeof...(Ns) - is - 1>(Ns...)...>>();
	}(std::make_index_sequence<SequenceLength<TSeq>>(), std::type_identity<TSeq>())
)::type;

template<typename TVal, int ord, typename TCont = std::vector<TVal>>
struct Tensor {
	template<typename TVal_, int ord_, typename TCont_>
	friend struct Tensor;
private:
	template<typename TNCont>
	Tensor(std::array<size_t, ord> dim, std::array<size_t, ord> stride, TNCont&& data):
		dim(dim),
		stride(stride),
		data(std::forward<TNCont>(data))
	{}

	template<size_t... is>
	constexpr TVal& at_impl(std::index_sequence<is...>, auto... inds) {
		return data[((stride[is] * inds) + ...)];
	}

	template<size_t... axis, typename TNCont>
	constexpr auto T_impl(TNCont&& ndata) {
		return Tensor<TVal, ord, TNCont>(
			std::array<size_t, ord>{dim[axis]...},
			std::array<size_t, ord>{stride[axis]...},
			std::forward<TNCont>(ndata)
		);
	}

	template<size_t... axis_>
	constexpr void ForEach_impl(std::index_sequence<axis_...>, auto const& func) {
		[&]<size_t... is, size_t... axis>(std::index_sequence<is...>, std::index_sequence<axis...>) {
			for (std::array<size_t, sizeof...(is)> i{};;) {
				func(i[sizeof...(is) - is - 1]...);
				if (([&]() {
					if (i[is] == dim[axis] - 1) {
						i[is] = 0;
						return true;
					} else {
						++i[is];
						return false;
					}
				}() && ... && true)) break;
			}
		}(InverseSequence<std::make_index_sequence<sizeof...(axis_)>>(), InverseSequence<std::index_sequence<axis_...>>());
	}

public:
	// TODO: initial value
	template<typename... Ts> requires (sizeof...(Ts) == ord)
	constexpr Tensor(Ts... dims) :
		dim{ size_t(dims)... },
		data((... * dims))
	{
		assert_(((dims > 0) && ...));
		size_t f = data.size();
		auto it = stride.begin();
		([&]() {
			f /= dims;
			*it++ = f;
		}(), ...);
	}

	constexpr TVal& at(auto... inds) {
		static_assert(sizeof...(inds) == ord);
		return at_impl(std::make_index_sequence<ord>(), inds...);
	}

	constexpr TVal const& at(auto... inds) const& {
		return const_cast<Tensor&>(*this).at(inds...);
	}

	template<size_t... axis>
	constexpr void ForEach(auto const& func) {
		if constexpr(sizeof...(axis) == 0) {
			ForEach_impl(std::make_index_sequence<ord>(), func);
		} else {
			ForEach_impl(std::index_sequence<axis...>(), func);
		}
	}

	template<size_t... axis>
	constexpr auto T() & {
		return T_impl<axis...>(data);
	}

	template<size_t... axis>
	constexpr auto T() && {
		return T_impl<axis...>(std::move(data));
	}

	std::array<size_t, ord> dim;
	std::array<size_t, ord> stride;
	TCont data;
};