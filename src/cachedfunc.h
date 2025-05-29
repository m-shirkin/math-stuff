#pragma once

#include "common.h"

namespace detail {
	template<typename TVal, typename... Ts>
	struct TupleMap;

	template<typename TVal, typename T>
	struct TupleMap<TVal, T> : std::unordered_map<T, TVal> {
		auto Access(T t) {
			return std::unordered_map<T, TVal>::try_emplace(t);
		}
	};

	template<typename TVal, typename T0, typename... Ts>
	struct TupleMap<TVal, T0, Ts...> : std::unordered_map<T0, TupleMap<TVal, Ts...>> {
		auto Access(T0 t0, Ts... ts) {
			auto itb = std::unordered_map<T0, TupleMap<TVal, Ts...>>::try_emplace(t0);
			return itb.first->second.Access(ts...);
		}
	};

	template<typename Fn, typename Fn2>
	struct CachedFunc;

	template<typename Fn2, typename TRet, typename... Args>
	struct CachedFunc<TRet(Args...), Fn2> {
		CachedFunc(Fn2 fn) :
			fn(std::move(fn))
		{}

		template<typename... NArgs>
		auto operator()(Args... args, NArgs&&... nargs) const& {
			auto [it, b] = cache.Access(args...);
			if (b) {
				it->second = fn(args..., std::forward<NArgs>(nargs)...);
			}
			return it->second;
		}

		Fn2 fn;
		mutable TupleMap<TRet, Args...> cache;
	};
}

template<typename Fn, typename Fn2>
auto MakeCached(Fn2 fn) {
	return detail::CachedFunc<Fn, Fn2>(std::move(fn));
}