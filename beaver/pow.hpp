#ifndef BEAVER_POW_HPP
#define BEAVER_POW_HPP
#include <cmath>

namespace beaver {

// Tunable: beyond this |n| we defer to std::pow for floating bases
inline constexpr int POW_STDPOW_CUTOFF = 64;

//---------------------------------------------
// 1) Compile-time exponent: pow<N>(x)
//    Exponentiation by squaring (balanced tree).
//---------------------------------------------
template<int N, class T>
inline __attribute__((always_inline))
constexpr T pow(T x) noexcept {
    if constexpr (N == 0) {
        return T{1};
    } else if constexpr (N > 0) {
        if constexpr ((N & 1) == 0) {           // even positive
            T y = pow<N/2>(x);
            return y * y;
        } else {                                 // odd positive
            return x * pow<N-1>(x);
        }
    } else { // N < 0
        if constexpr ((N & 1) == 0) {           // even negative
            T y = pow<N/2>(x);                  // still negative, halved
            return T{1} / (y * y);
        } else {                                 // odd negative
            return (T{1} / x) * pow<N+1>(x);    // moves toward zero without -N
        }
    }
}


//---------------------------------------------
// 2) Integer-base exact power (wraparound):
//    ipow(T base, unsigned n)
//    Only for integral T; uses binary exponentiation.
//    WARNING: may overflow T (intentional).
//---------------------------------------------
template<class T>
inline __attribute__((always_inline))
std::enable_if_t<std::is_integral_v<T>, T>
ipow(T x, unsigned n) noexcept {
    T res = 1;
    while (n) {
        if (n & 1) res = static_cast<T>(res * x);
        x = static_cast<T>(x * x);
        n >>= 1;
    }
    return res;
}

//---------------------------------------------
// 3) Runtime exponent: pow(x, n)
//    - Floating bases: binary exp for small |n|, else std::pow.
//    - Integral bases: promote to double/long double and use floating route,
//      unless the caller explicitly chose ipow().
//---------------------------------------------
namespace detail {
    template<class T>
    inline T pow_small_int(T x, int n) noexcept {
        // n assumed small magnitude here
        if (n == 0) return T{1};
        bool neg = (n < 0);
        unsigned uu = static_cast<unsigned>(neg ? -static_cast<long long>(n) : n);
        T res = T{1};
        while (uu) {
            if (uu & 1u) res *= x;
            x *= x;
            uu >>= 1u;
        }
        return neg ? T{1} / res : res;
    }
}

template<class T>
inline T pow(T x, int n) noexcept {
    if constexpr (std::is_floating_point_v<T>) {
        // For small |n|, manual binary exponentiation beats std::pow
        if (n >= -POW_STDPOW_CUTOFF && n <= POW_STDPOW_CUTOFF) {
            return detail::pow_small_int(x, n);
        }
        // Large |n|: delegate to libm for overflow/underflow/denormal behavior
        if constexpr (std::is_same_v<T, long double>) {
            return std::pow(x, static_cast<long double>(n));
        } else {
            return std::pow(static_cast<double>(x), static_cast<double>(n));
        }
    } else if constexpr (std::is_integral_v<T>) {
        // “Stable against large integers”: by default, promote and use floating.
        // If you truly want integer wraparound, call ipow() explicitly.
        using F = long double; // widest reasonable libm type
        F xf = static_cast<F>(x);
        if (n >= -POW_STDPOW_CUTOFF && n <= POW_STDPOW_CUTOFF) {
            return static_cast<T>(detail::pow_small_int(xf, n));
        }
        return static_cast<T>(std::pow(xf, static_cast<F>(n)));
    } else {
        // Fallback: try to treat as floating-like
        using F = long double;
        return static_cast<T>(std::pow(static_cast<F>(x), static_cast<F>(n)));
    }
}

} // namespace beaver
#endif