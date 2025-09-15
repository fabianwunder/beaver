#ifndef BEAVER_CONFIG_HPP
#define BEAVER_CONFIG_HPP

#if __cplusplus >= 201703L
  #define BEAVER_NODISCARD [[nodiscard]]
#else
  #define BEAVER_NODISCARD
#endif

#if defined(NDEBUG)
  #define BEAVER_ASSUME(cond) ((void)0)
#else
  #include <cassert>
  #define BEAVER_ASSUME(cond) assert(cond)
#endif

#endif // BEAVER_CONFIG_HPP
