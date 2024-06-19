#ifndef UTIL_HPP_
#define UTIL_HPP_

inline std::vector<size_t> get_combination(size_t i, std::vector<size_t> const &v) {

    std::vector<size_t> res(v.size());
    size_t x = i;
    for (size_t j = 0; j < v.size(); ++j) {
        res[j] = x % v[j];
        x /= v[j];
    }
    return res;
}

template <typename T1, typename T2, typename T3>
inline void exclusive_scan(T1 begin, T1 end, T1 out, T2 init, T3 op) {

    *out = init;
    for (auto it = begin; it != prev(end); ++it, ++out) {
        *(out + 1) = op(*it, *out);
    }
}

#define BREAKPOINT(MSG) do { puts(MSG); fflush(stdout); while (getchar() != '\n'); } while (0)

#endif
