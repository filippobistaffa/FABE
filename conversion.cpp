#include "conversion.hpp"

#include <cmath>    // fabs
#include <sstream>  // ostringstream
#include <iostream> // cout
#include <cstring>  // strlen
#include <numeric>  // accumulate

#include "libfa/fa.h"
#include "io.hpp"

static struct fa *fa_compile_minimise(std::vector<std::pair<std::vector<size_t>, value>>::const_iterator begin,
                                      std::vector<std::pair<std::vector<size_t>, value>>::const_iterator end) {

    struct fa *fa = fa_make_basic(FA_EMPTY);

    for (auto row = begin; row != end; ++row) {
        fa_add_word(fa, &(row->first[0]), row->first.size());
    }

    fa_minimize(fa);
    return fa;
}

__attribute__((always_inline)) inline
bool are_equal(value a, value b, value tolerance) {

    const value epsilon = std::numeric_limits<value>::epsilon();
    return fabs(a - b) <= tolerance + epsilon;
}

std::pair<automata, value> compute_automata(table const &t, value tolerance) {

    automata res = {
        std::vector<size_t>(t.vars),
        std::vector<size_t>(t.domains),
    };

    auto begin = t.rows.begin();
    auto end = begin;
    value max_error = 0;

    while (begin != t.rows.end()) {
        while (end != t.rows.end() && are_equal(end->second, begin->second, tolerance)) {
            end++;
        }
        max_error = std::max(max_error, prev(end)->second - begin->second);
        res.rows.insert({ prev(end)->second, fa_compile_minimise(begin, end) });
        begin = end;
    }

    return std::make_pair(res, max_error);
}

table compute_table(automata const &a) {

    table res = {
        std::vector<size_t>(a.vars),
        std::vector<size_t>(a.domains),
    };

    preallocate_rows(res, std::numeric_limits<value>::max());

    size_t *idx = new size_t[res.rows.size()];

    for (auto const &[ v, fa ] : a.rows) {
        auto n = fa_enumerate_idx(fa, &a.domains[0], idx);
        for (size_t i = 0; i < n; ++i) {
            res.rows[idx[i]].second = v;
        }
    }

    delete[] idx;

    return res;
}
