#include "conversion.hpp"

#include <cmath>        // fabs
#include <sstream>      // ostringstream
#include <iostream>     // cout
#include <cstring>      // strlen
#include <numeric>      // accumulate

#include "libfa/fa.h"

//#define OLD

#ifdef OLD

/*static struct fa *fa_compile_minimise(vector<pair<vector<size_t>, value>>::const_iterator begin,
                                      vector<pair<vector<size_t>, value>>::const_iterator end) {

        ostringstream oss;

        for (auto row = begin; row != end; ++row) {
                for (auto i : row->first) {
                        oss << ALPHABET[i];
                }
                oss << "|";
        }

        struct fa *fa;
        fa_compile(oss.str().c_str(), oss.str().length() - 1, &fa);
        fa_minimize(fa);
        return fa;
}*/

#else

static struct fa *fa_compile_minimise(vector<pair<vector<size_t>, value>>::const_iterator begin,
                                      vector<pair<vector<size_t>, value>>::const_iterator end) {

        struct fa *fa = fa_make_basic(FA_EMPTY);

        for (auto row = begin; row != end; ++row) {
                fa_add_word(fa, &(row->first[0]), row->first.size());
        }

        fa_minimize(fa);
        return fa;
}

#endif

__attribute__((always_inline)) inline
bool are_equal(value a, value b, value tolerance) {

        const value epsilon = numeric_limits<value>::epsilon();
        return fabs(a - b) <= tolerance + epsilon;
};

pair<automata, value> compute_automata(table const &t, value tolerance) {

        automata res = {
                vector<size_t>(t.vars),
                vector<size_t>(t.domains),
        };

        auto begin = t.rows.begin();
        auto end = begin;
        value max_error = 0;

        while (begin != t.rows.end()) {
                while (end != t.rows.end() && are_equal(end->second, begin->second, tolerance)) {
                        end++;
                }
                max_error = max(max_error, prev(end)->second - begin->second);
                res.rows.insert({ prev(end)->second, fa_compile_minimise(begin, end) });
                begin = end;
        }

        return make_pair(res, max_error);
}

table compute_table(automata const &a) {

        table res = {
                vector<size_t>(a.vars),
                vector<size_t>(a.domains),
        };

        const size_t max_rows = accumulate(a.domains.begin(), a.domains.end(), 1, multiplies<size_t>());

        for (auto const &[ v, fa ] : a.rows) {
                char **rows;
                const size_t n_rows = fa_enumerate(fa, max_rows, &rows);
                for (size_t r = 0; r < n_rows; ++r) {
                        vector<size_t> row = vector<size_t>(res.vars.size());
                        for (size_t v = 0; v < row.size(); ++v) {
                                row[v] = rows[r][v];
                        }
                        res.rows.push_back(make_pair(row, v));
                        free(rows[r]);
                }
                free(rows);
        }

        return res;
}
