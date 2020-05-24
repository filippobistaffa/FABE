#include "conversion.hpp"

#include <cmath>        // fabs
#include <sstream>      // ostringstream
#include <iostream>     // cout
#include <cstring>      // strlen
#include <numeric>      // accumulate

#include "libfa/fa.h"
#include "io.hpp"

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

//#define DEBUG_CLUSTERING

__attribute__((always_inline)) inline
bool are_equal(value a, value b, value tolerance) {

        const value epsilon = numeric_limits<value>::epsilon();
        #ifdef DEBUG_CLUSTERING
        cout << a << " == " << b << "? ";
        cout << (fabs(a - b) <= tolerance + epsilon ? "\033[1;32mYES\033[0m" : "\033[1;31mNO\033[0m") << endl;
        #endif
        return fabs(a - b) <= tolerance + epsilon;
};

pair<automata, value> compute_automata(table const &t, value tolerance, int type) {

        automata res = {
                vector<size_t>(t.vars),
                vector<size_t>(t.domains),
        };

        auto begin = t.rows.begin();
        auto end = begin;
        value max_error = 0;
        auto conv = [type](value val) { return type == WCSP ? val : exp(-val); };

        while (begin != t.rows.end()) {
                while (end != t.rows.end() && are_equal(conv(end->second), conv(begin->second), tolerance)) {
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

        preallocate_rows(res, numeric_limits<value>::max());
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
