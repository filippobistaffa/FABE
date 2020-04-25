#include "conversion.hpp"

#include <cmath>        // fabs
#include <sstream>      // ostringstream
#include <iostream>     // cout
#include <cstring>      // strlen
#include <numeric>      // accumulate
#include <algorithm>    // unique

#include "libfa/fa.h"
#include "kmeans.hpp"
//#include "Ckmeans.1d.dp.h"
#include "io.hpp"

//#define OLD

#ifdef OLD

static inline struct fa *fa_compile_minimise(vector<pair<vector<size_t>, value>>::const_iterator begin,
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
}

#else

static inline struct fa *fa_compile_minimise(vector<pair<vector<size_t>, value>>::const_iterator begin,
                                             vector<pair<vector<size_t>, value>>::const_iterator end) {

        struct fa *fa = NULL;

        for (auto row = begin; row != end; ++row) {
                ostringstream oss;
                for (auto i : row->first) {
                        oss << ALPHABET[i];
                }
                if (fa) {
                        fa_add_word(fa, oss.str().c_str(), oss.str().length());
                } else {
                        fa_compile(oss.str().c_str(), oss.str().length(), &fa);
                }
        }

        fa_minimize(fa);
        return fa;
}

#endif

template <typename T, typename PRED>
static inline size_t count_uniques(T first, T last, PRED equal) {

        if (first == last) {
                return 0;
        } else {
                size_t n = 1;
                for (T itr = first + 1; itr != last; ++itr) {
                        if (!equal(*itr, *(itr - 1))) {
                                n++;
                        }
                }
                return n;
        }
}

pair<automata, value> compute_automata(table const &t, size_t max_k) {

        automata res = {
                vector<size_t>(t.vars),
                vector<size_t>(t.domains),
        };

        vector<double> values(t.rows.size());

        for (size_t i = 0; i < t.rows.size(); ++i) {
                values[i] = t.rows[i].second;
        }

        size_t uniques = count_uniques(values.begin(), values.end(), 
                                       [](double const &x, double const &y)
                                       { return (fabs(x - y) <= std::numeric_limits<double>::epsilon()); });

        #ifdef DEBUG
        cout << endl << vec2str(values, "values") << endl;
        cout << "k = " << min(uniques, max_k) << endl;
        #endif
        size_t k = min(uniques, max_k);
        value max_error = 0;

        if (uniques > 1) {
                auto [ clusters, centers ] = kmeans_1d_dp(values, k);
                #ifdef DEBUG
                cout << vec2str(clusters, "clusters") << endl;
                cout << vec2str(centers, "centers") << endl;
                vector<value> errors(values.size());
                for (size_t i = 0; i < values.size(); ++i) {
                        errors[i] = fabs(centers[clusters[i]] - values[i]);
                }
                cout << vec2str(errors, "errors") << endl;
                #endif
                size_t begin = 0;
                size_t i = begin;
                while (begin != t.rows.size()) {
                        while (i != t.rows.size() && clusters[begin] == clusters[i]) {
                                max_error = max((value)fabs(centers[clusters[i]] - values[i]), max_error);
                                i++;
                        }
                        res.rows.insert({ centers[clusters[begin]],
                                          fa_compile_minimise(t.rows.begin() + begin, t.rows.begin() + i) });
                        begin = i;
                }
        } else {
                res.rows.insert({ t.rows[0].second, fa_compile_minimise(t.rows.begin(), t.rows.end()) });
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
                                row[v] = string(ALPHABET).find(rows[r][v]);
                        }
                        res.rows.push_back(make_pair(row, v));
                        free(rows[r]);
                }
                free(rows);
        }

        return res;
}
