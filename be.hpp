#ifndef BE_HPP_
#define BE_HPP_

#include <numeric>      // accumulate

#include "libfa/fa.h"
#include "types.hpp"
#include "bitset.hpp"
#include "util.hpp"
#include "order.hpp"
#include "io.hpp"

template <typename T>
size_t push_bucket(T const &f, vector<vector<T>> &buckets, vector<size_t> const &pos) {

        const size_t b = *max_element(f.vars.begin(), f.vars.end(), compare_pos(pos));
        buckets[b].push_back(f);
        return b;
}

template <typename T>
vector<vector<T>> compute_buckets(vector<T> const &functions, vector<size_t> const &pos) {

        vector<vector<T>> buckets(pos.size(), vector<T>());

        for (T const &f : functions) {
                push_bucket(f, buckets, pos);
        }

        return buckets;
}

__attribute__((always_inline)) inline
void free_bucket(vector<automata> &bucket) {

        for (auto& a : bucket) {
                for (auto& [ v, fa ] : a.rows) {
                        fa_free(fa);
                }
                a.rows.clear();
        }

        bucket.clear();
}

__attribute__((always_inline)) inline
automata copy_automata(automata const &a) {

        automata res = {
                a.vars,
                a.domains,
                unordered_map<value, struct fa *>()
        };

        for (auto& [ v, fa ] : a.rows) {
                res.rows.insert({ v, fa_clone(fa) });
        }

        return res;
}

automata join_bucket(vector<automata> const &bucket, vector<size_t> domains);

void reduce_var(automata &a, size_t var);

void bucket_elimination(vector<vector<automata>> &buckets, vector<size_t> const &order,
                        vector<size_t> const &pos, vector<size_t> const &domains,
                        size_t max_iter = numeric_limits<size_t>::max());

#endif /* BE_HPP_ */
