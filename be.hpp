#ifndef BE_HPP_
#define BE_HPP_

#include <numeric>      // accumulate
#include "libfa/fa.h"
#include "types.hpp"
#include "order.hpp"

__attribute__((always_inline)) inline
size_t push_bucket(automata const &a, vector<vector<automata>> &buckets, vector<size_t> const &pos) {

        const size_t b = *max_element(a.vars.begin(), a.vars.end(), compare_pos(pos));
        buckets[b].push_back(move(a));
        return b;
}

__attribute__((always_inline)) inline
vector<vector<automata>> compute_buckets(vector<automata> const &automatas, vector<size_t> const &pos) {

        vector<vector<automata>> buckets(pos.size(), vector<automata>());

        for (automata const &a : automatas) {
                push_bucket(a, buckets, pos);
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

value bucket_elimination(vector<vector<automata>> &buckets, vector<size_t> const &order,
                         vector<size_t> const &pos, vector<size_t> const &domains);

#endif /* BE_HPP_ */
