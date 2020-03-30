#ifndef BE_HPP_
#define BE_HPP_

#include <numeric>      // accumulate

#include "libfa/fa.h"
#include "types.hpp"
#include "bitset.hpp"
#include "util.hpp"
#include "order.hpp"

#include "io.hpp"

automata join_bucket(vector<automata> const &bucket, vector<size_t> domains);

void reduce_var(automata const &cost, size_t var);

template <typename T>
vector<vector<T>> compute_buckets(vector<T> const &costs, vector<size_t> const &pos) {

        vector<vector<T>> buckets(pos.size(), vector<T>());

        for (T c : costs) {
                const auto max_var = *max_element(c.vars.begin(), c.vars.end(), compare_pos(pos));
                buckets[max_var].push_back(c);
        }

        return buckets;
}

#endif /* BE_HPP_ */
