#ifndef BE_HPP_
#define BE_HPP_

#include "types.hpp"

enum be_inner_op {
        BE_SUM,
        BE_PROD
};

enum be_outer_op {
        BE_MIN,
        BE_MAX
};

vector<vector<automata>> compute_buckets(vector<automata> const &automatas, vector<size_t> const &pos);

value bucket_elimination(vector<vector<automata>> &buckets, int inner, int outer,
                         vector<size_t> const &order, vector<size_t> const &pos,
                         vector<size_t> const &domains, size_t ibound = 0);

#endif /* BE_HPP_ */
