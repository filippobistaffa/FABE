#ifndef ORDER_HPP_
#define ORDER_HPP_

#include "types.hpp"

enum order_heur {
        O_WEIGHTED_MIN_FILL,
        O_MIN_FILL,
        O_MIN_INDUCED_WIDTH,
        O_MIN_DEGREE,
        O_RANDOM
};

enum tie_heur {
        T_UNIQUENESS,
        T_RANDOM
};

#include <cmath>                // fabs

template <typename T>
struct compare_vec {
        compare_vec(vector<T> const &vec) : vec(vec) {};
        bool operator()(size_t const &x, size_t const &y) {
                if constexpr (is_integral_v<T>) {
                        return vec[x] < vec[y];
                } else if (is_floating_point_v<T>) {
                        return vec[x] < vec[y] - numeric_limits<weight>::epsilon();
                }
        }
        vector<T> vec;
};

using namespace std;

vector<size_t> greedy_order(vector<vector<weight>> const &adj, int order_heur = O_WEIGHTED_MIN_FILL,
                            int tie_heur = T_UNIQUENESS);

size_t induced_width(vector<vector<weight>> const &adj, vector<size_t> const &order);

vector<size_t> read_pseudotree_order(const char *filename, vector<size_t> const &domains, vector<vector<size_t>> &anc);

void export_order(vector<size_t> const &order, vector<size_t> const &domains, const char *output);

#endif /* ORDER_HPP_ */
