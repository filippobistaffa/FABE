#ifndef ORDER_HPP_
#define ORDER_HPP_

#include "types.hpp"

enum order_heur {
        WEIGHTED_MIN_FILL,
        MIN_FILL,
        MIN_INDUCED_WIDTH,
        MIN_DEGREE,
        RANDOM
};

template <typename T>
struct compare_pos {
        compare_pos(vector<T> const &pos) : pos(pos) {};
        bool operator()(size_t const &x, size_t const &y) {
                return pos[x] < pos[y];
        }
        vector<T> pos;
};

using namespace std;

vector<size_t> greedy_order(vector<vector<float>> const &adj, int order_heur = 0);

size_t induced_width(vector<vector<float>> const &adj, vector<size_t> const &order);

vector<size_t> read_pseudotree_order(const char *filename, vector<size_t> const &domains);

void export_order(vector<size_t> const &order, vector<size_t> const &domains, const char *output);

#endif /* ORDER_HPP_ */
