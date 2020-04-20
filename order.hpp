#ifndef ORDER_HPP_
#define ORDER_HPP_

#include <boost/dynamic_bitset.hpp>     // bitset
#include "bitset.hpp"
#include "types.hpp"
#include "io.hpp"

enum order_heur {
        WEIGHTED_MIN_FILL,
        MIN_FILL,
        MIN_INDUCED_WIDTH,
        MIN_DEGREE,
        RANDOM
};

struct compare_pos {
        compare_pos(vector<size_t> const &pos) : pos(pos) {};
        bool operator()(size_t const &x, size_t const &y) {
                return pos[x] < pos[y];
        }
        vector<size_t> pos;
};

using namespace std;

vector<size_t> greedy_order(int order_heur, vector<boost::dynamic_bitset<>> const &adj, vector<float> const &weights);

size_t induced_width(vector<boost::dynamic_bitset<>> const &adj, vector<size_t> const &order, vector<size_t> const &pos);

vector<size_t> read_pseudotree_order(const char *filename, vector<size_t> const &domains);

//void export_order(vector<size_t> const &order, vector<size_t> const &domains, const char *output);

#endif /* ORDER_HPP_ */
