#ifndef ORDER_HPP_
#define ORDER_HPP_

//#define MINWIDTH
//#define MININDUCEDWIDTH
#define MINFILL

#include <boost/dynamic_bitset.hpp>     // bitset

#include "bitset.hpp"
#include "types.hpp"
#include "io.hpp"

struct compare_pos {
        compare_pos(vector<size_t> const &pos) : pos(pos) {};
        bool operator()(size_t const &x, size_t const &y) {
                return pos[x] < pos[y];
        }
        vector<size_t> pos;
};

using namespace std;

vector<size_t> greedy_order(vector<boost::dynamic_bitset<>> const &adj);

size_t induced_width(vector<boost::dynamic_bitset<>> const &adj, vector<size_t> const &order, vector<size_t> const &pos);

vector<size_t> read_pseudotree_order(const char *filename);

//void export_order(vector<size_t> const &order, vector<size_t> const &domains, const char *output);

#endif /* ORDER_HPP_ */
