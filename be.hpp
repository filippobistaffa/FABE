#ifndef BE_HPP_
#define BE_HPP_

#include "types.hpp"
#include <cfloat>

#define QUANTISATION 1e10
//#define QUANTISATION 1e20
//#define QUANTISATION 1e30
//#define QUANTISATION 1e40
//#define QUANTISATION 1e50
//#define QUANTISATION 1e50

vector<vector<automata>> compute_buckets(vector<automata> const &automatas, vector<size_t> const &pos);

value bucket_elimination(vector<vector<automata>> &buckets, bool quant,
                         vector<size_t> const &order, vector<size_t> const &pos,
                         vector<size_t> const &domains, size_t ibound = 0);

#endif /* BE_HPP_ */
