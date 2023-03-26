#ifndef BE_HPP_
#define BE_HPP_

#include "types.hpp"

#define QUANTISATION 1e10

std::vector<std::vector<automata>> compute_buckets(std::vector<automata> const &automatas, std::vector<size_t> const &pos);

value bucket_elimination(std::vector<std::vector<automata>> &buckets, bool quant,
                         std::vector<size_t> const &order, std::vector<size_t> const &pos,
                         std::vector<size_t> const &domains, size_t ibound = 0);

#endif /* BE_HPP_ */
