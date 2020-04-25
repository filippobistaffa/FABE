#ifndef CONVERSION_HPP_
#define CONVERSION_HPP_

#include "types.hpp"

pair<automata, value> compute_automata(table const &t, size_t max_k = numeric_limits<size_t>::max());

table compute_table(automata const &a);

#endif /* CONVERSION_HPP_ */
