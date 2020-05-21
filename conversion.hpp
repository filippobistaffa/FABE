#ifndef CONVERSION_HPP_
#define CONVERSION_HPP_

#include "types.hpp"

pair<automata, value> compute_automata(table const &t, value tolerance = 0);

table compute_table(automata const &a);

pair<double *, size_t> compute_cpt(automata const &a);

#endif /* CONVERSION_HPP_ */
