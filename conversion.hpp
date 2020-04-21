#ifndef CONVERSION_HPP_
#define CONVERSION_HPP_

#include "types.hpp"

automata compute_automata(table const &t, value tolerance = 0);

table compute_table(automata const &a);

#endif /* CONVERSION_HPP_ */
