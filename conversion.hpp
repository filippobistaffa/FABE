#ifndef CONVERSION_HPP_
#define CONVERSION_HPP_

#include <cmath>        // fabs
#include <sstream>      // ostringstream
#include <iostream>     // cout
#include <cstring>      // strlen
#include <numeric>      // accumulate

#include "types.hpp"
#include "libfa/fa.h"

using namespace std;

automata compute_automata(table const &t, value tolerance = 0);

table compute_table(automata const &a);

#endif /* CONVERSION_HPP_ */
