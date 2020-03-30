#ifndef MINIMISE_HPP_
#define MINIMISE_HPP_

#include <cmath>                        // fabs
#include <sstream>                      // ostringstream

#include "types.hpp"
#include "libfa/fa.h"

using namespace std;

automata compute_automata(table const &t, value tolerance = 0);

#endif /* MINIMISE_HPP_ */
