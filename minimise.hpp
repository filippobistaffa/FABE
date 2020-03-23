#ifndef MINIMISE_HPP_
#define MINIMISE_HPP_

#include <sstream>
#include <iostream>
#include <cmath>
#include <string.h>

#include "types.hpp"
#include "fa.h"

using namespace std;

cost compress_clusters(cost const &in, value tolerance = 0);

#endif /* MINIMISE_HPP_ */
