#ifndef TYPES_HPP_
#define TYPES_HPP_

#include <vector>
#include <boost/dynamic_bitset.hpp>

using namespace std;

typedef float value;

struct row {
        vector<size_t> a;
        struct fa *fa;
        value v;
};

struct cost {
        vector<size_t> vars;
        vector<size_t> domains;
        vector<row> rows;
};

#endif /* TYPES_HPP_ */
