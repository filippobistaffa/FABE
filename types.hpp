#ifndef TYPES_HPP_
#define TYPES_HPP_

#include <vector>
#include <boost/dynamic_bitset.hpp>

using namespace std;

typedef float value;

struct row {
        boost::dynamic_bitset<> a;
        value v;
        struct fa *fa;
};

struct cost {
        vector<size_t> vars;
        vector<size_t> bin_vars;
        vector<row> rows;
};

#endif /* TYPES_HPP_ */
