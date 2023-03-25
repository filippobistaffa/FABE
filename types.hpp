#ifndef TYPES_HPP_
#define TYPES_HPP_

#include "libfa/fa.h"
#include <unordered_map>
#include <vector>

using namespace std;

typedef double value;

typedef double weight;

struct automata {
        vector<size_t> vars;
        vector<size_t> domains;
        unordered_map<value, struct fa *> rows {};
};

struct table {
        vector<size_t> vars;
        vector<size_t> domains;
        vector<pair<vector<size_t>, value>> rows {};
};

enum instance {
        WCSP,
        MPE
};

#endif /* TYPES_HPP_ */
