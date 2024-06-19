#ifndef TYPES_HPP_
#define TYPES_HPP_

#include "libfa/fa.h"
#include <unordered_map>
#include <vector>

typedef double value;

typedef double weight;

struct automata {
    std::vector<size_t> vars;
    std::vector<size_t> domains;
    std::unordered_map<value, struct fa *> rows {};
};

struct table {
    std::vector<size_t> vars;
    std::vector<size_t> domains;
    std::vector<std::pair<std::vector<size_t>, value>> rows {};
};

enum instance {
    WCSP,
    MPE
};

#endif
