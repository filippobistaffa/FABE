#ifndef TYPES_HPP_
#define TYPES_HPP_

#define ALPHABET "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"

#include <unordered_map>
#include <vector>
using namespace std;

typedef float value;

struct automata {
        vector<size_t> vars;
        vector<size_t> domains;
        unordered_map<value, struct fa *> rows;
};

struct table {
        vector<size_t> vars;
        vector<size_t> domains;
        vector<pair<vector<size_t>, value>> rows;
};

#endif /* TYPES_HPP_ */
