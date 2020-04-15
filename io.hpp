#ifndef IO_HPP_
#define IO_HPP_

// print variables' positions instead of variables
//#define PRINT_VAR_POS

#include <iostream>                     // cout
#include <sstream>                      // ostringstream
#include <fstream>                      // ifstream, getline
#include <string.h>                     // strdup, strdup
#include <math.h>                       // log10
#include <iomanip>                      // setw
#include <numeric>                      // accumulate
#include <boost/dynamic_bitset.hpp>     // bitset
#include <sys/stat.h>                   // filesystem

#include "types.hpp"
#include "util.hpp"
#include "libfa/fa.h"
#include "order.hpp"
#include "conversion.hpp"

using namespace std;

template <typename T>
__attribute__((always_inline)) inline
string vec2str(T const &vec, const char *name = nullptr, const char *sep = " ",
               const char *open = "[ ", const char *close = " ]") {

        ostringstream oss;
	if (name) { oss << name << " = "; }
        if (open) { oss << open; }
        if (vec.size() > 1) {
	        for (auto it = vec.begin(); it != prev(vec.end()); ++it) {
		        oss << *it << sep;
	        }
        }
        if (vec.size() > 0) { oss << vec.back(); }
        if (close) { oss << close; }
	return oss.str();
}

void print_table(table const &t);

void automata_dot(automata const &c, const char *root_dir = ".");

void print_adj(vector<boost::dynamic_bitset<>> const &adj);

vector<boost::dynamic_bitset<>> read_adj(const char *instance, int type);

pair<vector<size_t>, vector<table>> read_domains_tables(const char *instance, int type, vector<size_t> const &pos, value threshold);

void export_wcsp(vector<vector<automata>> buckets, vector<size_t> const &domains, const char *wcsp);

#endif /* IO_HPP_ */
