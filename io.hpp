#ifndef IO_H_
#define IO_H_

// threshold value to remove rows
#define THRESHOLD_VALUE (100000)
//#define THRESHOLD_VALUE (100)
//#define THRESHOLD_VALUE (10)

#ifndef THRESHOLD_VALUE
#define THRESHOLD_VALUE (numeric_limits<value>::max())
#endif

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

using namespace std;

template <typename T>
__attribute__((always_inline)) inline
string vec2str(T const &vec, const char *name = nullptr, const char *after = nullptr) {

        ostringstream oss;
	if (name) {
	        oss << name << " = ";
        }
        oss << "[ ";
	for (auto i : vec) {
		oss << i << " ";
	}
        oss << "]";
	if (after) {
                cout << after;
	}
	return oss.str();
}

void print_table(table const &t);

void automata_dot(automata const &c, const char *root_dir = ".");

void print_adj(vector<boost::dynamic_bitset<>> const &adj);

vector<boost::dynamic_bitset<>> read_adj(const char *wcsp);

pair<vector<size_t>, vector<table>> read_domains_tables(const char *wcsp, value threshold = THRESHOLD_VALUE);

#endif /* IO_H_ */
