#ifndef IO_H_
#define IO_H_

// threshold value to remove rows
#define THRESHOLD_VALUE (100000)

#ifndef THRESHOLD_VALUE
#define THRESHOLD_VALUE (numeric_limits<value>::max())
#endif

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <bits/stdc++.h>
#include <functional>
#include <sys/stat.h>

#include "types.hpp"
#include "order.hpp"
#include "libfa/fa.h"

using namespace std;

// Prints the content of an iterable type
template <typename T>
__attribute__((always_inline)) inline
void print_it(T const &vec, const char *name = nullptr, const char *format = nullptr, const char *after = nullptr) {

	if (name) printf("%s = [ ", name);
	else printf("[ ");
	for (auto i : vec) {
		if (format) { printf(format, i); printf(" "); }
		else cout << i << " ";
	}
	printf("]%s", (after) ? after : "\n");
}

void print_cost_table(cost const &c);

void cost_dot(cost const &c, const char *root_dir = ".");

void print_adj(vector<boost::dynamic_bitset<>> const &adj);

vector<boost::dynamic_bitset<>> read_adj(const char *wcsp);

pair<vector<cost>, vector<vector<size_t>>> read_costs_bin_vars(const char *wcsp, vector<size_t> const &pos = {}, value threshold = THRESHOLD_VALUE);

#endif /* IO_H_ */
