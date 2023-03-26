#ifndef IO_HPP_
#define IO_HPP_

#include "types.hpp"

#include <string>       // string
#include <sstream>      // ostringstream

template <typename T>
__attribute__((always_inline)) inline
std::string vec2str(T const &vec, const char *name = nullptr, const char *sep = " ",
               const char *open = "[ ", const char *close = " ]") {

    std::ostringstream oss;
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

void print_adj(std::vector<std::vector<weight>> const &adj);

void preallocate_rows(table &t, value def = 0);

std::pair<std::vector<size_t>, std::vector<std::vector<weight>>> read_domains_adj(const char *instance, int type);

std::vector<table> read_tables(const char *instance, int type, std::vector<size_t> const &pos, value threshold);

//void export_wcsp(std::vector<std::vector<automata>> buckets, std::vector<size_t> const &domains, const char *wcsp);

#endif /* IO_HPP_ */
