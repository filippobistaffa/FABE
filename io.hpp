#ifndef IO_HPP_
#define IO_HPP_

#include "types.hpp"

void print_table(table const &t);

void automata_dot(automata const &c, const char *root_dir = ".");

void print_adj(std::vector<std::vector<weight>> const &adj);

void preallocate_rows(table &t, value def = 0);

std::pair<std::vector<size_t>, std::vector<std::vector<weight>>> read_domains_adj(const char *instance, int type);

std::vector<table> read_tables(const char *instance, int type, std::vector<size_t> const &pos, value threshold);

//void export_wcsp(std::vector<std::vector<automata>> buckets, std::vector<size_t> const &domains, const char *wcsp);

#endif
