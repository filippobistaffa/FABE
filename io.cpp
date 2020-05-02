#include "io.hpp"

#include <iostream>             // cout
#include <sstream>              // ostringstream
#include <fstream>              // ifstream, getline
#include <string.h>             // strdup, strdup
#include <math.h>               // log10
#include <iomanip>              // setw
#include <numeric>              // accumulate
#include <sys/stat.h>           // filesystem
#include <algorithm>            // max_element, sort
#include <linux/limits.h>       // PATH_MAX
#include <unistd.h>             // getcwd
#include <cassert>              // assert

#include "util.hpp"
#include "libfa/fa.h"
#include "order.hpp"
#include "conversion.hpp"

void print_table(table const &t) {

        const size_t width = max(1.0, 1 + floor(log10(*max_element(t.vars.begin(), t.vars.end()))));

        for (auto var : t.vars) {
                cout << setw(width) << var << " ";
        }
        cout << endl;

        for (size_t i = 0; i < (width + 1) * t.vars.size() - 1; ++i) {
                cout << "-";
        }
        cout << endl;

        for (auto const &row : t.rows) {
                for (size_t j = 0; j < t.vars.size(); ++j) {
                        cout << setw(width) << ALPHABET[row.first[j]] << " ";
                }
                cout << "= " << row.second << endl;
        }
}

void automata_dot(automata const &a, const char *root_dir) {

        char cwd[PATH_MAX];
        getcwd(cwd, sizeof(cwd));
        chdir(root_dir);
        auto folder = vec2str(a.vars);
        mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        chdir(folder.c_str());

        for (auto const &[ v, fa ] : a.rows) {
                fa_make_dot(fa, "%.5f.dot", v);
        }

        chdir(cwd);
}

void print_adj(vector<vector<weight>> const &adj) {

        const size_t n_vars = adj.size();
        const int var = 1 + floor(log10(n_vars - 1));
        const int column = var;

        for (size_t i = 0; i < var + 3; ++i) {
                cout << " ";
        }

        for (size_t i = 0; i < n_vars; ++i) {
                cout << setw(column) << i << " ";
        }
        cout << endl;

        for (size_t i = 0; i < var + 1; ++i) {
                cout << " ";
        }

        cout << "+";
        for (size_t i = 0; i < (column + 1) * n_vars; ++i) {
                cout << "-";
        }
        cout << endl;

        for (size_t i = 0; i < n_vars; ++i) {
                cout << setw(var) << i << " | ";
                for (size_t j = 0; j < n_vars; ++j) {
                        cout << adj[i][j] << " ";
                }
                cout << endl;
        }
}

string trim(string const &str, string const &whitespace = " \t") {

        const auto strBegin = str.find_first_not_of(whitespace);
        if (strBegin == std::string::npos)
                return ""; // no content
        const auto strEnd = str.find_last_not_of(whitespace);
        const auto strRange = strEnd - strBegin + 1;
        return str.substr(strBegin, strRange);
}

template <typename T>
vector<T> tokenize(ifstream &f) {

        string str;
        getline(f, str);
        while (str.empty()) {
                getline(f, str);
        }
        str = trim(str);
        char *dup = strdup(str.c_str());
        const char *sep = (strstr(dup, " ") != NULL) ? " " : "\t";
        char *token = strtok(dup, sep);
        vector<T> v;

        while (token != NULL) {
                if constexpr (is_integral_v<T>) {
                        v.push_back(atoi(token));
                } else if (is_floating_point_v<T>) {
                        v.push_back(atof(token));
                }
                token = strtok(NULL, sep);
        }

        free(dup);
        return v;
}

template <typename T, size_t SKIP, size_t N>
array<T, N> tokenize(ifstream &f) {

        string str;
        getline(f, str);
        while (str.empty()) {
                getline(f, str);
        }
        str = trim(str);
        char *dup = strdup(str.c_str());
        const char *sep = (strstr(dup, " ") != NULL) ? " " : "\t";
        char *token = strtok(dup, sep);
        array<T, N> a;
        size_t skip = SKIP;
        size_t i = 0;

        while (token != NULL) {
                if (skip) {
                        skip--;
                } else {
                        if (i < N) {
                                if constexpr (is_integral_v<T>) {
                                        a[i++] = atoi(token);
                                } else if (is_floating_point_v<T>) {
                                        a[i++] = atof(token);
                                }
                        } else {
                                break;
                        }
                }
                token = strtok(NULL, sep);
        }

        free(dup);
        return a;
}

#define SKIP_LINE f.ignore(numeric_limits<streamsize>::max(), '\n')

pair<vector<size_t>, vector<vector<weight>>> read_domains_adj(const char *wcsp) {

        ifstream f(wcsp);
        const auto [ n_vars, max_domain, n_tables ] = tokenize<size_t, 1, 3>(f);
        assert(max_domain <= ALPHABET_LENGTH);
        const auto domains = tokenize<size_t>(f);
        vector<vector<weight>> adj(n_vars, vector<weight>(n_vars));

        for (size_t i = 0; i < n_tables; ++i) {
                auto temp = tokenize<value>(f);
                vector<size_t> vars(temp.begin() + 1, temp.begin() + temp[0] + 1);
                for (auto it = vars.begin(); it != vars.end(); ++it) {
                        for (auto it1 = it + 1; it1 != vars.end(); ++it1) {
                                adj[*it][*it1] = 1;
                                adj[*it1][*it] = 1;
                        }
                }
        }

        f.close();
        return make_pair(domains, adj);
}

vector<table> read_tables(const char *wcsp) {

        ifstream f(wcsp);
        const auto [ n_vars, max_domain, n_tables ] = tokenize<size_t, 1, 3>(f);
        const auto domains = tokenize<size_t>(f);
        vector<table> tables(n_tables);

        for (size_t i = 0; i < n_tables; ++i) {

                table t;
                auto temp = tokenize<value>(f);
                t.vars = vector<size_t>(temp.begin() + 1, temp.begin() + temp[0] + 1);

                for (auto var : t.vars) {
                        t.domains.push_back(domains[var]);
                }

                for (size_t j = 0; j < temp[temp[0] + 2]; ++j) {
                        auto row = tokenize<value>(f);
                        const value val = row[row.size() - 1];
                        t.rows.push_back(make_pair(vector<size_t>(row.begin(), row.end() - 1), val));
                }

                // sort according to values
                sort(t.rows.begin(), t.rows.end(), [](pair<vector<size_t>, value> const &x,
                                                      pair<vector<size_t>, value> const &y)
                                                      { return (x.second < y.second); });
                tables[i] = t;
        }

        f.close();
        return tables;
}
