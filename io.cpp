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
                        cout << setw(width) << row.first[j] << " ";
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

#define PRECISION 2

void print_adj(vector<vector<weight>> const &adj) {

        const size_t n_vars = adj.size();
        const int var = 1 + floor(log10(n_vars - 1));
        const int column = max(PRECISION + 2, var);

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
                        cout << fixed << setprecision(PRECISION) << adj[i][j] << " ";
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

template <typename T1, typename T2>
static inline void remove_positions(vector<T1> &vec, vector<T2> const &mask) {

        size_t i = 0;
        for (auto b : mask) {
                if (b) {
                        vec.erase(vec.begin() + i);
                } else {
                        i++;
                }
        }
}

static inline tuple<vector<size_t>, size_t, vector<size_t>> compute_evidence(vector<size_t> const &domains) {

        const size_t n_vars = domains.size();
        vector<size_t> ev(n_vars);
        for (size_t i = 0; i < n_vars; ++i) {
                if (domains[i] == 1) {
                        ev[i] = 1;
                }
        }
        vector<size_t> pfx(n_vars);
        partial_sum(ev.begin(), ev.end(), pfx.begin());
        vector<size_t> map(n_vars);
        iota(map.begin(), map.end(), 0);
        transform(map.begin(), map.end(), pfx.begin(), map.begin(), minus<size_t>());
        return make_tuple(ev, accumulate(ev.begin(), ev.end(), 0), map);
}

#define SKIP_LINE f.ignore(numeric_limits<streamsize>::max(), '\n')

static inline pair<vector<size_t>, vector<vector<weight>>> read_domains_adj_wcsp(const char *wcsp) {

        ifstream f(wcsp);
        const auto [ n_vars, max_domain, n_tables ] = tokenize<size_t, 1, 3>(f);
        auto domains = tokenize<size_t>(f);
        const auto [ ev, n_ev, map ] = compute_evidence(domains);
        vector<vector<weight>> adj(n_vars - n_ev, vector<weight>(n_vars - n_ev));
        vector<vector<weight>> tot(n_vars - n_ev, vector<weight>(n_vars - n_ev));

        for (size_t i = 0; i < n_tables; ++i) {
                auto temp = tokenize<value>(f);
                vector<size_t> vars(temp.begin() + 1, temp.begin() + temp[0] + 1);
                size_t n_rows = 1;
                for (auto it = vars.begin(); it != vars.end(); ++it) {
                        n_rows *= domains[*it];
                }
                vector<value> values(1, temp[temp[0] + 1]);
                for (size_t j = 0; j < temp[temp[0] + 2]; ++j) {
                        auto row = tokenize<value>(f);
                        values.push_back(row.back());
                }
                sort(values.begin(), values.end());
                const weight u = unique(values.begin(), values.end()) - values.begin();
                for (auto it = vars.begin(); it != vars.end(); ++it) {
                        if (!ev[*it]) {
                                for (auto it1 = it + 1; it1 != vars.end(); ++it1) {
                                        if (!ev[*it1]) {
                                                adj[map[*it]][map[*it1]] += u / n_rows;
                                                adj[map[*it1]][map[*it]] += u / n_rows;
                                                tot[map[*it]][map[*it1]]++;
                                                tot[map[*it1]][map[*it]]++;
                                        }
                                }
                        }
                }
        }

        for (size_t i = 0; i < n_vars - n_ev; ++i) {
                for (size_t j = 0; j < n_vars - n_ev; ++j) {
                        if (tot[i][j]) {
                                adj[i][j] /= tot[i][j];
                        } else {
                                adj[i][j] = 0;
                        }
                }
        }

        f.close();
        remove_positions(domains, ev);
        return make_pair(domains, adj);
}

static inline pair<vector<size_t>, vector<vector<weight>>> read_domains_adj_uai(const char *uai) {

        ifstream f(uai);
        SKIP_LINE;
        auto [ n_vars ] = tokenize<size_t, 0, 1>(f);
        auto domains = tokenize<size_t>(f);
        auto [ n_tables ] = tokenize<size_t, 0, 1>(f);
        const auto [ ev, n_ev, map ] = compute_evidence(domains);
        vector<vector<weight>> adj(n_vars - n_ev, vector<weight>(n_vars - n_ev));
        vector<vector<weight>> tot(n_vars - n_ev, vector<weight>(n_vars - n_ev));
        vector<vector<size_t>> vars(n_tables);

        for (size_t i = 0; i < n_tables; ++i) {
                auto temp = tokenize<value>(f);
                vars[i] = vector<size_t>(temp.begin() + 1, temp.begin() + temp[0] + 1);
        }

        for (size_t i = 0; i < n_tables; ++i) {
                auto [ n_rows ] = tokenize<size_t, 0, 1>(f);
                vector<value> values;
                while (values.size() < n_rows) {
                        auto temp = tokenize<value>(f);
                        values.insert(values.end(), temp.begin(), temp.end());
                }
                sort(values.begin(), values.end());
                const weight u = unique(values.begin(), values.end()) - values.begin();
                for (auto it = vars[i].begin(); it != vars[i].end(); ++it) {
                        if (!ev[*it]) {
                                for (auto it1 = it + 1; it1 != vars[i].end(); ++it1) {
                                        if (!ev[*it1]) {
                                                adj[map[*it]][map[*it1]] += u / n_rows;
                                                adj[map[*it1]][map[*it]] += u / n_rows;
                                                tot[map[*it]][map[*it1]]++;
                                                tot[map[*it1]][map[*it]]++;
                                        }
                                }
                        }
                }
        }

        for (size_t i = 0; i < n_vars - n_ev; ++i) {
                for (size_t j = 0; j < n_vars - n_ev; ++j) {
                        if (tot[i][j]) {
                                adj[i][j] /= tot[i][j];
                        } else {
                                adj[i][j] = 0;
                        }
                }
        }

        f.close();
        remove_positions(domains, ev);
        return make_pair(domains, adj);
}

pair<vector<size_t>, vector<vector<weight>>> read_domains_adj(const char *instance, int type) {

        if (type == WCSP) {
                return read_domains_adj_wcsp(instance);
        } else {
                return read_domains_adj_uai(instance);
        }
}

void preallocate_rows(table &t, value def) {

        const size_t n_rows = accumulate(t.domains.begin(), t.domains.end(), 1, multiplies<size_t>());
        t.rows.resize(n_rows);

        for (size_t i = 0; i < n_rows; ++i) {
                t.rows[i] = make_pair(get_combination(i, t.domains), def);
        }
}

__attribute__((always_inline)) inline
void remove_threshold(table &t, value threshold) {

        size_t r = 0;

        for (auto it = t.rows.rbegin(); it != t.rows.rend(); ++it) {
                if (it->second >= threshold) {
                        r++;
                } else {
                        break;
                }
        }

        if (r) {
                t.rows.resize(t.rows.size() - r);
        }
}

static inline value project_evidence(table &t, vector<size_t> const &ev, vector<size_t> const &map) {

        vector<size_t> mask(t.vars.size());
        for (size_t i = 0; i < t.vars.size(); ++i) {
                mask[i] = ev[t.vars[i]];
        }

        remove_positions(t.vars, mask);
        if (t.vars.size() == 0) {
                return t.rows.front().second;
        }

        remove_positions(t.domains, mask);
        for (auto &row : t.rows) {
                remove_positions(row.first, mask);                
        }

        for (size_t i = 0; i < t.vars.size(); ++i) {
                t.vars[i] = map[t.vars[i]];            
        }

        return 0;
}

static inline pair<vector<table>, value> read_tables_wcsp(const char *wcsp, value threshold) {

        ifstream f(wcsp);
        const auto [ n_vars, max_domain, n_tables ] = tokenize<size_t, 1, 3>(f);
        const auto domains = tokenize<size_t>(f);
        const auto [ ev, n_ev, map ] = compute_evidence(domains);
        vector<table> tables;
        value evid_value = 0;

        for (size_t i = 0; i < n_tables; ++i) {

                table t;
                auto temp = tokenize<value>(f);
                auto vars = vector<size_t>(temp.begin() + 1, temp.begin() + temp[0] + 1);
                t.vars = vars;
                sort(t.vars.begin(), t.vars.end(), greater<size_t>());

                for (auto var : t.vars) {
                        t.domains.push_back(domains[var]);
                }

                vector<size_t> pos(t.vars.size());
                for (size_t i = 0; i < t.vars.size(); ++i) {
                        pos[i] = find(t.vars.begin(), t.vars.end(), vars[i]) - t.vars.begin();
                }

                vector<size_t> pfx_prod(t.domains.size());
                exclusive_scan(t.domains.begin(), t.domains.end(), pfx_prod.begin(), 1, multiplies<>{});
                preallocate_rows(t, temp[temp[0] + 1]);

                for (size_t j = 0; j < temp[temp[0] + 2]; ++j) {
                        auto row = tokenize<value>(f);
                        const value val = row[row.size() - 1];
                        size_t idx = 0;
                        for (size_t k = 0; k < row.size() - 1; ++k) {
                                idx += row[k] * pfx_prod[pos[k]];
                        }
                        t.rows[idx].second = val;
                }

                // sort according to values
                sort(t.rows.begin(), t.rows.end(), [](pair<vector<size_t>, value> const &x,
                                                      pair<vector<size_t>, value> const &y)
                                                      { return (x.second < y.second); });
                remove_threshold(t, threshold);
                evid_value += project_evidence(t, ev, map);
                if (t.vars.size()) {
                        tables.push_back(t);
                }
        }

        f.close();
        return make_pair(tables, evid_value);
}

static inline pair<vector<table>, value> read_tables_uai(const char *uai, value threshold) {

        ifstream f(uai);
        SKIP_LINE;
        SKIP_LINE;
        const auto domains = tokenize<size_t>(f);
        auto [ n_tables ] = tokenize<size_t, 0, 1>(f);
        const auto [ ev, n_ev, map ] = compute_evidence(domains);
        vector<table> tables(n_tables);
        value evid_value = 0;

        for (size_t i = 0; i < n_tables; ++i) {
                table t;
                auto temp = tokenize<value>(f);
                t.vars = vector<size_t>(temp.begin() + 1, temp.begin() + temp[0] + 1);
                tables[i] = t;
        }

        for (size_t i = 0; i < n_tables; ++i) {

                auto vars = tables[i].vars;
                reverse(vars.begin(), vars.end());
                vector<size_t> orig_dom;
                sort(tables[i].vars.begin(), tables[i].vars.end(), greater<size_t>());

                for (auto var : vars) {
                        orig_dom.push_back(domains[var]);
                }

                for (auto var : tables[i].vars) {
                        tables[i].domains.push_back(domains[var]);
                }

                vector<size_t> pos(tables[i].vars.size());
                for (size_t j = 0; j < tables[i].vars.size(); ++j) {
                        pos[j] = find(tables[i].vars.begin(), tables[i].vars.end(), vars[j]) - tables[i].vars.begin();
                }

                vector<size_t> pfx_prod(tables[i].domains.size());
                exclusive_scan(tables[i].domains.begin(), tables[i].domains.end(), pfx_prod.begin(), 1, multiplies<>{});

                auto [ n_rows ] = tokenize<size_t, 0, 1>(f);
                vector<value> values;

                while (values.size() < n_rows) {
                        auto temp = tokenize<value>(f);
                        values.insert(values.end(), temp.begin(), temp.end());
                }

                preallocate_rows(tables[i]);
                //cout << vec2str(values, "values") << endl;

                for (size_t j = 0; j < n_rows; ++j) {
                        auto row = get_combination(j, orig_dom);
                        size_t idx = 0;
                        for (size_t k = 0; k < row.size(); ++k) {
                                idx += row[k] * pfx_prod[pos[k]];
                        }
                        tables[i].rows[idx].second = -log(values[j]);
                }

                // sort according to values
                sort(tables[i].rows.begin(), tables[i].rows.end(), [](pair<vector<size_t>, value> const &x,
                                                                      pair<vector<size_t>, value> const &y)
                                                                      { return (x.second < y.second); });
                remove_threshold(tables[i], threshold);
        }

        f.close();

        for (size_t i = 0; i < tables.size();) {
                evid_value += project_evidence(tables[i], ev, map);
                if (tables[i].vars.size() == 0) {
                        tables.erase(tables.begin() + i);
                        puts("Removed one table");
                } else {
                        i++;
                }
        }

        return make_pair(tables, evid_value);
}

pair<vector<table>, value> read_tables(const char *instance, int type, value threshold) {

        if (type == WCSP) {
                return read_tables_wcsp(instance, threshold);
        } else {
                return read_tables_uai(instance, threshold);
        }
}

/*void export_wcsp(vector<vector<automata>> buckets, vector<size_t> const &domains, const char *wcsp) {

        ostringstream oss;
        size_t n_funcs = 0;

        for (auto bucket : buckets) {
                n_funcs += bucket.size();
        }

        oss << wcsp << " ";
        oss << domains.size() << " ";
        oss << *max_element(domains.begin(), domains.end()) << " ";
        oss << n_funcs << " ";
        oss << numeric_limits<int>::max() << endl;
        oss << vec2str(domains, nullptr, " ", nullptr, nullptr) << endl;

        for (auto bucket : buckets) {
                for (auto a : bucket) {
                        const auto t = compute_table(a);
                        oss << t.vars.size() << " ";
                        oss << vec2str(t.vars, nullptr, " ", nullptr, nullptr) << " ";
                        oss << numeric_limits<int>::max() << " ";
                        oss << t.rows.size() << endl;
                        for (auto r : t.rows) {
                                oss << vec2str(r.first, nullptr, " ", nullptr, nullptr) << " " << r.second << endl;
                        }
                }
        }

        ofstream f(wcsp);
        f << oss.str();
        f.close();
}*/
