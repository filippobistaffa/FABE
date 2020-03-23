#include "io.hpp"

void print_cost_table(cost const &c) {

        const auto width = max(1.0, 1 + floor(log10(*max_element(c.bin_vars.begin(), c.bin_vars.end()))));

        for (auto var : c.bin_vars) {
                cout << setw(width) << var << " ";
        }
        cout << endl;

        for (auto i = 0; i < (width + 1) * c.bin_vars.size() - 1; ++i) {
                cout << "-";
        }
        cout << endl;

        for (auto row : c.rows) {
                for (auto j = 0; j < c.bin_vars.size(); ++j) {
                        for (auto i = 0; i < width - 1; ++i) {
                                cout << " ";
                        }
                        cout << row.a[j] << " ";
                }
                cout << "= " << row.v << endl;
        }
}

void print_adj(vector<boost::dynamic_bitset<>> const &adj) {

        const auto n_vars = adj.size();
        const auto width = 1 + floor(log10(n_vars - 1));

        for (auto i = 0; i < width + 3; ++i) {
                cout << " ";
        }

        for (auto i = 0; i < n_vars; ++i) {
                cout << setw(width) << i << " ";
        }
        cout << endl;

        for (auto i = 0; i < width + 1; ++i) {
                cout << " ";
        }

        cout << "+";
        for (auto i = 0; i < (width + 1) * n_vars; ++i) {
                cout << "-";
        }
        cout << endl;

        for (auto i = 0; i < n_vars; ++i) {
                cout << setw(width) << i << " | ";
                for (auto j = 0; j < n_vars; ++j) {
                        for (auto k = 0; k < width - 1; ++k) {
                                cout << " ";
                        }
                        cout << adj[i][j] << " ";
                }
                cout << endl;
        }
}

__attribute__((always_inline)) inline
void preallocate_rows(cost &c, value def) {

        if (c.bin_vars.size() > 64) {
                throw overflow_error("Cannot initialise more than 64 variables");
        }

        for (size_t i = 0; i < powl(2, c.bin_vars.size()); ++i) {
                row r = {
                        boost::dynamic_bitset<>(c.bin_vars.size(), i),
                        def
                };
                c.rows.push_back(r);
        }
}

template <typename T1, typename T2, typename T3>
__attribute__((always_inline)) inline
void exclusive_scan(T1 begin, T1 end, T1 out, T2 init, T3 op) {

        *out = init;
        for (auto it = begin; it != prev(end); ++it, ++out) {
                *(out + 1) = op(*it, *out);
        }
}

template <typename T>
vector<T> tokenize(ifstream &f, const char *sep = " ", size_t skip = 0, size_t max_n = numeric_limits<size_t>::max()) {

        string str;
        getline(f, str);
        char *dup = strdup(str.c_str());
        char *token = strtok(dup, sep);
        vector<T> v;

        while (token != NULL) {
                if (skip) {
                        skip--;
                } else {
                        if (max_n) {
                                if constexpr (is_integral_v<T>) {
                                        v.push_back(atoi(token));
                                } else if (is_floating_point_v<T>) {
                                        v.push_back(atof(token));
                                }
                                max_n--;
                        } else {
                                break;
                        }
                }
                token = strtok(NULL, sep);
        }

        free(dup);
        return v;
}

vector<boost::dynamic_bitset<>> read_adj(const char *wcsp) {

        ifstream f(wcsp);
        auto vars_costs = tokenize<size_t>(f, " ", 1, 3);
        f.ignore(numeric_limits<streamsize>::max(), '\n');
        vector<boost::dynamic_bitset<>> adj(vars_costs[0]);

        for (auto i = 0; i < vars_costs[0]; ++i) {
                adj[i] = boost::dynamic_bitset<>(vars_costs[0]);
        }

        for (auto i = 0; i < vars_costs[2]; ++i) {
                auto temp = tokenize<value>(f);
                vector<size_t> vars;
                for (auto it = temp.begin() + 1; it != temp.begin() + temp[0] + 1; ++it) {
                        vars.push_back(*it);
                }
                for (auto it = vars.begin(); it != vars.end(); ++it) {
                        for (auto it1 = it + 1; it1 != vars.end(); ++it1) {
                                adj[*it].set(*it1);
                                adj[*it1].set(*it);
                        }
                }
                for (auto j = 0; j < temp[temp[0] + 2]; ++j) {
                        f.ignore(numeric_limits<streamsize>::max(), '\n');
                }
        }

        f.close();
        return adj;
}

__attribute__((always_inline)) inline
void remove_threshold(cost &c, value threshold) {

        size_t r = 0;

        for (auto it = c.rows.rbegin(); it != c.rows.rend(); ++it) {
                if (it->v >= threshold) {
                        r++;
                } else {
                        break;
                }
        }

        if (r) {
                c.rows.resize(c.rows.size() - r);
        }
}

pair<vector<cost>, vector<vector<size_t>>> read_costs_bin_vars(const char *wcsp, vector<size_t> const &pos, value threshold) {

        ifstream f(wcsp);
        auto vars_costs = tokenize<size_t>(f, " ", 1, 3);
        auto domains = tokenize<size_t>(f);
        vector<size_t> bin_per_var;
        vector<cost> costs;

        // compute number of bits needed for each domain
        for (auto domain : domains) {
                bin_per_var.push_back(1 + floor(log2(domain - 1)));
        }

        // compute exclusive prefix sum
        vector<size_t> pfx(bin_per_var.size());
        exclusive_scan(bin_per_var.begin(), bin_per_var.end(), pfx.begin(), 0, plus<>{});
        vector<vector<size_t>> bin_vars(bin_per_var.size());

        for (auto i = 0; i < vars_costs[0]; ++i) {
                vector<size_t> tmp;
                for (auto j = pfx[i]; j != pfx[i] + bin_per_var[i]; ++j) {
                        tmp.push_back(j);
                }
                bin_vars[i] = tmp;
        }

        for (auto i = 0; i < vars_costs[2]; ++i) {

                cost c;
                auto temp = tokenize<value>(f);
                vector<size_t> vars, binary_domains;

                for (auto it = temp.begin() + 1; it != temp.begin() + temp[0] + 1; ++it) {
                        vars.push_back(*it);
                }

                c.vars = vars;
                #ifdef SORT_INPUT_TABLES
                sort(c.vars.begin(), c.vars.end(), compare_pos(pos));
                #endif
                vector<size_t> map;

                for (auto var : vars) {
                        map.push_back(find(c.vars.begin(), c.vars.end(), var) - c.vars.begin());
                }

                //print_it(vars);
                //print_it(c.vars);
                //print_it(map);

                for (auto var : c.vars) {
                        binary_domains.push_back(pow(2, bin_per_var[var]));
                        for (auto bin_var : bin_vars[var]) {
                                c.bin_vars.push_back(bin_var);
                        }
                }

                vector<size_t> pfx_prod(binary_domains.size());
                exclusive_scan(binary_domains.begin(), binary_domains.end(), pfx_prod.begin(), 1, multiplies<>{});
                //print_it(binary_domains);
                //print_it(pfx_prod);
                //print_it(c.variables);
                preallocate_rows(c, temp[temp[0] + 1]);

                for (auto j = 0; j < temp[temp[0] + 2]; ++j) {
                        auto row = tokenize<value>(f);
                        const value val = row[row.size() - 1];
                        auto idx = 0;
                        for (auto k = 0; k < row.size() - 1; ++k) {
                                idx += row[k] * pfx_prod[map[k]];
                        }
                        c.rows[idx].v = val;
                        //cout << idx << endl;
                        //print_it(row);
                }

                // sort according to values
                sort(c.rows.begin(), c.rows.end(), [](const row &x, const row &y) { return (x.v < y.v); });
                remove_threshold(c, threshold);
                costs.push_back(c);
        }

        f.close();
        return make_pair(costs, bin_vars);
}
