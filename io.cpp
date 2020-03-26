#include "io.hpp"

void print_cost_table(cost const &c) {

        const auto width = max(1.0, 1 + floor(log10(*max_element(c.vars.begin(), c.vars.end()))));

        for (auto var : c.vars) {
                cout << setw(width) << var << " ";
        }
        cout << endl;

        for (auto i = 0; i < (width + 1) * c.vars.size() - 1; ++i) {
                cout << "-";
        }
        cout << endl;

        for (auto row : c.rows) {
                for (auto j = 0; j < c.vars.size(); ++j) {
                        cout << setw(width) << ALPHABET[row.a[j]] << " ";
                }
                cout << "= " << row.v << endl;
        }
}

void cost_dot(cost const &c, const char *root_dir) {

        char cwd[PATH_MAX];
        getcwd(cwd, sizeof(cwd));
        chdir(root_dir);
        ostringstream oss;
        oss << c.vars[0];

        for (auto i = 1; i < c.vars.size(); ++i) {
                oss << "-" << c.vars[i];
        }

        mkdir(oss.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        chdir(oss.str().c_str());

        for (auto row : c.rows) {
                ostringstream filename;
                filename << row.v << ".dot";
                auto fp = fopen(filename.str().c_str(), "w");
                fa_dot(fp, row.fa);
                fclose(fp);
        }

        chdir(cwd);
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

        const auto n_rows = accumulate(c.domains.begin(), c.domains.end(), 1, std::multiplies<size_t>());

        for (size_t i = 0; i < n_rows; ++i) {
                row r = {
                        vector<size_t>(c.vars.size()),
                        fa_make_basic(FA_EMPTY),
                        def
                };
                auto x = i;
                for (size_t j = 0; j < c.vars.size(); ++j) {
                        r.a[j] = x % c.domains[j];
                        x /= c.domains[j];
                }
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
vector<T> tokenize(ifstream &f, const char *sep = " ") {

        string str;
        getline(f, str);
        char *dup = strdup(str.c_str());
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

template <typename T, size_t SKIP, size_t MAX_N>
array<T, MAX_N + 1 - SKIP> tokenize(ifstream &f, const char *sep = " ") {

        string str;
        getline(f, str);
        char *dup = strdup(str.c_str());
        char *token = strtok(dup, sep);
        array<T, MAX_N + 1 - SKIP> a;
        auto skip = SKIP;
        auto max_n = MAX_N;
        auto i = 0;

        while (token != NULL) {
                if (skip) {
                        skip--;
                } else {
                        if (max_n) {
                                if constexpr (is_integral_v<T>) {
                                        a[i++] = atoi(token);
                                } else if (is_floating_point_v<T>) {
                                        a[i++] = atof(token);
                                }
                                max_n--;
                        } else {
                                break;
                        }
                }
                token = strtok(NULL, sep);
        }

        free(dup);
        return a;
}

vector<boost::dynamic_bitset<>> read_adj(const char *wcsp) {

        ifstream f(wcsp);
        const auto [ n_vars, max_domain, n_costs ] = tokenize<size_t, 1, 3>(f);
        f.ignore(numeric_limits<streamsize>::max(), '\n');
        vector<boost::dynamic_bitset<>> adj(n_vars);

        for (auto i = 0; i < n_vars; ++i) {
                adj[i] = boost::dynamic_bitset<>(n_vars);
        }

        for (auto i = 0; i < n_costs; ++i) {
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

vector<cost> read_costs(const char *wcsp, value threshold) {

        ifstream f(wcsp);
        const auto [ n_vars, max_domain, n_costs ] = tokenize<size_t, 1, 3>(f);
        const auto all_domains = tokenize<size_t>(f);
        print_it(all_domains);
        vector<cost> costs(n_costs);

        for (auto i = 0; i < n_costs; ++i) {

                cost c;
                auto temp = tokenize<value>(f);

                for (auto it = temp.begin() + 1; it != temp.begin() + temp[0] + 1; ++it) {
                        c.vars.push_back(*it);
                        c.domains.push_back(all_domains[*it]);
                }

                vector<size_t> pfx_prod(c.domains.size());
                exclusive_scan(c.domains.begin(), c.domains.end(), pfx_prod.begin(), 1, multiplies<>{});
                preallocate_rows(c, temp[temp[0] + 1]);

                for (auto j = 0; j < temp[temp[0] + 2]; ++j) {
                        auto row = tokenize<value>(f);
                        const value val = row[row.size() - 1];
                        auto idx = 0;
                        for (auto k = 0; k < row.size() - 1; ++k) {
                                idx += row[k] * pfx_prod[k];
                        }
                        c.rows[idx].v = val;
                }

                // sort according to values
                sort(c.rows.begin(), c.rows.end(), [](const row &x, const row &y) { return (x.v < y.v); });
                remove_threshold(c, threshold);
                costs[i] = c;
        }

        f.close();
        return costs;
}
