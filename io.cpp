#include "io.hpp"

void print_table(table const &t) {

        const auto width = max(1.0, 1 + floor(log10(*max_element(t.vars.begin(), t.vars.end()))));

        for (auto var : t.vars) {
                cout << setw(width) << var << " ";
        }
        cout << endl;

        for (auto i = 0; i < (width + 1) * t.vars.size() - 1; ++i) {
                cout << "-";
        }
        cout << endl;

        for (auto row : t.rows) {
                for (auto j = 0; j < t.vars.size(); ++j) {
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

        for (auto& [ v, fa ] : a.rows) {
                fa_make_dot(fa, "%.0f.dot", v);
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
void preallocate_rows(table &t, value def) {

        const auto n_rows = accumulate(t.domains.begin(), t.domains.end(), 1, multiplies<size_t>());
        t.rows.resize(n_rows);

        for (size_t i = 0; i < n_rows; ++i) {
                t.rows[i] = make_pair(get_combination(i, t.domains), def);
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
        const auto [ n_vars, max_domain, n_tables ] = tokenize<size_t, 1, 3>(f);
        f.ignore(numeric_limits<streamsize>::max(), '\n');
        vector<boost::dynamic_bitset<>> adj(n_vars);

        for (auto i = 0; i < n_vars; ++i) {
                adj[i] = boost::dynamic_bitset<>(n_vars);
        }

        for (auto i = 0; i < n_tables; ++i) {
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

pair<vector<size_t>, vector<table>> read_domains_tables(const char *wcsp, value threshold) {

        ifstream f(wcsp);
        const auto [ n_vars, max_domain, n_tables ] = tokenize<size_t, 1, 3>(f);
        const auto all_domains = tokenize<size_t>(f);
        vector<table> tables(n_tables);

        for (auto i = 0; i < n_tables; ++i) {

                table t;
                auto temp = tokenize<value>(f);

                for (auto it = temp.begin() + 1; it != temp.begin() + temp[0] + 1; ++it) {
                        t.vars.push_back(*it);
                        t.domains.push_back(all_domains[*it]);
                }

                vector<size_t> pfx_prod(t.domains.size());
                exclusive_scan(t.domains.begin(), t.domains.end(), pfx_prod.begin(), 1, multiplies<>{});
                preallocate_rows(t, temp[temp[0] + 1]);

                for (auto j = 0; j < temp[temp[0] + 2]; ++j) {
                        auto row = tokenize<value>(f);
                        const value val = row[row.size() - 1];
                        auto idx = 0;
                        for (auto k = 0; k < row.size() - 1; ++k) {
                                idx += row[k] * pfx_prod[k];
                        }
                        t.rows[idx].second = val;
                }

                // sort according to values
                sort(t.rows.begin(), t.rows.end(), [](pair<vector<size_t>, value> const &x,
                                                      pair<vector<size_t>, value> const &y)
                                                      { return (x.second < y.second); });
                remove_threshold(t, threshold);
                tables[i] = t;
        }

        f.close();
        return make_pair(all_domains, tables);
}
