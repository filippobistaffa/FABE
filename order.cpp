#include "order.hpp"

#include <cmath>            // fabs
#include <iomanip>          // ostringstream
#include <numeric>          // accumulate
#include <algorithm>        // count_if
#include <fstream>          // ifstream, getline
#include <unordered_set>    // unordered_set

// fmt library
#define FMT_HEADER_ONLY
#include <fmt/core.h>
#include <fmt/ranges.h>

#include "io.hpp"

#define GET_MACRO(_1,_2,_3,NAME,...) NAME
#define EACH_NONZERO(...) GET_MACRO(__VA_ARGS__, EACH_NONZERO_3, EACH_NONZERO_2)(__VA_ARGS__)
#define EACH_NONZERO_2(B, I) (size_t I = 0; I < (B).size(); ++I) if ((B)[I] != 0)
#define EACH_NONZERO_3(B, I, S) (size_t I = S; I < (B).size(); ++I) if ((B)[I] != 0)
#define AVOID_ZERO(X) ((X) == 0 ? 1 : (X))

//#define DEBUG_GREEDY_ORDER

/*static inline weight new_edge_value(std::vector<std::vector<weight>> const &adj, size_t i, size_t j) {

        const weight avg1 = accumulate(adj[i].begin(), adj[i].end(), 0.0) / adj[i].size();
        const weight avg2 = accumulate(adj[j].begin(), adj[j].end(), 0.0) / adj[j].size();
        return (avg1 + avg2) / 2;
}*/

#define EDGE_VAL(AVG_W, I, J) ((AVG_W[I] + AVG_W[J]) / 2)

static inline weight metric(int order_heur, std::vector<std::vector<weight>> const &adj, size_t node, std::vector<weight> avg_w) {

        if (order_heur == O_MIN_DEGREE || order_heur == O_MIN_INDUCED_WIDTH) {
                return avg_w[node];
        } else {
                weight fill = 0;
                for EACH_NONZERO(adj[node], i) {
                        for EACH_NONZERO(adj[node], j, i + 1) {
                                if (!adj[i][j]) {
                                        if (order_heur == O_WEIGHTED_MIN_FILL) {
                                                #ifdef DEBUG_GREEDY_ORDER
                                                fmt::print("edge ({}, {}) not present, adding {}\n", i, j, EDGE_VAL(avg_w, i, j));
                                                #endif
                                                fill += EDGE_VAL(avg_w, i, j);
                                        } else {
                                                #ifdef DEBUG_GREEDY_ORDER
                                                fmt::print("edge ({}, {}) not present, adding 1\n", i, j, EDGE_VAL(avg_w, i, j));
                                                #endif
                                                fill++;
                                        }
                                }
                        }
                }
                return fill;
        }
}

static inline void connect_neighbours(int order_heur, std::vector<std::vector<weight>> &adj, size_t node, std::vector<weight> avg_w) {

        for EACH_NONZERO(adj[node], i) {
                for EACH_NONZERO(adj[node], j, i + 1) {
                        if (!adj[i][j]) {
                                if (order_heur == O_WEIGHTED_MIN_FILL) {
                                        adj[i][j] = adj[j][i] = EDGE_VAL(avg_w, i, j);
                                } else {
                                        adj[i][j] = adj[j][i] = 1;
                                }
                                #ifdef DEBUG_GREEDY_ORDER
                                fmt::print("edge ({}, {}) <- {}\n", i, j, adj[i][j]);
                                fmt::print("edge ({}, {}) <- {}\n", j, i, adj[i][j]);
                                #endif
                        }
                }
        }
}

std::vector<size_t> greedy_order(std::vector<std::vector<weight>> const &adj, int order_heur, int tie_heur) {

        std::vector<size_t> order;
        std::vector<std::vector<weight>> tmp_adj(adj);
        std::unordered_set<size_t> not_marked(adj.size());
        std::vector<weight> avg_w(adj.size());

        for (size_t i = 0; i < adj.size(); ++i) {
                not_marked.insert(i);
                avg_w[i] = accumulate(tmp_adj[i].begin(), tmp_adj[i].end(), 0.0) / tmp_adj.size();
        }

        //std::vector<weight> avg_w_in = avg_w;

        while (!not_marked.empty()) {
                #ifdef DEBUG_GREEDY_ORDER
                fmt::print("avg_w: {}\n", avg_w);
                #endif
                std::vector<size_t> cand;
                weight min_met = std::numeric_limits<weight>::max();
                for (auto i : not_marked) {
                        weight met = metric(order_heur, tmp_adj, i, avg_w);
                        #ifdef DEBUG_GREEDY_ORDER
                        fmt::print("metric({}) = {} (min = {})\n", i, met, min_met);
                        #endif
                        if (fabs(met - min_met) <= std::numeric_limits<weight>::epsilon()) {
                                cand.push_back(i);
                        } else if (met < min_met) {
                                min_met = met;
                                cand.clear();
                                cand.push_back(i);
                        }
                }
                #ifdef DEBUG_GREEDY_ORDER
                fmt::print("candidates: {}\n", cand);
                for (auto c : cand) {
                        fmt::print("{} ", avg_w[c]);
                }
                fmt::print("\n");
                #endif
                size_t sel;
                if (tie_heur == T_RANDOM) {
                        sel = cand[rand() % cand.size()];
                } else {
                        sel = *min_element(cand.begin(), cand.end(), compare_vec(avg_w));
                }
                #ifdef DEBUG_GREEDY_ORDER
                fmt::print("selected = {}\n", sel);
                #endif
                // connect neighbours
                if (order_heur != O_MIN_DEGREE) {
                        connect_neighbours(order_heur, tmp_adj, sel, avg_w);
                }
                // remove node from graph
                not_marked.erase(sel);
                for (size_t i = 0; i < adj.size(); ++i) {
                        if (tmp_adj[i][sel]) {
                                tmp_adj[i][sel] = 0;
                                avg_w[i] = accumulate(tmp_adj[i].begin(), tmp_adj[i].end(), 0.0) / tmp_adj.size();
                        }
                }
                order.push_back(sel);
        }

        return order;
}

size_t induced_width(std::vector<std::vector<weight>> const &adj, std::vector<size_t> const &order) {

        std::vector<std::vector<weight>> tmp_adj(adj);
        size_t w = 0;

        for (auto i = order.rbegin(); i != order.rend(); ++i) {
                const size_t deg = count_if(tmp_adj[*i].begin(), tmp_adj[*i].end(), [](weight i) { return i > 0; });
                w = std::max(w, deg);
                for EACH_NONZERO(tmp_adj[*i], n1) {
                        tmp_adj[n1][*i] = 0;
                        for EACH_NONZERO(tmp_adj[*i], n2, n1 + 1) {
                                tmp_adj[n2][*i] = 0;
                                if (!tmp_adj[n1][n2]) {
                                        tmp_adj[n1][n2] = 1;
                                        tmp_adj[n2][n1] = 1;
                                }
                        }
                }
        }

	return w;
}

static inline void parse(std::string str, std::vector<size_t> &order) {

        size_t var;
        size_t i = str.find("(");

        if (i != std::string::npos) {
                var = atoi(str.substr(0, i).c_str());
                auto children = str.substr(i);
                while (children.size() > 0) {
                        size_t j = 1;
                        size_t b = 1;
                        while (b) {
                                if (children[j] == '(') {
                                        b++;
                                } else if (children[j] == ')') {
                                        b--;
                                }
                                j++;
                        }
                        parse(children.substr(1, j - 2), order);
                        children = children.substr(j);
                }
        } else {
                var = atoi(str.c_str());
        }

        order.push_back(var);
}

std::vector<size_t> read_pseudotree_order(const char *filename, std::vector<size_t> const &domains) {

        std::vector<size_t> vars;
        std::vector<size_t> ev(domains.size());

        for (size_t i = 0; i < domains.size(); ++i) {
                if (domains[i] == 1) {
                        vars.push_back(i);
                        ev[i] = 1;
                }
        }

        std::vector<size_t> offset(ev.size());
        partial_sum(ev.begin(), ev.end(), offset.begin());

        for (int i = offset.size(); i --> 0; ) {
                if (ev[i]) {
                        offset.erase(offset.begin() + i);
                }
        }

        std::vector<size_t> order;
        std::ifstream f(filename);
        std::string str;
        getline(f, str);
        parse(str.substr(1, str.size() - 2), order);
        order.pop_back();
        f.close();

        for (size_t i = 0; i < order.size(); ++i) {
                order[i] += offset[order[i]];
        }

        order.insert(order.begin(), vars.begin(), vars.end());

        return order;
}

void export_order(std::vector<size_t> const &order, std::vector<size_t> const &domains, const char *output) {

        std::vector<size_t> ev(order.size());
        size_t n = 0;

        for (auto var : order) {
                if (domains[var] > 1) {
                        n++;
                } else {
                        ev[var] = 1;
                }
        }

        std::vector<size_t> rem(order.size());
        partial_sum(ev.begin(), ev.end(), rem.begin());

        std::ostringstream oss;
        oss << "# exported in aolib format: first line is the number of variables, ";
        oss << "variables with domain 1 are considered evidence and removed, ";
        oss << "remaining variables are re-indexed." << '\n';
        oss << n << '\n';

        for (auto it = order.rbegin(); it != order.rend(); ++it) {
                if (domains[*it] > 1) {
                        oss << *it - rem[*it] << '\n';
                }
        }

        std::ofstream f(output);
        f << oss.str();
        f.close();
}
