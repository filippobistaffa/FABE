#include "order.hpp"

#ifdef MINWIDTH
#define METRIC degree
#endif

#ifdef MININDUCEDWIDTH
#define METRIC degree
#define CONNECT
#endif

#ifdef MINFILL
#define METRIC fill
#define CONNECT
#endif

#ifdef WEIGHTEDMINFILL
#define METRIC fill
#define CONNECT
#endif

__attribute__((always_inline)) inline
auto degree(size_t row, vector<boost::dynamic_bitset<>> const &adj, boost::dynamic_bitset<> const &mask, vector<float> const &weights) {

        return (adj[row] & mask).count();
}

__attribute__((always_inline)) inline
auto fill(size_t row, vector<boost::dynamic_bitset<>> const &adj, boost::dynamic_bitset<> const &mask, vector<float> const &weights) {

        float fill = 0;
        auto tmp = adj[row] & mask;

        for EACH_SET_BIT(tmp, i) {
                for EACH_SET_BIT(tmp, j, i) {
                        if (!(adj[i] & mask).test(j)) {
                                #ifdef WEIGHTEDMINFILL
                                fill += weights[i] + weights[j];
                                #else
                                fill++;
                                #endif
                        }
                }
        }

        return fill;
}

__attribute__((always_inline)) inline
void connect_neighbours(size_t row, vector<boost::dynamic_bitset<>> &adj, boost::dynamic_bitset<> const &mask) {

        auto tmp = adj[row] & mask;

        for EACH_SET_BIT(tmp, i) {
                for EACH_SET_BIT(tmp, j, i) {
                        adj[i].set(j);
                        adj[j].set(i);
                }
        }
}

vector<size_t> greedy_order(vector<boost::dynamic_bitset<>> const &adj, vector<float> const &weights) {

        vector<size_t> order;

        const auto n_bin_vars = adj.size();
        boost::dynamic_bitset<> not_assigned(n_bin_vars);
        not_assigned.set();

        /*for (auto i = 0; i < adj.size(); ++i) {
                if (METRIC(i, adj, not_assigned) == 0) {
                        cout << i << " -> " << 0 << endl;
                        order.push_back(i);
                        not_assigned.reset(i);
                }
        }*/

        vector<boost::dynamic_bitset<>> tmp_adj(adj);

        while (not_assigned.any()) {
                //cout << not_assigned << endl;
                size_t min_node;
                auto min_metric = numeric_limits<float>::max();
                for EACH_SET_BIT(not_assigned, i) {
                        auto metric = METRIC(i, tmp_adj, not_assigned, weights);
                        //cout << i << " -> " << metric << endl;
                        if (metric < min_metric) {
                                min_metric = metric;
                                min_node = i;
                        }
                }
                order.push_back(min_node);
                #ifdef CONNECT
                connect_neighbours(min_node, tmp_adj, not_assigned);
                #endif
                not_assigned.reset(min_node);
                //print_it(order);
        }

        return order;
}

size_t induced_width(vector<boost::dynamic_bitset<>> const &adj, vector<size_t> const &order, vector<size_t> const &pos) {

        vector<boost::dynamic_bitset<>> tmp_adj(adj);
	size_t w = 0;

        for (auto i = order.rbegin(); i != order.rend(); ++i) {
                w = max(w, tmp_adj[*i].count());
                for EACH_SET_BIT(tmp_adj[*i], n1) {
                        tmp_adj[n1].reset(*i);
                        for EACH_SET_BIT(tmp_adj[*i], n2, n1) {
                                tmp_adj[n2].reset(*i);
                                if (!tmp_adj[n1].test(n2)) {
                                        tmp_adj[n1].set(n2);
                                        tmp_adj[n2].set(n1);
                                }
                        }
                }
        }

	return w;
}

static inline void parse(string str, vector<size_t> &order) {

        size_t var;
        size_t i = str.find("(");

        if (i != string::npos) {
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

vector<size_t> read_pseudotree_order(const char *filename, vector<size_t> const &domains) {

        vector<size_t> vars;
        vector<size_t> ev(domains.size());

        for (size_t i = 0; i < domains.size(); ++i) {
                if (domains[i] == 1) {
                        vars.push_back(i);
                        ev[i] = 1;
                }
        }

        vector<size_t> offset(ev.size());
        partial_sum(ev.begin(), ev.end(), offset.begin());

        for (int i = offset.size(); i --> 0; ) {
                if (ev[i]) {
                        offset.erase(offset.begin() + i);
                }
        }

        vector<size_t> order;
        ifstream f(filename);
        string str;
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

/*void export_order(vector<size_t> const &order, vector<size_t> const &domains, const char *output) {

        vector<size_t> ev(order.size());
        size_t n = 0;

        for (auto var : order) {
                if (domains[var] > 1) {
                        n++;
                } else {
                        ev[var] = 1;
                }
        }

        vector<size_t> rem(order.size());
        partial_sum(ev.begin(), ev.end(), rem.begin());

        ostringstream oss;
        oss << "# exported in aolib format: first line is the number of variables, ";
        oss << "variables with domain 1 are considered evidence and removed, ";
        oss << "remaining variables are re-indexed." << endl;
        oss << n << endl;

        for (auto it = order.rbegin(); it != order.rend(); ++it) {
                if (domains[*it] > 1) {
                        oss << *it - rem[*it] << endl;
                }
        }

        ofstream f(output);
        f << oss.str();
        f.close();
}*/
