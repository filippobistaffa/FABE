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

__attribute__((always_inline)) inline
auto degree(size_t row, vector<boost::dynamic_bitset<>> const &adj, boost::dynamic_bitset<> const &mask) {

        return (adj[row] & mask).count();
}

__attribute__((always_inline)) inline
auto fill(size_t row, vector<boost::dynamic_bitset<>> const &adj, boost::dynamic_bitset<> const &mask) {

        size_t fill = 0;
        auto tmp = adj[row] & mask;

        for EACH_SET_BIT(tmp, i) {
                for EACH_SET_BIT(tmp, j, i) {
                        if (!(adj[i] & mask).test(j)) {
                                fill++;
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

vector<size_t> greedy_order(vector<boost::dynamic_bitset<>> const &adj) {

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
                size_t min_node, min_metric = numeric_limits<size_t>::max();
                for EACH_SET_BIT(not_assigned, i) {
                        auto metric = METRIC(i, tmp_adj, not_assigned);
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

        //cout << "str = " << str << endl;
        size_t var;
        size_t i = str.find("(");
        //cout << "i = " << i << endl;
        if (i != string::npos) {
                var = atoi(str.substr(0, i).c_str());
                auto children = str.substr(i);
                //cout << children << endl;
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
                        //cout << "\"" << children.substr(1, j - 2) << "\"" <<endl;
                        //ret.ch.push_back(parse(children.substr(1, j - 2), order));
                        parse(children.substr(1, j - 2), order);
                        //cout << "\"" << children.substr(j) << "\"" <<endl;
                        children = children.substr(j);
                }
        } else {
                var = atoi(str.c_str());
        }
        //cout << ret.var << endl;
        order.push_back(var);
}

vector<size_t> read_pseudotree_order(const char *filename) {

        vector<size_t> order;
        ifstream f(filename);
        string str;
        getline(f, str);
        parse(str.substr(1, str.size() - 2), order);
        order.pop_back();
        f.close();
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
