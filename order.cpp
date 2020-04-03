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
