#include "main.hpp"

int main(int argc, char *argv[]) {

        if (argc != 2) {
                cout << "Usage: " << argv[0] << " wcsp_instance" << endl;
                return EXIT_FAILURE;
        }

        auto adj = read_adj(argv[1]);
        print_adj(adj);
        cout << endl;

        auto order = greedy_order(adj);
        vector<size_t> pos(order.size());

        for (auto i = 0; i < order.size(); ++i) {
                pos[order[i]] = i;
        }

        auto costs = read_costs(argv[1]);
 
        for (auto cost : costs) {
                print_cost_table(cost);
                cout << endl;
        }

        vector<cost> fas(costs.size());

        #pragma omp parallel for
        for (auto i = 0; i < costs.size(); ++i) {
                fas[i] = compress_clusters(costs[i]);
        }

        for (auto cost : fas) {
                cost_dot(cost, "dot");
        }

        //cost_dot(fas[0], "dot");

        //reduce_var(fas[0], 2);

        /*auto buckets = compute_buckets(costs, pos);

        for (auto i : order) {
                cout << "Bucket " << i << " ";
                print_it(bin_vars[i]);
                cout << endl;
                for (auto c : buckets[i]) {
                        print_cost_table(c);
                        cout << endl;
                }
        }

        print_it(order, "Ord.");
        print_it(pos, "Pos.");
        cout << "I.W. = " << induced_width(adj, order, pos) << endl;*/

        for (auto cost : costs) {
                for (auto row : cost.rows) {
                        fa_free(row.fa);
                }
        }

        return EXIT_SUCCESS;
}
