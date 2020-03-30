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

        auto [ domains, tables ] = read_domains_tables(argv[1]);
 
        for (auto table : tables) {
                print_table(table);
                cout << endl;
        }

        vector<automata> automatas(tables.size());

        #pragma omp parallel for
        for (auto i = 0; i < tables.size(); ++i) {
                automatas[i] = compute_automata(tables[i]);
        }

        auto t_buckets = compute_buckets(tables, pos);
        auto buckets = compute_buckets(automatas, pos);

        for (auto i : order) {
                cout << "Bucket " << i << endl << endl;
                for (auto t : t_buckets[i]) {
                        print_table(t);
                        cout << endl;
                }
        }

        cout << vec2str(order, "Ord.") << endl;
        cout << vec2str(pos, "Pos.") << endl;
        cout << "I.W. = " << induced_width(adj, order, pos) << endl;

        fa_minimization_algorithm = FA_MIN_BRZOZOWSKI;

        auto j = join_bucket(buckets[order.back()], domains);
        automata_dot(j, "dot");

        return EXIT_SUCCESS;
}
