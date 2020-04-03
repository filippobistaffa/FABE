#include "main.hpp"

#ifdef PRINT_VAR_POS
extern vector<size_t> var_map;
#endif

int main(int argc, char *argv[]) {

        if (argc != 2 && argc != 3) {
                cout << "Usage: " << argv[0] << " wcsp_instance [max_iter]" << endl;
                return EXIT_FAILURE;
        }

        size_t max_iter = numeric_limits<size_t>::max();

        if (argc == 3) {
                max_iter = max(1, atoi(argv[2]));
        }

        auto adj = read_adj(argv[1]);
        //print_adj(adj);
        //cout << endl;

        auto order = greedy_order(adj);
        reverse(order.begin(), order.end());
        vector<size_t> pos(order.size());

        for (auto i = 0; i < order.size(); ++i) {
                pos[order[i]] = i;
        }

        #ifdef PRINT_VAR_POS
        var_map = pos;
        #endif

        // look for a known threshold to remove rows
        value threshold = numeric_limits<value>::max();

        for (auto i = 0; i < N_DATASETS; ++i) {
                if (strstr(argv[1], datasets[i]) != NULL) {
                        threshold = thresholds[i];
                }
        }

        auto [ domains, tables ] = read_domains_tables(argv[1], pos, threshold);

        /*for (auto const &table : tables) {
                print_table(table);
                cout << endl;
        }*/

        vector<automata> automatas(tables.size());
        // change minimisation algorithm
        fa_minimization_algorithm = FA_MIN_BRZOZOWSKI;

        #pragma omp parallel for
        for (auto i = 0; i < tables.size(); ++i) {
                automatas[i] = compute_automata(tables[i]);
        }

        auto buckets = compute_buckets(automatas, pos);

        /*for (auto i : order) {
                cout << "Bucket " << i << endl << endl;
                for (auto a : buckets[i]) {
                        print_table(compute_table(a));
                        cout << endl;
                }
        }*/

        //cout << vec2str(order, "Ord.") << endl;
        //cout << vec2str(pos, "Pos.") << endl;
        cout << "I.W. = " << induced_width(adj, order, pos) << endl << endl;

        #ifdef PROFILE
        ProfilerStart(PROFILE);
        #endif

        const auto optimal = bucket_elimination(buckets, order, pos, domains, max_iter);

        #ifdef PROFILE
        ProfilerStop();
        #endif

        cout << endl << optimal << endl;

        return EXIT_SUCCESS;
}
