#include "main.hpp"

#ifdef PRINT_VAR_POS
extern vector<size_t> var_map;
#endif

int main(int argc, char *argv[]) {

        if (argc != 2 && argc != 3) {
                cout << "Usage: " << argv[0] << " wcsp_instance [i-bound]" << endl;
                return EXIT_FAILURE;
        }

        const size_t ibound = (argc == 3) ? max(0, atoi(argv[2])) : 0;
        auto start_t = chrono::high_resolution_clock::now();

        if (ibound) {
                cout << "I-bound = " << ibound << endl;
        }

        // look for a known threshold to remove rows
        value threshold = numeric_limits<value>::max();

        for (auto i = 0; i < N_DATASETS; ++i) {
                if (strstr(argv[1], datasets[i]) != NULL) {
                        threshold = thresholds[i];
                        cout << "Thresholds value = " << threshold << endl;
                }
        }

        auto adj = read_adj(argv[1]);
        //print_adj(adj);
        //cout << endl;

        cout << "Computing variable order..." << endl;
        auto order = greedy_order(adj);
        reverse(order.begin(), order.end());
        vector<size_t> pos(order.size());

        for (auto i = 0; i < order.size(); ++i) {
                pos[order[i]] = i;
        }

        #ifdef PRINT_VAR_POS
        var_map = pos;
        #endif

        auto [ domains, tables ] = read_domains_tables(argv[1], pos, threshold);

        /*for (auto const &table : tables) {
                print_table(table);
                cout << endl;
        }*/

        vector<automata> automatas(tables.size());
        // change minimisation algorithm
        //fa_minimization_algorithm = FA_MIN_BRZOZOWSKI;
        fa_minimization_algorithm = FA_MIN_BUBENZER;

        double total_rows = 0;
        double actual_rows = 0;

        //#pragma omp parallel for
        for (auto i = 0; i < tables.size(); ++i) {
                automatas[i] = compute_automata(tables[i]);
                actual_rows += automatas[i].rows.size();
                total_rows += accumulate(automatas[i].domains.begin(),
                                         automatas[i].domains.end(),
                                         1, multiplies<size_t>());
        }

        cout << "Redundancy = " << 1 - actual_rows / total_rows << endl << endl;

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
        //cout << "I.W. = " << induced_width(adj, order, pos) << endl << endl;

        const auto optimal = bucket_elimination(buckets, order, pos, domains, ibound);

        chrono::duration<double> runtime = chrono::high_resolution_clock::now() - start_t;
        cout << endl << "Time elapsed = " << runtime.count() << endl;
        cout << optimal << endl;

        return EXIT_SUCCESS;
}
