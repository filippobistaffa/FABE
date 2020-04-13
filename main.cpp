#include "main.hpp"

#ifdef PRINT_VAR_POS
extern vector<size_t> var_map;
#endif

bool parallel;

static void print_usage(const char *bin) {

        cout << "Usage: " << bin << " [-h] [-p] [-a bub|brz|hop] [-i bound] -f instance" << endl;
}

int main(int argc, char *argv[]) {

        fa_minimization_algorithm = FA_MIN_BUBENZER;
        size_t ibound = 0;
        char *instance = NULL;
        int opt;

        while ((opt = getopt(argc, argv, "a:i:f:ph")) != -1) {
                switch (opt) {
                        case 'a':
                                if (strcmp(optarg, "bub") == 0) {
                                        fa_minimization_algorithm = FA_MIN_BUBENZER;
                                } else if (strcmp(optarg, "brz") == 0) {
                                        fa_minimization_algorithm = FA_MIN_BRZOZOWSKI;
                                } else if (strcmp(optarg, "hop") == 0) {
                                        fa_minimization_algorithm = FA_MIN_HOPCROFT;
                                } else {
                                        cerr << argv[0] << ": invalid minimisation algorithm -- '" << optarg << "'" << endl;
                                        print_usage(argv[0]);
                                        return EXIT_FAILURE;
                                }
                                continue;
                        case 'p':
                                parallel = true;
                                continue;
                        case 'i':
                                ibound = max(0, atoi(optarg));
                                continue;
                        case 'f':
                                FILE *file;
                                if (!(file = fopen(optarg, "r"))) {
                                        cerr << argv[0] << ": file not found -- '" << optarg << "'" << endl;
                                        print_usage(argv[0]);
                                        return EXIT_FAILURE;
                                } else {
                                        instance = optarg;
                                        fclose(file);
                                }
                                continue;
                        case 'h':
                        default :
                                print_usage(argv[0]);
                                return EXIT_FAILURE;
                }
        }

        if (!instance) {
                cerr << argv[0] << ": instance not specified!" << endl;
                print_usage(argv[0]);
                return EXIT_FAILURE;
        }

        int inst_type;

        if (strstr(instance, "wcsp")) {
                cout << "Found WCSP instance: executing MIN-SUM algorithm" << endl;
                inst_type = BE_WCSP;
        } else {
                cout << "Found MPE instance: executing MAX-PROD algorithm" << endl;
                inst_type = BE_MPE;
        }

        if (ibound) {
                cout << "I-bound = " << ibound << endl;
        }

        if (parallel) {
                cout << "Parallelism enabled" << endl;
        }

        string algorithms[] = { "Hopcroft", "Brzozowski", "Bubenzer" };
        cout << "Minimisation algorithm = " << algorithms[fa_minimization_algorithm] << endl;
        auto start_t = chrono::high_resolution_clock::now();

        // look for a known threshold to remove rows
        value threshold = numeric_limits<value>::max();

        for (auto i = 0; i < N_DATASETS; ++i) {
                if (strstr(instance, datasets[i]) != NULL) {
                        threshold = thresholds[i];
                        cout << "Thresholds value = " << threshold << endl;
                }
        }

        auto adj = read_adj(instance);
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

        auto [ domains, tables ] = read_domains_tables(instance, pos, threshold);

        /*for (auto const &table : tables) {
                print_table(table);
                cout << endl;
        }*/

        vector<automata> automatas(tables.size());
        double total_rows = 0;
        double actual_rows = 0;
        cout << "Computing automata..." << endl;

        #pragma omp parallel for schedule(dynamic) if (parallel)
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

        const int inner = (inst_type == BE_WCSP) ? BE_SUM : BE_PROD;
        const int outer = (inst_type == BE_WCSP) ? BE_MIN : BE_MAX;
        const auto optimal = bucket_elimination(buckets, inner, outer, order, pos, domains, ibound);

        chrono::duration<double> runtime = chrono::high_resolution_clock::now() - start_t;
        cout << endl << "Time elapsed = " << runtime.count() << endl;
        cout << optimal << endl;

        return EXIT_SUCCESS;
}
