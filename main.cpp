#include "main.hpp"

bool parallel = false;

static inline void print_usage(const char *bin) {

        cout << "Usage: " << bin << " [-h] [-a bub|brz|hop] [-i bound] [-o wmf|mf|miw|md|rand|*.pt] -f instance" << endl;
}

static inline bool exists(const char *filename) {

        FILE *file = fopen(filename, "r");
        if (!file) {
                return false;
        } else {
                fclose(file);
                return true;
        }
}

int main(int argc, char *argv[]) {

        fa_minimization_algorithm = FA_MIN_BUBENZER;
        size_t ibound = 0;
        char *instance = NULL;
        int order_heur = 0;
        char *pseudotree = NULL;
        int opt;

        while ((opt = getopt(argc, argv, "a:i:f:o:h")) != -1) {
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
                        case 'i':
                                ibound = max(0, atoi(optarg));
                                continue;
                        case 'f':
                                if (exists(optarg)) {
                                        instance = optarg;
                                } else {
                                        cerr << argv[0] << ": file not found -- '" << optarg << "'" << endl;
                                        print_usage(argv[0]);
                                        return EXIT_FAILURE;
                                }
                                continue;
                        case 'o':
                                if (strcmp(optarg, "wmf") == 0) {
                                        order_heur = WEIGHTED_MIN_FILL;
                                } else if (strcmp(optarg, "mf") == 0) {
                                        order_heur = MIN_FILL;
                                } else if (strcmp(optarg, "miw") == 0) {
                                        order_heur = MIN_INDUCED_WIDTH;
                                } else if (strcmp(optarg, "md") == 0) {
                                        order_heur = MIN_DEGREE;
                                } else if (strcmp(optarg, "rand") == 0) {
                                        order_heur = RANDOM;
                                } else if (exists(optarg)) {
                                        pseudotree = optarg;
                                } else {
                                        cerr << argv[0] << ": file not found -- '" << optarg << "'" << endl;
                                        print_usage(argv[0]);
                                        return EXIT_FAILURE;
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
                inst_type = WCSP;
        } else {
                cout << "Found MPE instance: executing MAX-PROD algorithm" << endl;
                inst_type = MPE;
        }

        if (ibound) {
                cout << "I-bound = " << ibound << endl;
        }

        if (parallel) {
                cout << "Parallelism enabled" << endl;
        }

        string algorithms[] = { "Hopcroft", "Brzozowski", "Bubenzer" };
        cout << "Minimisation algorithm = " << algorithms[fa_minimization_algorithm] << endl;

        // look for a known threshold to remove rows
        value threshold = numeric_limits<value>::max();

        for (size_t i = 0; i < N_DATASETS; ++i) {
                if (strstr(instance, datasets[i]) != NULL) {
                        threshold = thresholds[i];
                        cout << "Thresholds value = " << threshold << endl;
                }
        }

        auto [ domains, adj, weights ] = read_domains_adj_weights(instance, inst_type);
        //print_adj(adj);
        //cout << endl;
        //cout << vec2str(weights, "Weights") << endl;

        string heuristics[] = { "WEIGHTED-MIN-FILL", "MIN-FILL", "MIN-INDUCED-WIDTH", "MIN-DEGREE" };
        vector<size_t> order;

        if (pseudotree) {
                cout << "Reading order from " << pseudotree << endl;
                order = read_pseudotree_order(pseudotree, domains);
        } else {
                if (order_heur == RANDOM) {
                        cout << "Computing RANDOM variable order..." << endl;
                        order.resize(domains.size());
                        iota(order.begin(), order.end(), 0);
                        srand(unsigned (std::time(0)));
                        random_shuffle(order.begin(), order.end());
                } else {
                        cout << "Computing " << heuristics[order_heur] << " variable order..." << endl;
                        order = greedy_order(order_heur, adj, weights);
                }
        }

        cout << vec2str(order, "Order") << endl;
        reverse(order.begin(), order.end());
        vector<size_t> pos(order.size());
        for (size_t i = 0; i < order.size(); ++i) {
                pos[order[i]] = i;
        }

        cout << "Induced width = " << induced_width(adj, order, pos) << endl;
        auto tables = read_tables(instance, inst_type, pos, threshold);

        /*for (auto const &table : tables) {
                print_table(table);
                cout << endl;
        }*/

        auto start_t = chrono::high_resolution_clock::now();
        vector<automata> automatas(tables.size());
        double total_rows = 0;
        double actual_rows = 0;
        cout << "Computing automata..." << endl;

        #pragma omp parallel for schedule(dynamic) if (parallel)
        for (size_t i = 0; i < tables.size(); ++i) {
                automatas[i] = compute_automata(tables[i]);
                actual_rows += automatas[i].rows.size();
                total_rows += accumulate(automatas[i].domains.begin(),
                                         automatas[i].domains.end(),
                                         1, multiplies<size_t>());
        }

        cout << "Redundancy = " << 1 - actual_rows / total_rows << endl << endl;

        /*for (auto const &a : automatas) {
                automata_dot(a, "dot");
        }*/

        auto buckets = compute_buckets(automatas, pos);

        /*for (auto i : order) {
                cout << "Bucket " << i << endl << endl;
                for (auto a : buckets[i]) {
                        print_table(compute_table(a));
                        cout << endl;
                }
        }*/

        //const int inner = (inst_type == WCSP) ? BE_SUM : BE_PROD;
        //const int outer = (inst_type == WCSP) ? BE_MIN : BE_MAX;
        const auto optimal = bucket_elimination(buckets, BE_SUM, BE_MIN, order, pos, domains, ibound);

        chrono::duration<double> runtime = chrono::high_resolution_clock::now() - start_t;
        cout << endl << "Time elapsed = " << runtime.count() << endl;
        cout << optimal << endl;

        return EXIT_SUCCESS;
}
