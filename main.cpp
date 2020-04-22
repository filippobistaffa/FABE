#include "main.hpp"

#include <iostream>     // cout
#include <chrono>       // time measurement
#include <numeric>      // accumulate
#include <string.h>     // strcmp
#include <unistd.h>     // getopt
#include <algorithm>    // reverse, random_shuffle

#include "be.hpp"
#include "io.hpp"
#include "log.hpp"
#include "order.hpp"
#include "conversion.hpp"

bool parallel = false;

static inline void print_usage(const char *bin) {

        cout << "Usage: " << bin << " [-h] [-a bub|brz|hop] [-i bound] [-o wmf|mf|miw|md|random|*.pt] [-t unique|random] [-s seed] -f instance" << endl;
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
        int ord_heur = O_WEIGHTED_MIN_FILL;
        int tie_heur = T_UNIQUENESS;
        char *instance = NULL;
        char *pseudotree = NULL;
        size_t seed = time(NULL);
        int opt;

        while ((opt = getopt(argc, argv, "a:i:f:o:t:s:h")) != -1) {
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
                                        ord_heur = O_WEIGHTED_MIN_FILL;
                                } else if (strcmp(optarg, "mf") == 0) {
                                        ord_heur = O_MIN_FILL;
                                } else if (strcmp(optarg, "miw") == 0) {
                                        ord_heur = O_MIN_INDUCED_WIDTH;
                                } else if (strcmp(optarg, "md") == 0) {
                                        ord_heur = O_MIN_DEGREE;
                                } else if (strcmp(optarg, "random") == 0) {
                                        ord_heur = O_RANDOM;
                                } else if (exists(optarg)) {
                                        pseudotree = optarg;
                                } else {
                                        cerr << argv[0] << ": file not found -- '" << optarg << "'" << endl;
                                        print_usage(argv[0]);
                                        return EXIT_FAILURE;
                                }
                                continue;
                        case 't':
                                if (strcmp(optarg, "unique") == 0) {
                                        tie_heur = T_UNIQUENESS;
                                } else if (strcmp(optarg, "random") == 0) {
                                        tie_heur = T_RANDOM;
                                } else {
                                        cerr << argv[0] << ": tie-breaking heuristic not valid -- '" << optarg << "'" << endl;
                                        print_usage(argv[0]);
                                        return EXIT_FAILURE;
                                }
                                continue;
                        case 's':
                                seed = atoi(optarg);
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

        log_line();
        log_value("Instance", instance);
        int inst_type;

        if (strstr(instance, "wcsp")) {
                log_value("Instance type", "WCSP");
                log_value("Bucket elimination algorithm", "MIN-SUM");
                inst_type = WCSP;
        } else {
                log_value("Instance type", "UAI");
                log_value("Bucket elimination algorithm", "MIN-SUM with -log()");
                inst_type = MPE;
        }

        log_value("I-bound", (ibound == 0) ? "inf" : to_string(ibound));
        log_value("Parallel mode", parallel);

        string algorithms[] = { "Hopcroft", "Brzozowski", "Bubenzer" };
        log_value("Automata minimisation algorithm", algorithms[fa_minimization_algorithm]);

        // look for a known threshold to remove rows
        value threshold = numeric_limits<value>::max();

        for (size_t i = 0; i < N_DATASETS; ++i) {
                if (strstr(instance, datasets[i]) != NULL) {
                        threshold = thresholds[i];
                }
        }

        log_value("Thresholds value", threshold);

        auto [ domains, adj ] = read_domains_adj(instance, inst_type);
        //print_adj(adj);
        //cout << endl;

        log_value("Seed", seed);
        string ord_heur_names[] = { "WEIGHTED-MIN-FILL", "MIN-FILL", "MIN-INDUCED-WIDTH", "MIN-DEGREE" };
        string tie_heur_names[] = { "MIN-UNIQUENESS", "RANDOM" };
        vector<size_t> order;
        auto start_t = chrono::high_resolution_clock::now();

        if (pseudotree) {
                log_value("Variable order heuristic", pseudotree);
                order = read_pseudotree_order(pseudotree, domains);
        } else {
                srand(seed);
                if (ord_heur == O_RANDOM) {
                        log_value("Variable order heuristic", "Random");
                        order.resize(domains.size());
                        iota(order.begin(), order.end(), 0);
                        srand(unsigned (std::time(0)));
                        random_shuffle(order.begin(), order.end());
                } else {
                        log_value("Variable order heuristic", ord_heur_names[ord_heur]);
                        log_value("Tie-breaking heuristic", tie_heur_names[tie_heur]);
                        order = greedy_order(adj, ord_heur, tie_heur);
                }
        }

        //cout << vec2str(order, "Order") << endl;
        reverse(order.begin(), order.end());
        vector<size_t> pos(order.size());
        for (size_t i = 0; i < order.size(); ++i) {
                pos[order[i]] = i;
        }

        chrono::duration<double> runtime = chrono::high_resolution_clock::now() - start_t;
        log_value("Order computation runtime", runtime.count());
        start_t = chrono::high_resolution_clock::now();
        //export_order(order, domains, "order.vo");
        log_value("Induced width", induced_width(adj, order));
        auto tables = read_tables(instance, inst_type, pos, threshold);

        /*for (auto const &table : tables) {
                print_table(table);
                cout << endl;
        }*/

        vector<automata> automatas(tables.size());
        double total_rows = 0;
        double actual_rows = 0;

        #pragma omp parallel for schedule(dynamic) if (parallel)
        for (size_t i = 0; i < tables.size(); ++i) {
                automatas[i] = compute_automata(tables[i]);
                actual_rows += automatas[i].rows.size();
                total_rows += accumulate(automatas[i].domains.begin(),
                                         automatas[i].domains.end(),
                                         1, multiplies<size_t>());
        }

        log_value("Value redundancy", 1 - actual_rows / total_rows);

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

        log_line();

        //const int inner = (inst_type == WCSP) ? BE_SUM : BE_PROD;
        //const int outer = (inst_type == WCSP) ? BE_MIN : BE_MAX;
        const auto optimal = bucket_elimination(buckets, BE_SUM, BE_MIN, order, pos, domains, ibound);

        runtime = chrono::high_resolution_clock::now() - start_t;
        log_value("Solution", optimal);
        log_value("Bucket elimination runtime", runtime.count());
        log_line();

        return EXIT_SUCCESS;
}
