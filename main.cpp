#include "main.hpp"

#include <iostream>     // cout
#include <chrono>       // time measurement
#include <numeric>      // accumulate
#include <string.h>     // strcmp
#include <unistd.h>     // getopt
#include <algorithm>    // reverse, random_shuffle
#include <sstream>      // for ostringstream

#include "be.hpp"
#include "io.hpp"
#include "log.hpp"
#include "order.hpp"
#include "conversion.hpp"
#include "export.hpp"

// print tables after parsing
//#define PRINT_TABLES

// export dot representation of created automata
//#define EXPORT_AUTOMATA_DOT

bool parallel = false;

static inline void print_usage(const char *bin) {

        cerr << "Usage: " << bin << " [-h] [-a bub|brz|hop] [-i bound] [-o wmf|mf|miw|md|random|*.pt] ";
        cerr << "[-t unique|random] [-s seed] [-e tolerance] [-r] -f instance" << endl;
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
        value tolerance = 0;
        size_t ibound = 0;
        int ord_heur = O_WEIGHTED_MIN_FILL;
        int tie_heur = T_UNIQUENESS;
        char *instance = NULL;
        char *pseudotree = NULL;
        size_t seed = time(NULL);
        bool print_red = false;
        struct mbe_data mbe;
        int opt;

        while ((opt = getopt(argc, argv, "a:i:f:o:t:s:e:m:hr")) != -1) {
                switch (opt) {
                        case 'a':
                                if (strcmp(optarg, "bub") == 0) {
                                        fa_minimization_algorithm = FA_MIN_BUBENZER;
                                } else if (strcmp(optarg, "brz") == 0) {
                                        fa_minimization_algorithm = FA_MIN_BRZOZOWSKI;
                                } else if (strcmp(optarg, "hop") == 0) {
                                        fa_minimization_algorithm = FA_MIN_HOPCROFT;
                                } else {
                                        cerr << argv[0] << ": invalid minimisation algorithm -- '";
                                        cerr << optarg << "'" << endl;
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
                                        cerr << argv[0] << ": file not found -- '";
                                        cerr << optarg << "'" << endl;
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
                                        cerr << argv[0] << ": file not found -- '";
                                        cerr << optarg << "'" << endl;
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
                                        cerr << argv[0] << ": tie-breaking heuristic not valid -- '";
                                        cerr << optarg << "'" << endl;
                                        print_usage(argv[0]);
                                        return EXIT_FAILURE;
                                }
                                continue;
                        case 's':
                                seed = atoi(optarg);
                                continue;
                        case 'e':
                                tolerance = atof(optarg);
                                continue;
                        case 'm':
                                mbe.filename = optarg;
                                continue;
                        case 'r':
                                print_red = true;
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

        string algorithms[] = { "Hopcroft", "Brzozowski", "Bubenzer" };
        log_value("Automata minimisation algorithm", algorithms[fa_minimization_algorithm]);

        log_value("I-bound", (ibound == 0) ? "inf" : to_string(ibound));
        log_value("Tolerance", tolerance);
        //log_value("Parallel mode", parallel ? "Enabled" : "Disabled");

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
        vector<vector<size_t>> anc(domains.size(), vector<size_t>());
        auto start_t = chrono::high_resolution_clock::now();

        if (pseudotree) {
                log_value("Variable order heuristic", pseudotree);
                order = read_pseudotree_order(pseudotree, domains, anc);
                //for (size_t i = 0; i < anc.size(); ++i) {
                //        cout << i << " : " << vec2str(anc[i]) << endl;
                //}
        } else {
                srand(seed);
                if (ord_heur == O_RANDOM) {
                        log_value("Variable order heuristic", "RANDOM");
                        order.resize(domains.size());
                        iota(order.begin(), order.end(), 0);
                        random_shuffle(order.begin(), order.end());
                } else {
                        log_value("Variable order heuristic", ord_heur_names[ord_heur]);
                        log_value("Tie-breaking heuristic", tie_heur_names[tie_heur]);
                        order = greedy_order(adj, ord_heur, tie_heur);
                }
        }

        if (mbe.filename) {
                mbe.augmented = vector<vector<uchar>>(domains.size() + 1, vector<uchar>());
                mbe.n_augmented = vector<size_t>(domains.size());
                mbe.intermediate = vector<vector<uchar>>(domains.size() + 1, vector<uchar>());
                mbe.n_intermediate = vector<size_t>(domains.size());
                mbe.evid_offset = vector<size_t>(domains.size());
                mbe.anc = anc;
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

        #ifdef PRINT_TABLES
        cout << endl;
        for (auto const &table : tables) {
                print_table(table);
                cout << endl;
        }
        #endif

        vector<automata> automatas(tables.size());
        double total_rows = 0;
        double actual_rows = 0;
        value tot_error = 0;

        #pragma omp parallel for schedule(dynamic) if (parallel)
        for (size_t i = 0; i < tables.size(); ++i) {
                auto [ a, error ] = compute_automata(tables[i], tolerance);
                automatas[i] = a;
                tot_error += error;
                actual_rows += automatas[i].rows.size();
                total_rows += accumulate(automatas[i].domains.begin(),
                                         automatas[i].domains.end(),
                                         1, multiplies<size_t>());
        }

        ostringstream oss;
        oss << total_rows - actual_rows << "/" << total_rows << " (" << 1 - actual_rows / total_rows << ")";
        log_value("Value redundancy", oss.str());

        #ifdef EXPORT_AUTOMATA_DOT
        for (auto const &a : automatas) {
                automata_dot(a, "dot");
        }
        #endif

        auto buckets = compute_buckets(automatas, pos);

        /*for (auto i : order) {
                cout << "Bucket " << i << endl << endl;
                for (auto a : buckets[i]) {
                        print_table(compute_table(a));
                        cout << endl;
                }
        }*/

        log_line();

        if (print_red) {
                exit(0);
        }

        //const int inner = (inst_type == WCSP) ? BE_SUM : BE_PROD;
        //const int outer = (inst_type == WCSP) ? BE_MIN : BE_MAX;
        const auto optimal = bucket_elimination(buckets, BE_SUM, BE_MIN, order, pos, domains, ibound, mbe);
        runtime = chrono::high_resolution_clock::now() - start_t;
        log_value("Solution value", optimal);
        //log_value("Maximum optimality gap", tolerance * tables.size());
        //log_value("Maximum optimality gap (%)", 100 * (tolerance * tables.size()) / optimal);
        //log_value("Optimality gap", tot_error);
        oss.str(string());
        oss << tot_error << " (" << setprecision(3) << 100 * (tot_error) / optimal << "%)";
        log_value("Optimality gap", oss.str());
        log_value("Bucket elimination runtime", runtime.count());
        log_line();

        return EXIT_SUCCESS;
}
