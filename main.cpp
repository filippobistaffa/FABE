#include "main.hpp"

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

// fmt library
#define FMT_HEADER_ONLY
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <fmt/chrono.h>

// print tables after parsing
//#define PRINT_TABLES

// export dot representation of created automata
//#define EXPORT_AUTOMATA_DOT

bool parallel = false;

static inline void print_usage(const char *bin) {

        fmt::print(stderr, "Usage: {} [-h] [-a bub|brz|hop] [-i bound] [-o wmf|mf|miw|md|random|*.pt] ", bin);
        fmt::print(stderr, "[-t unique|random] [-s seed] [-O order] [-r] -f instance\n");
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
        char *order_file = NULL;
        size_t seed = time(NULL);
        bool print_red = false;
        int opt;

        while ((opt = getopt(argc, argv, "a:i:f:o:t:s:O:hr")) != -1) {
                switch (opt) {
                        case 'a':
                                if (strcmp(optarg, "bub") == 0) {
                                        fa_minimization_algorithm = FA_MIN_BUBENZER;
                                } else if (strcmp(optarg, "brz") == 0) {
                                        fa_minimization_algorithm = FA_MIN_BRZOZOWSKI;
                                } else if (strcmp(optarg, "hop") == 0) {
                                        fa_minimization_algorithm = FA_MIN_HOPCROFT;
                                } else {
                                        fmt::print(stderr, "{}: invalid minimisation algorithm -- '{}'\n", argv[0], optarg);
                                        print_usage(argv[0]);
                                        return EXIT_FAILURE;
                                }
                                continue;
                        case 'i':
                                ibound = std::max(0, atoi(optarg));
                                continue;
                        case 'f':
                                if (exists(optarg)) {
                                        instance = optarg;
                                } else {
                                        fmt::print(stderr, "{}: file not found -- '{}'\n", argv[0], optarg);
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
                                        fmt::print(stderr, "{}: file not found -- '{}'\n", argv[0], optarg);
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
                                        fmt::print(stderr, "{}: tie-breaking heuristic not valid -- '{}'\n", argv[0], optarg);
                                        print_usage(argv[0]);
                                        return EXIT_FAILURE;
                                }
                                continue;
                        case 's':
                                seed = atoi(optarg);
                                continue;
                        case 'r':
                                print_red = true;
                                continue;
                        case 'O':
                                order_file = optarg;
                                continue;
                        case 'h':
                        default :
                                print_usage(argv[0]);
                                return EXIT_FAILURE;
                }
        }

        if (!instance) {
                fmt::print(stderr, "{}: instance not specified!\n", argv[0]);
                print_usage(argv[0]);
                return EXIT_FAILURE;
        }

        log_line();
        log_title("Finite-State Automata Bucket Elimination (FABE)");
        log_title("https://github.com/filippobistaffa/FABE");
        log_line();
        log_fmt("Instance", instance);
        int inst_type;

        if (strstr(instance, "wcsp")) {
                log_fmt("Instance type", "WCSP");
                //log_fmt("Bucket elimination algorithm", "MIN-SUM");
                inst_type = WCSP;
        } else {
                log_fmt("Instance type", "UAI");
                //log_fmt("Bucket elimination algorithm", "MIN-SUM with -log()");
                inst_type = MPE;
        }

        std::string algorithms[] = { "Hopcroft", "Brzozowski", "Bubenzer" };
        log_fmt("Automata minimisation algorithm", algorithms[fa_minimization_algorithm]);
        log_fmt("I-bound", (ibound == 0) ? "inf" : fmt::format("{}", ibound));

        if constexpr (QUANTISATION) {
                if (inst_type == WCSP) {
                        log_fmt("Precision", "-");
                } else {
                        log_fmt("Precision", 1 / QUANTISATION);
                }
        } else {
                log_fmt("Precision", "-");
        }

        //log_fmt("Tolerance", tolerance);
        //log_fmt("Parallel mode", parallel ? "Enabled" : "Disabled");

        // look for a known threshold to remove rows
        value threshold = std::numeric_limits<value>::max();

        for (size_t i = 0; i < N_DATASETS; ++i) {
                if (strstr(instance, datasets[i]) != NULL) {
                        threshold = thresholds[i];
                }
        }

        log_fmt("Thresholds value", threshold);

        auto [ domains, adj ] = read_domains_adj(instance, inst_type);
        //print_adj(adj);
        //fmt::print("\n");

        log_fmt("Seed", seed);
        std::string ord_heur_names[] = { "WEIGHTED-MIN-FILL", "MIN-FILL", "MIN-INDUCED-WIDTH", "MIN-DEGREE" };
        std::string tie_heur_names[] = { "MIN-UNIQUENESS", "RANDOM" };
        std::vector<size_t> order;
        auto start_t = std::chrono::high_resolution_clock::now();

        if (pseudotree) {
                log_fmt("Variable order heuristic", pseudotree);
                order = read_pseudotree_order(pseudotree, domains);
        } else {
                srand(seed);
                if (ord_heur == O_RANDOM) {
                        log_fmt("Variable order heuristic", "RANDOM");
                        order.resize(domains.size());
                        iota(order.begin(), order.end(), 0);
                        random_shuffle(order.begin(), order.end());
                } else {
                        log_fmt("Variable order heuristic", ord_heur_names[ord_heur]);
                        log_fmt("Tie-breaking heuristic", tie_heur_names[tie_heur]);
                        order = greedy_order(adj, ord_heur, tie_heur);
                }
        }

        //fmt::print("Order: {}\n", order);
        reverse(order.begin(), order.end());
        std::vector<size_t> pos(order.size());
        for (size_t i = 0; i < order.size(); ++i) {
                pos[order[i]] = i;
        }

        std::chrono::duration<double> runtime = std::chrono::high_resolution_clock::now() - start_t;
        log_fmt("Order computation runtime", fmt::format("{:%T}", runtime));
        start_t = std::chrono::high_resolution_clock::now();
        log_fmt("Induced width", induced_width(adj, order));

        if (order_file) {
                export_order(order, domains, order_file);
                log_fmt("Order file", order_file);
                log_line();
                return EXIT_SUCCESS;
        }

        auto tables = read_tables(instance, inst_type, pos, threshold);

        #ifdef PRINT_TABLES
        fmt::print("\n");
        for (auto const &table : tables) {
                print_table(table);
                fmt::print("\n");
        }
        #endif

        std::vector<automata> automatas(tables.size());
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
                                         1, std::multiplies<size_t>());
        }

        log_string("Value redundancy", fmt::format("{}/{} ({:.3f})", total_rows - actual_rows, total_rows, 1 - actual_rows / total_rows));

        #ifdef EXPORT_AUTOMATA_DOT
        for (auto const &a : automatas) {
                automata_dot(a, "dot");
        }
        #endif

        auto buckets = compute_buckets(automatas, pos);

        log_line();

        if (print_red) {
                exit(0);
        }

        const auto optimal = bucket_elimination(buckets, inst_type == MPE, order, pos, domains, ibound);
        runtime = std::chrono::high_resolution_clock::now() - start_t;
        if (inst_type == WCSP) {
                log_fmt("Solution value", optimal);
        } else {
                log_fmt("Solution value (-log)", fmt::format("{:.3f} ({:.3f})", exp(-(optimal)), optimal));
        }
        log_string("Optimality gap", fmt::format("{} ({:3f}%)", tot_error, 100 * (tot_error) / optimal));
        log_fmt("Bucket elimination runtime", fmt::format("{:%T}", runtime));
        log_line();

        return EXIT_SUCCESS;
}
