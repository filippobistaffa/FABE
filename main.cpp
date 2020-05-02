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

// print tables after parsing
//#define PRINT_TABLES

// export dot representation of created automata
//#define EXPORT_AUTOMATA_DOT

bool parallel = false;

static inline void print_usage(const char *bin) {

        cerr << "Usage: " << bin << " [-h] [-a bub|brz|hop] [-i bound] [-e tolerance] -f instance" << endl;
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
        char *instance = NULL;
        int opt;

        while ((opt = getopt(argc, argv, "a:i:f:e:h")) != -1) {
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
                        case 'e':
                                tolerance = atof(optarg);
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

        string algorithms[] = { "Hopcroft", "Brzozowski", "Bubenzer" };
        log_value("Automata minimisation algorithm", algorithms[fa_minimization_algorithm]);

        log_value("I-bound", (ibound == 0) ? "inf" : to_string(ibound));
        log_value("Tolerance", tolerance);
        //log_value("Parallel mode", parallel ? "Enabled" : "Disabled");

        auto [ domains, adj ] = read_domains_adj(instance);
        //print_adj(adj);
        //cout << endl;

        vector<size_t> order(domains.size());
        iota(order.begin(), order.end(), 0);
        reverse(order.begin(), order.end());
        vector<size_t> pos(order.size());
        for (size_t i = 0; i < order.size(); ++i) {
                pos[order[i]] = i;
        }

        auto start_t = chrono::high_resolution_clock::now();
        //export_order(order, domains, "order.vo");
        log_value("Induced width", induced_width(adj, order));
        auto tables = read_tables(instance);

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

        log_value("Value redundancy", 1 - actual_rows / total_rows);

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
        //const int inner = (inst_type == WCSP) ? BE_SUM : BE_PROD;
        //const int outer = (inst_type == WCSP) ? BE_MIN : BE_MAX;
        const auto optimal = bucket_elimination(buckets, BE_SUM, BE_MIN, order, pos, domains, ibound);
        chrono::duration<double> runtime = chrono::high_resolution_clock::now() - start_t;
        log_value("Solution value", optimal);
        //log_value("Maximum optimality gap", tolerance * tables.size());
        //log_value("Maximum optimality gap (%)", 100 * (tolerance * tables.size()) / optimal);
        //log_value("Optimality gap", tot_error);
        ostringstream oss;
        oss << tot_error << " (" << setprecision(3) << ((tot_error) ? 100 * (tot_error) / optimal : tot_error) << "%)";
        log_value("Optimality gap", oss.str());
        log_value("Bucket elimination runtime", runtime.count());
        log_line();

        return EXIT_SUCCESS;
}
