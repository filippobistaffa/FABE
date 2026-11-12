#include <limits>       // std::numeric_limits
#include <filesystem>   // std::filesystem::exists
#include <chrono>       // time measurement
#include <numeric>      // accumulate
#include <string.h>     // strcmp
#include <unistd.h>     // getopt
#include <algorithm>    // reverse
#include <random>       // shuffle

#include "types.hpp"
#include "be.hpp"
#include "io.hpp"
#include "log.hpp"
#include "order.hpp"
#include "conversion.hpp"

// fmt library
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <fmt/color.h>
#include <fmt/chrono.h>

#include <cstddef>
#include <functional>
#include <algorithm>
#include <utility>
#include <vector>
#include <random>
#include <map>

inline table random_table(unsigned seed, size_t max_vars, size_t max_domain, size_t max_value) {

    std::mt19937_64 rng(seed);
    auto uniform = [&](size_t lo, size_t hi) {
        return std::uniform_int_distribution<size_t>(lo, hi)(rng);
    };

    table t;
    const size_t n = uniform(2, max_vars);
    for (size_t i = 0; i < n; ++i) {
        t.vars.push_back(i);
        t.domains.push_back(uniform(2, max_domain));
    }

    std::vector<size_t> a(n, 0);
    while (true) {
        t.rows.emplace_back(a, static_cast<value>(uniform(0, max_value)));
        size_t i = n;
        while (i > 0 && ++a[i - 1] == t.domains[i - 1]) a[--i] = 0;
        if (i == 0) break;
    }
    return t;
}

inline table example_table() {
    return table{
        /* vars    */ {0, 1, 2},
        /* domains */ {3, 3, 3},
        /* rows    */ {
            {{0, 0, 0}, 3}, {{0, 0, 1}, 3}, {{0, 0, 2}, 3},
            {{0, 1, 0}, 3}, {{0, 1, 1}, 3}, {{0, 1, 2}, 3},
            {{0, 2, 0}, 3}, {{0, 2, 1}, 3}, {{0, 2, 2}, 3},
            {{1, 0, 0}, 3}, {{1, 0, 1}, 3}, {{1, 0, 2}, 3},
            {{1, 1, 0}, 3}, {{1, 1, 1}, 3}, {{1, 1, 2}, 7},
            {{1, 2, 0}, 7}, {{1, 2, 1}, 7}, {{1, 2, 2}, 7},
            {{2, 0, 0}, 0}, {{2, 0, 1}, 1}, {{2, 0, 2}, 5},
            {{2, 1, 0}, 4}, {{2, 1, 1}, 4}, {{2, 1, 2}, 4},
            {{2, 2, 0}, 4}, {{2, 2, 1}, 2}, {{2, 2, 2}, 2},
        }
    };
}

template <typename Op = std::plus<value>>
table reduce_last(const table& t, Op op = {}) {

    table out;
    out.vars.assign(t.vars.begin(), t.vars.end() - 1);
    out.domains.assign(t.domains.begin(), t.domains.end() - 1);
    std::map<std::vector<size_t>, value> acc;

    for (const auto& [assignment, v] : t.rows) {
        std::vector<size_t> key(assignment.begin(), assignment.end() - 1);
        auto [it, inserted] = acc.try_emplace(std::move(key), v);
        if (!inserted) it->second = op(it->second, v);
    }

    out.rows.reserve(acc.size());
    for (auto& [key, v] : acc) out.rows.emplace_back(key, v);
    return out;
}

// print tables during execution
//#define PRINT_TABLES

// compare against reference table implementation
//#define CHECK_CORRECT

int main(int argc, char** argv) {

    if (argc < 2) {
        return 1;
    }

    const value tolerance = 0;
    const unsigned seed = static_cast<unsigned>(std::stoul(argv[1]));
    const unsigned n = static_cast<unsigned>(std::stoul(argv[2]));

    auto table = random_table(seed, n, n, 10);
    //auto table = example_table();

    #ifdef PRINT_TABLES
    print_table(table);
    fmt::print("\n");
    #endif

    sort(table.rows.begin(), table.rows.end(), [](
        std::pair<std::vector<size_t>, value> const &x, std::pair<std::vector<size_t>, value> const &y
    ) { return (x.second < y.second); });

    std::vector<size_t> order = table.vars;
    std::vector<size_t> domains = table.domains;
    std::vector<size_t> pos(order.size());
    for (size_t i = 0; i < order.size(); ++i) {
        pos[order[i]] = i;
    }

    #ifdef CHECK_CORRECT
    auto reduced_table = reduce_last(table);
    sort(reduced_table.rows.begin(), reduced_table.rows.end());
    #ifdef PRINT_TABLES
    print_table(reduced_table);
    fmt::print("\n");
    #endif
    fmt::print("SEED = {}, ", seed);
    #endif
    auto start_t = std::chrono::high_resolution_clock::now();
    auto [ a, error ] = compute_automata(table, tolerance);
    auto reduce_t = std::chrono::high_resolution_clock::now();
    auto reduced_automata = map_reduce_last_var(a, false, pos, domains);
    auto reduce_runtime = std::chrono::high_resolution_clock::now() - reduce_t;
    auto constr_runtime = std::chrono::high_resolution_clock::now() - start_t - reduce_runtime;
    #ifdef CHECK_CORRECT
    auto reduced_automata_table = compute_table(reduced_automata);
    sort(reduced_automata_table.rows.begin(), reduced_automata_table.rows.end());
    #ifdef PRINT_TABLES
    print_table(reduced_automata_table);
    fmt::print("\n");
    #endif
    if (reduced_table.rows == reduced_automata_table.rows) {
        fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "SUCCESS\n");
    } else {
        fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, "FAILED\n");
        return EXIT_FAILURE;
    }
    #else
    fmt::print("Seed = {}\n", seed);
    fmt::print("Construction = {:%T}\n", constr_runtime);
    fmt::print("Filter = {:%T}\n", reduce_runtime);
    #endif
    return EXIT_SUCCESS;
}
