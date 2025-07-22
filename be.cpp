#include "be.hpp"

#include <numeric>      // accumulate
#include <algorithm>    // max_element, sort

#include "libfa/fa.h"
#include "order.hpp"
#include "log.hpp"

// print tables during execution
//#define PRINT_TABLES

// print bucket debug messages
//#define DEBUG_BUCKETS

// enable profiling
//#define CPU_PROFILER
//#define HEAP_PROFILER

#ifdef CPU_PROFILER
#define CPU_PROFILER_OUTPUT "trace.prof"
#include <gperftools/profiler.h>
#endif

#ifdef HEAP_PROFILER
#define HEAP_PROFILER_PREFIX "memory/memory"
#include <gperftools/heap-profiler.h>
#endif

#ifdef PRINT_TABLES
#include "conversion.hpp"
#include "io.hpp"
#endif

#define OP_FREE_OLD(OP, FREE, X, Y) { auto __TMP = (X); (X) = OP(X, Y); FREE(__TMP); }
#define SET_OP(OP, X, Y, R, CMP) (OP((X).begin(), (X).end(), (Y).begin(), (Y).end(), inserter((R), (R).begin()), CMP))

extern bool parallel;
size_t tot_states;
size_t tot_keys;

value quantise(value val) {

    if constexpr (QUANTISATION) {
        const value q = floor(exp(-val) * QUANTISATION);
        return -log(q / QUANTISATION);
    } else {
        return val;
    }
}

static inline automata join(automata &a1, automata &a2, bool quant, std::vector<size_t> const &pos,
                            std::vector<size_t> const &domains) {

    #ifdef PRINT_TABLES
    fmt::print("\nJoining:\n\n");
    print_table(compute_table(a1));
    fmt::print("\n");
    print_table(compute_table(a2));
    fmt::print("\n");
    #endif

    automata join;

    // compute variables (and their domains) of join function
    SET_OP(set_union, a1.vars, a2.vars, join.vars, compare_vec(pos));

    for (auto var : join.vars) {
        join.domains.push_back(domains[var]);
    }

    // variables to be added
    std::vector<size_t> add1;
    std::vector<size_t> add2;
    SET_OP(set_difference, a2.vars, a1.vars, add1, compare_vec(pos));
    SET_OP(set_difference, a1.vars, a2.vars, add2, compare_vec(pos));

    // non-shared variables' positions in union
    std::vector<size_t> padd1;
    std::vector<size_t> padd2;

    for (auto var : add1) {
        padd1.push_back(lower_bound(join.vars.begin(), join.vars.end(), var, compare_vec(pos)) - join.vars.begin());
    }

    for (auto var : add2) {
        padd2.push_back(lower_bound(join.vars.begin(), join.vars.end(), var, compare_vec(pos)) - join.vars.begin());
    }

    /*
    fmt::print("{}\n", a1.vars, "v1");
    fmt::print("{}\n", a2.vars, "v2");
    fmt::print("{}\n", join.vars, "vj");
    fmt::print("{}\n", add1, "add1");
    fmt::print("{}\n", add2, "add2");
    fmt::print("{}\n", padd1, "padd1");
    fmt::print("{}\n", padd2, "padd2");
    */

    for (auto &[ v, fa ] : a1.rows) {
        for (size_t i = 0; i < add1.size(); ++i) {
            fa_add_level(fa, padd1[i], domains[add1[i]]);
        }
    }


    for (auto &[ v, fa ] : a2.rows) {
        for (size_t i = 0; i < add2.size(); ++i) {
            fa_add_level(fa, padd2[i], domains[add2[i]]);
        }
    }

    std::vector<value> keys;

    for (auto &[ v1, fa1 ] : a1.rows) {
        for (auto &[ v2, fa2 ] : a2.rows) {
            const value v1v2 = (quant) ? quantise(v1 + v2) : v1 + v2;
            if (std::isinf(v1v2)) {
                continue;
            }
            auto in = fa_intersect(fa1, fa2);
            if (fa_is_basic(in, FA_EMPTY)) { // these two rows do not have any common variable
                fa_free(in);         // assignment of shared variables
            } else {
                auto it = join.rows.find(v1v2);
                if (it == join.rows.end()) {
                    join.rows.insert({ v1v2, in });
                    keys.push_back(v1v2);
                } else {
                    fa_union_in_place(it->second, &in);
                }
            }
        }
    }

    //sort(keys.begin(), keys.end());
    //fmt::print("{}\n", keys, "keys");
    tot_keys += keys.size();

    //#pragma omp parallel for schedule(dynamic) if (parallel)
    for (size_t i = 0; i < keys.size(); ++i) {
        fa_minimize(join.rows[keys[i]]);
    }

    #ifdef PRINT_TABLES
    fmt::print("Result:\n\n");
    print_table(compute_table(join));
    fmt::print("\n");
    #endif
    return join;
}

static inline void free_automata(automata &a) {

    for (auto& [ v, fa ] : a.rows) {
        fa_free(fa);
    }
    a.rows.clear();
}

static inline automata join_bucket(std::vector<automata> &bucket, bool quant, std::vector<size_t> const &pos, std::vector<size_t> const &domains) {

    auto res = bucket.front();

	for (auto it = next(bucket.begin()); it != bucket.end(); ++it) {
        auto old = res;
        #ifdef HEAP_PROFILER
        HeapProfilerDump("pre-join");
        #endif
        res = join(old, *it, quant, pos, domains);
        #ifdef HEAP_PROFILER
        HeapProfilerDump("post-join");
        #endif
        free_automata(old);
        free_automata(*it);
        #ifdef HEAP_PROFILER
        HeapProfilerDump("post-free");
        #endif
        #ifdef DEBUG_BUCKETS
        fmt::print("Joined {} functions\n", (it - bucket.begin()) + 1);
        #endif
	}

    for (auto &[ v, fa ] : res.rows) {
        tot_states += fa_n_states(fa);
    }

    return res;
}

static inline value reduce_last_var(automata &a) {

    #ifdef PRINT_TABLES
    fmt::print("\nMinimising:\n\n");
    print_table(compute_table(a));
    fmt::print("\n");
    #endif

    // remove last variable and domain
    a.vars.pop_back();
    a.domains.pop_back();
    std::vector<value> keys;

    for (auto &[ v, fa ] : a.rows) {
        if (fa_remove_last(fa) > 1) {
            fa_merge_accept(fa);
        }
        keys.push_back(v);
    }

    sort(keys.begin(), keys.end());

    if (a.vars.size() == 0) {
        return keys.front();
    }

    std::vector<struct fa *> pfx_union(keys.size());
    pfx_union[0] = fa_make_basic(FA_EMPTY);

    for (size_t i = 1; i < keys.size(); ++i) {
        pfx_union[i] = fa_union(a.rows[keys[i - 1]], pfx_union[i - 1]);
        fa_minimize(pfx_union[i]);
    }

    std::vector<value> empty;

    for (size_t i = 1; i < keys.size(); ++i) {
        OP_FREE_OLD(fa_minus, fa_free, a.rows[keys[i]], pfx_union[i]);
        if (fa_is_basic(a.rows[keys[i]], FA_EMPTY)) {
            empty.push_back(keys[i]);
        }
    }

    for (auto e : empty) {
        fa_free(a.rows[e]);
        a.rows.erase(e);
    }

    for (auto a : pfx_union) {
        fa_free(a);
    }

    #ifdef PRINT_TABLES
    print_table(compute_table(a));
    fmt::print("\n");
    #endif
    return 0;
}

static inline size_t push_bucket(automata const &a, std::vector<std::vector<automata>> &buckets, std::vector<size_t> const &pos) {

    const size_t b = *max_element(a.vars.begin(), a.vars.end(), compare_vec(pos));
    buckets[b].push_back(std::move(a));
    return b;
}

static inline value process_bucket(std::vector<automata> &bucket, std::vector<std::vector<automata>> &buckets, bool quant,
                                   std::vector<size_t> const &pos, std::vector<size_t> const &domains) {

    if (bucket.size()) {
        auto h = join_bucket(bucket, quant, pos, domains);
        if (h.rows.size() > 0) {
            //automata_dot(h, "dot");
            const value res = reduce_last_var(h);
            if (h.vars.size() > 0) {
                //automata_dot(h, "dot");
                #ifdef DEBUG_BUCKETS
                fmt::print("Result placed in bucket {}\n", push_bucket(h, buckets, pos));
                #else
                push_bucket(h, buckets, pos);
                #endif
            }
            return res;
        } else { // no assignment is in all functions, some constraint is violated
            return std::numeric_limits<value>::infinity();
        }
    } else {
        return 0;
    }
}

//#define COMPARE_MBE

static inline std::vector<std::vector<automata>> mini_buckets(std::vector<automata> &bucket, size_t ibound, std::vector<size_t> const &pos) {

    // doing a simple FFD bin packing
    sort(bucket.begin(), bucket.end(), [](automata const &x, automata const &y) {
        #ifdef COMPARE_MBE
        if (x.vars.size() == y.vars.size()) {
            auto vx = x.vars;
            auto vy = y.vars;
            sort(vx.begin(), vx.end());
            sort(vy.begin(), vy.end());
            return lexicographical_compare(vx.begin(), vx.end(), vy.begin(), vy.end());
        } else
        #endif
        return (x.vars.size() > y.vars.size());
    });

    // initialise first mini-bucket with first (larger) function
    std::vector<std::vector<size_t>> bucket_vars = { bucket.front().vars };
    std::vector<std::vector<automata>> mini_buckets = { { std::move(bucket.front()) } };

    for (auto it = next(bucket.begin()); it != bucket.end(); ++it) {
        size_t mb = 0;
        for (; mb < mini_buckets.size(); ++mb) {
            std::vector<size_t> tmp;
            SET_OP(set_union, it->vars, bucket_vars[mb], tmp, compare_vec(pos));
            if (tmp.size() == bucket_vars[mb].size() || tmp.size() <= ibound + 1) {
                // function fits in this mini-bucket
                mini_buckets[mb].push_back(std::move(*it));
                bucket_vars[mb] = tmp;
                break;
            }
        }
        if (mb == mini_buckets.size()) { // could not find any bucket
            bucket_vars.push_back({ it->vars });
            mini_buckets.push_back({ std::move(*it) });
        }
    }

    return mini_buckets;
}

std::vector<std::vector<automata>> compute_buckets(std::vector<automata> const &automatas, std::vector<size_t> const &pos) {

    std::vector<std::vector<automata>> buckets(pos.size(), std::vector<automata>());

    for (automata const &a : automatas) {
        push_bucket(a, buckets, pos);
    }

    return buckets;
}

value bucket_elimination(std::vector<std::vector<automata>> &buckets, bool quant,
                         std::vector<size_t> const &order, std::vector<size_t> const &pos,
                         std::vector<size_t> const &domains, size_t ibound) {

    #ifdef CPU_PROFILER
    ProfilerStart(CPU_PROFILER_OUTPUT);
    #endif

    #ifdef HEAP_PROFILER
    HeapProfilerStart(HEAP_PROFILER_PREFIX);
    #endif

    value optimal = 0;

    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        #ifdef DEBUG_BUCKETS
        fmt::print("Processing bucket {} with {} functions\n", *it, buckets[*it].size());
        #endif
        if (ibound && buckets[*it].size() > 1) {
            auto mb = mini_buckets(buckets[*it], ibound, pos);
            #ifdef DEBUG_BUCKETS
            fmt::print("There are {} mini-buckets\n", mb.size());
            #endif
            for (auto &bucket : mb) {
                optimal += process_bucket(bucket, buckets, quant, pos, domains);
            }
        } else {
            optimal += process_bucket(buckets[*it], buckets, quant, pos, domains);
        }
        #ifndef DEBUG_BUCKETS
        log_progress_increase(1, order.size());
        #endif
    }

    #ifdef CPU_PROFILER
    ProfilerStop();
    #endif

    #ifdef HEAP_PROFILER
    HeapProfilerStop();
    #endif

    log_line();
    log_fmt("Total number of automata states", tot_states);
    log_fmt("Total number of keys", tot_keys);

    return optimal;
}
