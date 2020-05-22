#include "be.hpp"

#include <numeric>      // accumulate
#include <algorithm>    // max_element, sort

#include "libfa/fa.h"
#include "order.hpp"
#include "log.hpp"
#include "util.hpp"
#include "export.hpp"

// print tables during execution
//#define PRINT_TABLES

// print bucket debug messages
//#define DEBUG_BUCKETS

// enable profiling
//#define CPU_PROFILER
//#define HEAP_PROFILER
#define COUNT_STATES

#ifdef CPU_PROFILER
#define CPU_PROFILER_OUTPUT "trace.prof"
#include <gperftools/profiler.h>
#endif

#ifdef HEAP_PROFILER
#define HEAP_PROFILER_PREFIX "memory/memory"
#include <gperftools/heap-profiler.h>
#endif

#ifdef COUNT_STATES
size_t tot_states;
#endif

#ifdef PRINT_TABLES
#include "conversion.hpp"
#include "io.hpp"
#endif

#define OP_FREE_OLD(OP, FREE, X, Y) { auto __TMP = (X); (X) = OP(X, Y); FREE(__TMP); }
#define SET_OP(OP, X, Y, R, CMP) (OP((X).begin(), (X).end(), (Y).begin(), (Y).end(), inserter((R), (R).begin()), CMP))

extern bool parallel;
bool debug = false;

static inline automata join(automata &a1, automata &a2, int inner, vector<size_t> const &domains) {

        #ifdef PRINT_TABLES
        cout << endl << "Joining:" << endl << endl;
        print_table(compute_table(a1));
        cout << endl;
        print_table(compute_table(a2));
        cout << endl;
        #endif

        automata join;

        // compute variables (and their domains) of join function
        SET_OP(set_union, a1.vars, a2.vars, join.vars, greater<size_t>());

        for (auto var : join.vars) {
                join.domains.push_back(domains[var]);
        }

        // variables to be added
        vector<size_t> add1;
        vector<size_t> add2;
        SET_OP(set_difference, a2.vars, a1.vars, add1, greater<size_t>());
        SET_OP(set_difference, a1.vars, a2.vars, add2, greater<size_t>());

        // non-shared variables' positions in union
        vector<size_t> padd1;
        vector<size_t> padd2;

        for (auto var : add1) {
                padd1.push_back(lower_bound(join.vars.begin(), join.vars.end(), var, greater<size_t>()) - join.vars.begin());
        }

        for (auto var : add2) {
                padd2.push_back(lower_bound(join.vars.begin(), join.vars.end(), var, greater<size_t>()) - join.vars.begin());
        }

        /*
        cout << vec2str(a1.vars, "v1") << endl;
        cout << vec2str(a2.vars, "v2") << endl;
        cout << vec2str(join.vars, "vj") << endl;
        cout << vec2str(add1, "add1") << endl;
        cout << vec2str(add2, "add2") << endl;
        cout << vec2str(padd1, "padd1") << endl;
        cout << vec2str(padd2, "padd2") << endl;
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

        vector<value> keys;

        for (auto &[ v1, fa1 ] : a1.rows) {
                for (auto &[ v2, fa2 ] : a2.rows) {
                        auto in = fa_intersect(fa1, fa2);
                        if (fa_is_basic(in, FA_EMPTY)) { // these two rows do not have any common variable
                                fa_free(in);             // assignment of shared variables
                        } else {
                                const value v1v2 = (inner == BE_SUM) ? v1 + v2 : v1 * v2;
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

        //#pragma omp parallel for schedule(dynamic) if (parallel)
        for (size_t i = 0; i < keys.size(); ++i) {
                fa_minimize(join.rows[keys[i]]);
        }

        #ifdef PRINT_TABLES
        cout << "Result:" << endl << endl;
        print_table(compute_table(join));
        cout << endl;
        #endif
        return join;
}

static inline void free_automata(automata &a) {

        for (auto& [ v, fa ] : a.rows) {
                fa_free(fa);
        }
        a.rows.clear();
}

static inline automata join_bucket(vector<automata> &bucket, int inner, vector<size_t> const &pos, vector<size_t> const &domains) {

        auto res = bucket.front();

	for (auto it = next(bucket.begin()); it != bucket.end(); ++it) {
                auto old = res;
                #ifdef HEAP_PROFILER
                HeapProfilerDump("pre-join");
                #endif
                res = join(old, *it, inner, domains);
                #ifdef HEAP_PROFILER
                HeapProfilerDump("post-join");
                #endif
                free_automata(old);
                free_automata(*it);
                #ifdef HEAP_PROFILER
                HeapProfilerDump("post-free");
                #endif
                #ifdef DEBUG_BUCKETS
                cout << "Joined " << (it - bucket.begin()) + 1 << " functions" << endl;
                #endif
	}

        #ifdef COUNT_STATES
        for (auto &[ v, fa ] : res.rows) {
                tot_states += fa_n_states(fa);
        }
        #endif

        return res;
}

static inline value reduce_var(automata &a, size_t var, int outer) {

        #ifdef PRINT_TABLES
        cout << endl << "Minimising:" << endl << endl;
        print_table(compute_table(a));
        cout << endl;
        #endif

        const size_t idx = lower_bound(a.vars.begin(), a.vars.end(), var, greater<size_t>()) - a.vars.begin();

        if (debug) {
                //automata_dot(a, "dot");
        }

        a.vars.erase(a.vars.begin() + idx);
        a.domains.erase(a.domains.begin() + idx);
        vector<value> keys;

        for (auto &[ v, fa ] : a.rows) {
                if (idx == a.vars.size()) {
                        if (fa_remove_last(fa) > 1) {
                                fa_merge_accept(fa);
                        }
                } else {
                        fa_remove_level(fa, idx);
                }
                keys.push_back(v);
        }

        if (debug) {
                //automata_dot(a, "dot");
        }

        if (outer == BE_MIN) {
                sort(keys.begin(), keys.end());
        } else {
                sort(keys.begin(), keys.end(), greater<value>());
        }

        if (a.vars.size() == 0) {
                return keys.front();
        }

        vector<struct fa *> pfx_union(keys.size());
        pfx_union[0] = fa_make_basic(FA_EMPTY);

        for (size_t i = 1; i < keys.size(); ++i) {
                pfx_union[i] = fa_union(a.rows[keys[i - 1]], pfx_union[i - 1]);
                fa_minimize(pfx_union[i]);
        }

        vector<value> empty;

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
        cout << endl;
        #endif
        return 0;
}

static inline size_t push_bucket(automata const &a, vector<vector<automata>> &buckets, vector<size_t> const &pos) {

        const size_t b = *max_element(a.vars.begin(), a.vars.end(), compare_vec(pos));
        buckets[b].push_back(move(a));
        return b;
}

static inline value process_bucket(size_t var, vector<automata> &bucket, vector<vector<automata>> &buckets, int inner,
                                   int outer, vector<size_t> const &pos, vector<size_t> const &domains,
                                   struct bin_data &mbe) {

        value res = 0;

        if (bucket.size()) {
                auto h = join_bucket(bucket, inner, pos, domains);
                if (h.rows.size() > 0) {
                        const value v_h = reduce_var(h, var, outer);
                        res += v_h;
                        //automata_dot(h, "dot");
                        if (h.vars.size() > 0) {
                                size_t b = push_bucket(h, buckets, pos);
                                #ifdef DEBUG_BUCKETS
                                cout << "Result placed in bucket " << b << endl;
                                #endif
                                if (mbe.filename) {
                                        export_function(h, var, b, mbe);
                                }
                        } else {
                                if (mbe.filename) {
                                        export_root_function(v_h, var, mbe);
                                }
                        }
                }
        }

        return res;
}

#define COMPARE_MBE

static inline vector<vector<automata>> mini_buckets(vector<automata> &bucket, size_t ibound,
                                                    vector<size_t> const &pos) {

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
        vector<vector<size_t>> bucket_vars = { bucket.front().vars };
        vector<vector<automata>> mini_buckets = { { move(bucket.front()) } };

        for (auto it = next(bucket.begin()); it != bucket.end(); ++it) {
                size_t mb = 0;
                for (; mb < mini_buckets.size(); ++mb) {
                        vector<size_t> tmp;
                        SET_OP(set_union, it->vars, bucket_vars[mb], tmp, greater<size_t>());
                        if (tmp.size() == bucket_vars[mb].size() || tmp.size() <= ibound + 1) {
                                // function fits in this mini-bucket
                                mini_buckets[mb].push_back(move(*it));
                                bucket_vars[mb] = tmp;
                                break;
                        }
                }
                if (mb == mini_buckets.size()) { // could not find any bucket
                        bucket_vars.push_back({ it->vars });
                        mini_buckets.push_back({ move(*it) });
                }
        }

        return mini_buckets;
}

vector<vector<automata>> compute_buckets(vector<automata> const &automatas, vector<size_t> const &pos) {

        vector<vector<automata>> buckets(pos.size(), vector<automata>());

        for (automata const &a : automatas) {
                push_bucket(a, buckets, pos);
        }

        return buckets;
}

value bucket_elimination(vector<vector<automata>> &buckets, int inner, int outer,
                         vector<size_t> const &order, vector<size_t> const &pos,
                         vector<size_t> const &domains, size_t ibound, struct bin_data &mbe) {

        #ifdef CPU_PROFILER
        ProfilerStart(CPU_PROFILER_OUTPUT);
        #endif

        #ifdef HEAP_PROFILER
        HeapProfilerStart(HEAP_PROFILER_PREFIX);
        #endif

        value optimal = 0;

        size_t i = 0;

        for (auto it = order.rbegin(); it != order.rend(); ++it) {
                if (i == 11) {
                        debug = true;
                }
                #ifdef DEBUG_BUCKETS
                cout << "Processing bucket " << *it << " with " << buckets[*it].size() << " functions" << endl;
                #endif
                if (ibound && buckets[*it].size() > 1) {
                        auto mb = mini_buckets(buckets[*it], ibound, pos);
                        #ifdef DEBUG_BUCKETS
                        cout << "There are " << mb.size() << " mini-buckets" << endl;
                        #endif
                        for (auto &bucket : mb) {
                                optimal += process_bucket(*it, bucket, buckets, inner, outer, pos, domains, mbe);
                        }
                } else {
                        optimal += process_bucket(*it, buckets[*it], buckets, inner, outer, pos, domains, mbe);
                }
                #ifndef DEBUG_BUCKETS
                log_progress_increase(1, order.size());
                #endif
                
                if (i == 11) {
                        //break;
                }
                i++;
        }

        #ifdef CPU_PROFILER
        ProfilerStop();
        #endif

        #ifdef HEAP_PROFILER
        HeapProfilerStop();
        #endif

        log_line();

        #ifdef COUNT_STATES
        log_value("Total number of automata states", tot_states);
        #endif

        // export to binary file
        if (mbe.filename) {
                write_binary(mbe, optimal, ibound, order.front());
        }

        return optimal;
}
