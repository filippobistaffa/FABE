#include "be.hpp"

#define REDUCTION_MIN
//#define REDUCTION_MAX

// print tables during execution
//#define PRINT_TABLES

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

#define JOIN_OP(X, Y) (X + Y)
#define OP_FREE_OLD(OP, FREE, X, Y) { auto __TMP = (X); (X) = OP(X, Y); FREE(__TMP); }
#define SET_OP(OP, X, Y, R, CMP) (OP((X).begin(), (X).end(), (Y).begin(), (Y).end(), inserter((R), (R).begin()), CMP))

extern bool parallel;

static inline automata join(automata &a1, automata &a2, vector<size_t> const &pos, vector<size_t> const &domains) {

        #ifdef PRINT_TABLES
        cout << endl << "Joining:" << endl << endl;
        print_table(compute_table(a1));
        cout << endl;
        print_table(compute_table(a2));
        cout << endl;
        #endif

        automata join;

        // compute variables (and their domains) of join function
        SET_OP(set_union, a1.vars, a2.vars, join.vars, compare_pos(pos));

        for (auto var : join.vars) {
                join.domains.push_back(domains[var]);
        }

        // variables to be added
        vector<size_t> add1;
        vector<size_t> add2;
        SET_OP(set_difference, a2.vars, a1.vars, add1, compare_pos(pos));
        SET_OP(set_difference, a1.vars, a2.vars, add2, compare_pos(pos));

        // non-shared variables' positions in union
        vector<size_t> padd1;
        vector<size_t> padd2;

        for (auto var : add1) {
                padd1.push_back(lower_bound(join.vars.begin(), join.vars.end(), var, compare_pos(pos)) - join.vars.begin());
        }

        for (auto var : add2) {
                padd2.push_back(lower_bound(join.vars.begin(), join.vars.end(), var, compare_pos(pos)) - join.vars.begin());
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
                for (auto i = 0; i < add1.size(); ++i) {
                        fa_add_level(fa, padd1[i], ALPHABET[domains[add1[i]] - 1]);
                }
        }

        for (auto &[ v, fa ] : a2.rows) {
                for (auto i = 0; i < add2.size(); ++i) {
                        fa_add_level(fa, padd2[i], ALPHABET[domains[add2[i]] - 1]);
                }
        }

        for (auto &[ v1, fa1 ] : a1.rows) {
                for (auto &[ v2, fa2 ] : a2.rows) {
                        auto in = fa_intersect(fa1, fa2);
                        if (fa_is_basic(in, FA_EMPTY)) { // these two rows do not have any common variable
                                fa_free(in);             // assignment of shared variables
                        } else {
                                const value v1v2 = JOIN_OP(v1, v2);
                                auto it = join.rows.find(v1v2);
                                if (it == join.rows.end()) {
                                        join.rows.insert({ v1v2, in });
                                } else {
                                        fa_union_in_place(it->second, &in);
                                }
                        }
                }
        }

        for (auto &[ v, fa ] : join.rows) {
                fa_minimize(fa);
        }

        #ifdef PRINT_TABLES
        cout << "Result:" << endl << endl;
        print_table(compute_table(join));
        #endif
        return join;
}

static inline automata join_bucket(vector<automata> &bucket, vector<size_t> const &pos, vector<size_t> const &domains) {

        auto res = bucket.front();

	for (auto it = next(bucket.begin()); it != bucket.end(); ++it) {
                auto old = res;
                #ifdef HEAP_PROFILER
                HeapProfilerDump("pre-join");
                #endif
                res = join(old, *it, pos, domains);
                #ifdef HEAP_PROFILER
                HeapProfilerDump("post-join");
                #endif
                free_automata(old);
                free_automata(*it);
                #ifdef HEAP_PROFILER
                HeapProfilerDump("post-free");
                #endif
                cout << "Joined " << (it - bucket.begin()) + 1 << " functions" << endl;
	}

        return res;
}

static inline value reduce_last_var(automata &a) {

        #ifdef PRINT_TABLES
        cout << endl << "Minimising:" << endl << endl;
        print_table(compute_table(a));
        cout << endl;
        #endif

        // remove last variable and domain
        a.vars.pop_back();
        a.domains.pop_back();
        vector<value> keys;

        for (auto &[ v, fa ] : a.rows) {
                if (fa_remove_last(fa) > 1) {
                        fa_merge_accept(fa);
                }
                keys.push_back(v);
        }

        #ifdef REDUCTION_MIN
        sort(keys.begin(), keys.end());
        #else
        sort(keys.begin(), keys.end(), greater<value>());
        #endif

        if (a.vars.size() == 0) {
                return keys.front();
        }

        //START_CLOCK;

        vector<struct fa *> pfx_union(keys.size());
        pfx_union[0] = fa_make_basic(FA_EMPTY);

        for (auto i = 1; i < keys.size(); ++i) {
                pfx_union[i] = fa_union(a.rows[keys[i - 1]], pfx_union[i - 1]);
                fa_minimize(pfx_union[i]);
        }

        //STOP_CLOCK;

        vector<value> empty;

        for (auto i = 1; i < keys.size(); ++i) {
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

static inline value process_bucket(vector<automata> &bucket, vector<vector<automata>> &buckets,
                                   vector<size_t> const &pos, vector<size_t> const &domains) {

        value res = 0;

        if (bucket.size()) {
                auto h = join_bucket(bucket, pos, domains);
                //automata_dot(h, "dot");
                res += reduce_last_var(h);
                if (h.vars.size() > 0) {
                        //automata_dot(h, "dot");
                        cout << "Result placed in bucket " << push_bucket(h, buckets, pos) << endl << endl;
                }
        }

        return res;
}

static inline vector<vector<automata>> mini_buckets(vector<automata> &bucket, size_t ibound,
                                                    vector<size_t> const &pos) {

        // doing a simple FFD bin packing
        sort(bucket.begin(), bucket.end(), [](automata const &x, automata const &y) {
                return (x.vars.size() > y.vars.size());
        });

        // initialise first mini-bucket with first (larger) function
        vector<vector<size_t>> bucket_vars = { bucket.front().vars };
        vector<vector<automata>> mini_buckets = { { move(bucket.front()) } };

        for (auto it = next(bucket.begin()); it != bucket.end(); ++it) {
                auto mb = 0;
                for (; mb < mini_buckets.size(); ++mb) {
                        vector<size_t> tmp;
                        SET_OP(set_union, it->vars, bucket_vars[mb], tmp, compare_pos(pos));
                        if (tmp.size() <= ibound) { // function fits in this mini-bucket
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

value bucket_elimination(vector<vector<automata>> &buckets, vector<size_t> const &order,
                         vector<size_t> const &pos, vector<size_t> const &domains, size_t ibound) {

        #ifdef CPU_PROFILER
        ProfilerStart(CPU_PROFILER_OUTPUT);
        #endif

        #ifdef HEAP_PROFILER
        HeapProfilerStart(HEAP_PROFILER_PREFIX);
        #endif

        value optimal = 0;

        for (auto it = order.rbegin(); it != order.rend(); ++it) {
                cout << "Processing bucket " << *it << " with " << buckets[*it].size() << " functions" << endl;
                if (ibound && buckets[*it].size() > 1) {
                        auto mb = mini_buckets(buckets[*it], ibound, pos);
                        for (auto &bucket : mb) {
                                optimal += process_bucket(bucket, buckets, pos, domains);
                        }
                } else {
                        optimal += process_bucket(buckets[*it], buckets, pos, domains);
                }
        }

        #ifdef CPU_PROFILER
        ProfilerStop();
        #endif

        #ifdef HEAP_PROFILER
        HeapProfilerStop();
        #endif

        //cout << endl << "Section time = " << TOTAL_TIME << endl;

        return optimal;
}
