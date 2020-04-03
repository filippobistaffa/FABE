#include "be.hpp"

//#define PRINT_TABLES

#define REDUCTION_MIN
//#define REDUCTION_MAX

#define JOIN_OP(X, Y) (X + Y)
#define OP_FREE_OLD(OP, FREE, X, Y) { auto __TMP = (X); (X) = OP(X, Y); FREE(__TMP); }
#define SET_OP(OP, X, Y, R, CMP) (OP((X).begin(), (X).end(), (Y).begin(), (Y).end(), inserter((R), (R).begin()), CMP))

// measure sections of code
/*
#include <chrono>
static chrono::duration<double> total_t;
#define START_CLOCK auto start_t = chrono::high_resolution_clock::now()
#define STOP_CLOCK total_t += chrono::high_resolution_clock::now() - start_t
#define TOTAL_TIME total_t.count()
*/

static automata join(automata &a1, automata &a2, vector<size_t> const &pos, vector<size_t> const &domains) {

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

static automata join_bucket(vector<automata> &bucket, vector<size_t> const &pos, vector<size_t> const &domains) {

        auto res = bucket.front();

	for (auto it = next(bucket.begin()); it != bucket.end(); ++it) {
	        auto old = res;
	        res = join(old, *it, pos, domains);
	}

        return res;
}

static value reduce_last_var(automata &a) {

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

        vector<value> empty;

        for (auto it = next(keys.begin()); it != keys.end(); ++it) {
                for (auto prev = keys.begin(); prev != it; ++prev) {
                        if (!fa_is_basic(a.rows[*prev], FA_EMPTY)) {
                                OP_FREE_OLD(fa_minus, fa_free, a.rows[*it], a.rows[*prev]);
                        }
                }
                if (fa_is_basic(a.rows[*it], FA_EMPTY)) {
                        empty.push_back(*it);
                }
        }

        for (auto e : empty) {
                fa_free(a.rows[e]);
                a.rows.erase(e);
        }

        #ifdef PRINT_TABLES
        print_table(compute_table(a));
        cout << endl;
        #endif
        return 0;
}

value bucket_elimination(vector<vector<automata>> &buckets, vector<size_t> const &order,
                         vector<size_t> const &pos, vector<size_t> const &domains, size_t max_iter) {

        #ifdef PROFILE
        ProfilerStart(PROFILE);
        #endif

        value optimal = 0;

        for (auto it = order.rbegin(); it != order.rend(); ++it) {
                if (buckets[*it].size()) {
                        cout << "Processing bucket " << *it << " with " << buckets[*it].size() << " functions" << endl;
                        auto h = join_bucket(buckets[*it], pos, domains);
                        //automata_dot(h, "dot");
                        optimal += reduce_last_var(h);
                        if (h.vars.size() > 0) {
                                //automata_dot(h, "dot");
                                cout << "Result placed in bucket " << push_bucket(h, buckets, pos) << endl << endl;
                        }
                }
                if (--max_iter == 0) {
                        return optimal;
                }
        }

        #ifdef PROFILE
        ProfilerStop();
        #endif

        //cout << endl << "Section time = " << TOTAL_TIME << endl;

        return optimal;
}
