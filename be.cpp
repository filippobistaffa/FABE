#include "be.hpp"

#define REDUCTION_MIN
//#define REDUCTION_MAX

#define JOIN_OP(X, Y) (X + Y)
#define OP_FREE_OLD(OP, FREE, X, Y) { auto __TMP = (X); (X) = OP(X, Y); FREE(__TMP); }
#define SET_OP(OP, X, Y, R, CMP) (OP((X).begin(), (X).end(), (Y).begin(), (Y).end(), inserter((R), (R).begin()), CMP))

static automata join(automata &a1, automata &a2, vector<size_t> const &pos, vector<size_t> const &domains) {

        #ifdef PRINT_TABLES
        cout << endl << "Joining:" << endl << endl;
        print_table(compute_table(a1));
        cout << endl;
        print_table(compute_table(a2));
        cout << endl;
        #endif

        automata join;
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

        cout << "Cloning first function..." << endl;
        auto a1c = clone_automata(a1);
        cout << "Cloning second function..." << endl;
        auto a2c = clone_automata(a2);
        a1c.vars = a2c.vars = join.vars;
        a1c.domains = a2c.domains = join.domains;
        cout << "Filling levels..." << endl;

        for (auto &[ v, fa ] : a1c.rows) {
                for (auto i = 0; i < add1.size(); ++i) {
                        fa_add_level(fa, padd1[i], ALPHABET[domains[add1[i]] - 1]);
                }
        }

        for (auto &[ v, fa ] : a2c.rows) {
                for (auto i = 0; i < add2.size(); ++i) {
                        fa_add_level(fa, padd2[i], ALPHABET[domains[add2[i]] - 1]);
                }
        }

        //print_table(compute_table(a1c));
        //print_table(compute_table(a2c));
        cout << "Joining..." << endl;

        for (auto &[ v1, fa1 ] : a1c.rows) {
                for (auto &[ v2, fa2 ] : a2c.rows) {
                        auto in = fa_intersect(fa1, fa2);
                        if (fa_is_basic(in, FA_EMPTY)) { // these two rows do not have any common variable
                                fa_free(in);             // assignment of shared variables
                        } else {
                                //fa_make_dot(fa1, "dot/v1=%.0f.dot", v1);
                                //fa_make_dot(fa2, "dot/v2=%.0f.dot", v2);
                                const value v1v2 = JOIN_OP(v1, v2);
                                auto it = join.rows.find(v1v2);
                                if (it == join.rows.end()) {
                                        join.rows.insert({ v1v2, in });
                                        //fa_make_dot(in, "dot/r-v1v2=%.0f.dot", v1v2);
                                } else {
                                        OP_FREE_OLD(fa_union, fa_free, it->second, in);
                                        fa_minimize(it->second);
                                        //fa_make_dot(it->second, "dot/v1v2=%.0f.dot", v1v2);
                                }
                        }
                }
        }

        cout << "Done joining." << endl;
        //print_table(compute_table(join));
        return join;
}

static automata join_bucket(vector<automata> &bucket, vector<size_t> const &pos, vector<size_t> const &domains) {

        auto res = clone_automata(bucket.front());

	for (auto it = next(bucket.begin()); it != bucket.end(); ++it) {
	        auto old = res;
	        res = join(old, *it, pos, domains);
	}

        cout << "Joined all functions." << endl;
        return res;
}

static value reduce_last_var(automata &a) {

        #ifdef PRINT_TABLES
        cout << endl;
        print_table(compute_table(a));
        cout << endl;
        #endif

        // level to remove
        const auto i = a.vars.size() - 1;

        // remove last variable and domain
        a.vars.pop_back();
        a.domains.pop_back();
        vector<value> keys;
        cout << "Collapsing levels..." << endl;

        for (auto &[ v, fa ] : a.rows) {
                fa_collapse_level(fa, i);
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
        cout << "Minimizing..." << endl;

        for (auto it = next(keys.begin()); it != keys.end(); ++it) {
                for (auto prev = keys.begin(); prev != it; ++prev) {
                        OP_FREE_OLD(fa_minus, fa_free, a.rows[*it], a.rows[*prev]);
                }
                //fa_minimize(a.rows[*it]);
                if (fa_is_basic(a.rows[*it], FA_EMPTY)) {
                        empty.push_back(*it);
                }
        }

        for (auto e : empty) {
                a.rows.erase(e);
        }

        cout << "Done minimizing." << endl;
        #ifdef PRINT_TABLES
        print_table(compute_table(a));
        cout << endl;
        #endif
        return 0;
}

value bucket_elimination(vector<vector<automata>> &buckets, vector<size_t> const &order,
                         vector<size_t> const &pos, vector<size_t> const &domains, size_t max_iter) {

        value optimal = 0;

        for (auto it = order.rbegin(); it != order.rend(); ++it) {
                if (buckets[*it].size()) {
                        cout << "Processing bucket " << *it << " with " << buckets[*it].size() << " functions" << endl;
                        auto h = join_bucket(buckets[*it], pos, domains);
                        //automata_dot(h, "dot");
                        optimal += reduce_last_var(h);
                        if (h.vars.size() > 0) {
                                //automata_dot(h, "dot");
                                cout << "Placed in bucket " << push_bucket(h, buckets, pos) << endl << endl;
                        }
                }
                if (--max_iter == 0) {
                        return optimal;
                }
        }

        return optimal;
}
