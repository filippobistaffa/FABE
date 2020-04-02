#include "be.hpp"

#define REDUCTION_MIN
//#define REDUCTION_MAX

#define JOIN_OP(X, Y) (X + Y)

#define OP_FREE_OLD(OP, FREE, X, Y) { auto __TMP = (X); (X) = OP(X, Y); FREE(__TMP); }

/*static automata join(automata const &a1, automata const &a2, vector<size_t> domains) {

        automata join;

        // compute shared variables
        vector<size_t> shared;
        auto v1 = vector<size_t>(a1.vars);
        auto v2 = vector<size_t>(a2.vars);
        sort(v1.begin(), v1.end());
        sort(v2.begin(), v2.end());
        set_intersection(v1.begin(),v1.end(), v2.begin(), v2.end(), back_inserter(shared));

        vector<pair<size_t, size_t>> pos;
        vector<size_t> sd;

        for (auto var : shared) {
                pos.push_back(make_pair(find(a1.vars.begin(), a1.vars.end(), var) - a1.vars.begin(),
                                        find(a2.vars.begin(), a2.vars.end(), var) - a2.vars.begin()));
                sd.push_back(domains[var]);
        }

        auto nsv1 = vector<size_t>(a1.vars);
        auto nsv2 = vector<size_t>(a2.vars);

        for (auto var : shared) {
                nsv1.erase(find(nsv1.begin(), nsv1.end(), var));
                nsv2.erase(find(nsv2.begin(), nsv2.end(), var));
        }

        join.vars.insert(join.vars.end(), shared.begin(), shared.end());
        join.vars.insert(join.vars.end(), nsv1.begin(), nsv1.end());
        join.vars.insert(join.vars.end(), nsv2.begin(), nsv2.end());

        for (auto var : join.vars) {
                join.domains.push_back(domains[var]);
        }

        const auto n_comb = accumulate(sd.begin(), sd.end(), 1, multiplies<size_t>());
        cout << vec2str(join.vars) << endl;
        cout << vec2str(shared) << " " << vec2str(sd) << " -> " << n_comb << endl;

        for (auto const &r1 : a1.rows) {
                for (auto const &r2 : a2.rows) {
                        const value join_val = JOIN_OP(r1.first, r2.first);
                        for (size_t i = 0; i < n_comb; ++i) {
                                cout << i + 1 << "/" << n_comb << endl;
                                struct fa *faj;
                                auto fa1 = fa_clone(r1.second);
                                auto fa2 = fa_clone(r2.second);
                                unordered_map<value, struct fa *>::iterator it;
                                auto comb = get_combination(i, sd);
                                char* re = new char[shared.size() + 1];
                                re[shared.size()] = 0;
                                for (int v = shared.size(); v --> 0;) {
                                        re[v] = ALPHABET[comb[v]];
                                        fa_filter_letter(fa1, pos[v].first, ALPHABET[comb[v]]);
                                        fa_filter_letter(fa2, pos[v].second, ALPHABET[comb[v]]);
                                        if (fa_is_basic(fa1, FA_EMPTY) || fa_is_basic(fa2, FA_EMPTY)) {
                                                goto next_comb;
                                        }
                                        fa_collapse_level(fa1, pos[v].first);
                                        fa_collapse_level(fa2, pos[v].second);
                                }
                                //fa_make_dot(fa1, "dot/a1-val=%.0f-%s=%s.dot", r1.first, vec2str(shared).c_str(), vec2str(comb).c_str());
                                //fa_make_dot(fa2, "dot/a2-val=%.0f-%s=%s.dot", r2.first, vec2str(shared).c_str(), vec2str(comb).c_str());
                                fa_compile(re, shared.size(), &faj);
                                OP_FREE_OLD(fa_concat, fa_free, faj, fa1);
                                OP_FREE_OLD(fa_concat, fa_free, faj, fa2);
                                it = join.rows.find(join_val);
                                if (it == join.rows.end()) {
                                        join.rows.insert({ join_val, faj });
                                } else {
                                        OP_FREE_OLD(fa_union, fa_free, it->second, faj);
                                        fa_minimize(it->second);
                                }
                            next_comb:
                                fa_free(fa1);
                                fa_free(fa2);
                                delete[] re;
                        }
                }
        }

        return join;
}*/

/*static automata join(automata &a1, automata &a2, vector<size_t> const &domains) {

        #ifdef PRINT_TABLES
        cout << endl << "Joining:" << endl << endl;
        print_table(compute_table(a1));
        cout << endl;
        print_table(compute_table(a2));
        cout << endl;
        #endif

        // result function
        automata join;

        // shared variables (sv)
        vector<size_t> sv;
        auto v1 = a1.vars;
        auto v2 = a2.vars;
        sort(v1.begin(), v1.end());
        sort(v2.begin(), v2.end());
        set_intersection(v1.begin(),v1.end(), v2.begin(), v2.end(), back_inserter(sv));

        // shared variables' positions (psv1, psv2)
        auto psv1 = vector<size_t>(sv.size());
        auto psv2 = vector<size_t>(sv.size());

        for (auto i = 0; i < sv.size(); ++i) {
                psv1[i] = find(a1.vars.begin(), a1.vars.end(), sv[i]) - a1.vars.begin();
                psv2[i] = find(a2.vars.begin(), a2.vars.end(), sv[i]) - a2.vars.begin();
        }

        // non-shared variables (nsv1, nsv2)
        auto nsv1 = a1.vars;
        auto nsv2 = a2.vars;

        for (auto var : sv) {
                nsv1.erase(find(nsv1.begin(), nsv1.end(), var));
                nsv2.erase(find(nsv2.begin(), nsv2.end(), var));
        }

        // non-shared variables' positions (pnsv1, pnsv2)
        auto pnsv1 = vector<size_t>(nsv1.size());
        auto pnsv2 = vector<size_t>(nsv2.size());

        for (auto i = 0; i < nsv1.size(); ++i) {
                pnsv1[i] = find(a1.vars.begin(), a1.vars.end(), nsv1[i]) - a1.vars.begin();
        }

        for (auto i = 0; i < nsv2.size(); ++i) {
                pnsv2[i] = find(a2.vars.begin(), a2.vars.end(), nsv2[i]) - a2.vars.begin();
        }

        sort(psv1.begin(), psv1.end());
        sort(psv2.begin(), psv2.end());
        sort(pnsv1.begin(), pnsv1.end());
        sort(pnsv2.begin(), pnsv2.end());

        cout << vec2str(a1.vars, "v1") << endl;
        cout << vec2str(a2.vars, "v2") << endl;
        cout << vec2str(sv, "sv") << endl;
        cout << vec2str(psv1, "psv1") << endl;
        cout << vec2str(psv2, "psv2") << endl;
        cout << vec2str(nsv1, "nsv1") << endl;
        cout << vec2str(nsv2, "nsv2") << endl;
        cout << vec2str(pnsv1, "pnsv1") << endl;
        cout << vec2str(pnsv2, "pnsv2") << endl;

        // join function's variables and domains
        join.vars.insert(join.vars.end(), a1.vars.begin(), a1.vars.end());
        join.vars.insert(join.vars.end(), nsv2.begin(), nsv2.end());

        for (auto var : join.vars) {
                join.domains.push_back(domains[var]);
        }

        // functions with all non-shared variables filtered out (a1c, a2c)
        auto a1c = clone_automata(a1);
        auto a2c = clone_automata(a2);

        for (auto it = pnsv1.rbegin(); it != pnsv1.rend(); ++it) {
                for (auto &[ v, fa ] : a1c.rows) {
                        fa_collapse_level(fa, *it);
                }
        }

        for (auto it = pnsv2.rbegin(); it != pnsv2.rend(); ++it) {
                for (auto &[ v, fa ] : a2c.rows) {
                        fa_collapse_level(fa, *it);
                }
        }

        for (auto &[ v1, fa1 ] : a1c.rows) {
                //cout << v1 << endl;
                for (auto &[ v2, fa2 ] : a2c.rows) {
                        if (v1 + v2 == 15)
                                BREAKPOINT("DELETE DOTS");
                        //cout << v2 << endl;
                        auto in = fa_intersect(fa1, fa2);
                        if (fa_is_basic(in, FA_EMPTY)) { // these two rows do not have any common variable
                                fa_free(in);             // assignment of shared variables
                        } else {
                                fa_make_dot(a1.rows[v1], "dot/0-a1-v1=%.0f.dot", v1);
                                fa_make_dot(a2.rows[v2], "dot/0-a2-v2=%.0f.dot", v2);
                                // save intersection automata for later
                                struct fa *in2;
                                if (nsv2.size() > 0) {
                                        in2 = fa_clone(in);
                                }
                                // if "a1" has extra variables, remove the paths that don't go through
                                // any common assignment of shared variables 
                                if (nsv1.size() > 0) {
                                        fa_make_dot(fa1, "dot/1-v1=%.0f-v2=%.0f-2.dot", v1, v2);
                                        fa_make_dot(in, "dot/1-v1=%.0f-v2=%.0f-3.dot", v1, v2);
                                        for (auto i = 0; i < nsv1.size(); ++i) {
                                                fa_add_level(in, pnsv1[i], ALPHABET[domains[nsv1[i]] - 1]);
                                        }
                                        fa_make_dot(in, "dot/1-v1=%.0f-v2=%.0f-4.dot", v1, v2);
                                        OP_FREE_OLD(fa_intersect, fa_free, in, a1.rows[v1]);
                                        fa_make_dot(in, "dot/1-v1=%.0f-v2=%.0f-5.dot", v1, v2);
                                }
                                if (nsv2.size() > 0) {
                                        fa_make_dot(fa2, "dot/2-v1=%.0f-v2=%.0f-2.dot", v1, v2);
                                        fa_make_dot(in2, "dot/2-v1=%.0f-v2=%.0f-3.dot", v1, v2);
                                        for (auto i = 0; i < nsv2.size(); ++i) {
                                                fa_add_level(in2, pnsv2[i], ALPHABET[domains[nsv2[i]] - 1]);
                                        }
                                        fa_make_dot(in2, "dot/2-v1=%.0f-v2=%.0f-4.dot", v1, v2);
                                        OP_FREE_OLD(fa_intersect, fa_free, in2, a2.rows[v2]);
                                        fa_make_dot(in2, "dot/2-v1=%.0f-v2=%.0f-5.dot", v1, v2);
                                        // remove shared variables from second automata (to append it)
                                        for (auto it = psv2.rbegin(); it != psv2.rend(); ++it) {
                                                fa_collapse_level(in2, *it);
                                                fa_make_dot(in2, "dot/collapsed=%zu.dot", *it);
                                        }
                                        fa_make_dot(in2, "dot/2-v1=%.0f-v2=%.0f-6.dot", v1, v2);
                                        OP_FREE_OLD(fa_concat, fa_free, in, in2);
                                        fa_make_dot(in, "dot/2-v1=%.0f-v2=%.0f-7.dot", v1, v2);
                                }
                                const value v1v2 = JOIN_OP(v1, v2);
                                auto it = join.rows.find(v1v2);
                                if (it == join.rows.end()) {
                                        join.rows.insert({ v1v2, in });
                                        fa_make_dot(in, "dot/r-v1v2=%.0f.dot", v1v2);
                                } else {
                                        OP_FREE_OLD(fa_union, fa_free, it->second, in);
                                        fa_minimize(it->second);
                                        fa_make_dot(it->second, "dot/r-v1v2=%.0f.dot", v1v2);
                                }
                                if (v1 + v2 == 15)
                                        BREAKPOINT("CHECK DOTS");
                        }
                }
        }

        #ifdef PRINT_TABLES
        cout << endl;
        print_table(compute_table(join));
        #endif

        return join;
}*/

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

        cout << vec2str(a1.vars, "v1") << endl;
        cout << vec2str(a2.vars, "v2") << endl;
        cout << vec2str(join.vars, "vj") << endl;
        cout << vec2str(add1, "add1") << endl;
        cout << vec2str(add2, "add2") << endl;
        cout << vec2str(padd1, "padd1") << endl;
        cout << vec2str(padd2, "padd2") << endl;

        cout << "cloning1" << endl;

        auto a1c = clone_automata(a1);
        auto a2c = clone_automata(a2);
        a1c.vars = a2c.vars = join.vars;
        a1c.domains = a2c.domains = join.domains;

        cout << "filling levels" << endl;

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

        cout << "joining" << endl;
        //print_table(compute_table(a1c));
        //print_table(compute_table(a2c));

        for (auto &[ v1, fa1 ] : a1c.rows) {
                for (auto &[ v2, fa2 ] : a2c.rows) {
                        auto in = fa_intersect(fa1, fa2);
                        if (fa_is_basic(in, FA_EMPTY)) { // these two rows do not have any common variable
                                fa_free(in);             // assignment of shared variables
                        } else {
                                fa_make_dot(fa1, "dot/v1=%.0f.dot", v1);
                                fa_make_dot(fa2, "dot/v2=%.0f.dot", v2);
                                const value v1v2 = JOIN_OP(v1, v2);
                                auto it = join.rows.find(v1v2);
                                if (it == join.rows.end()) {
                                        join.rows.insert({ v1v2, in });
                                        //fa_make_dot(in, "dot/r-v1v2=%.0f.dot", v1v2);
                                } else {
                                        OP_FREE_OLD(fa_union, fa_free, it->second, in);
                                        fa_minimize(it->second);
                                        fa_make_dot(it->second, "dot/v1v2=%.0f.dot", v1v2);
                                }
                        }
                }
        }

        cout << "done" << endl;

        //print_table(compute_table(join));

        //BREAKPOINT("");
        
        return join;
}

static automata join_bucket(vector<automata> &bucket, vector<size_t> const &pos, vector<size_t> const &domains) {

        auto res = clone_automata(bucket.front());

	for (auto it = next(bucket.begin()); it != bucket.end(); ++it) {
	        auto old = res;
	        res = join(old, *it, pos, domains);
	}

        cout << "added all" << endl;

        return res;
}

static value reduce_var(automata &a, size_t var) {

        #ifdef PRINT_TABLES
        cout << endl;
        cout << "Reducing var " << var << ":" << endl;
        cout << endl;
        print_table(compute_table(a));
        cout << endl;
        #endif

        const auto i = find(a.vars.begin(), a.vars.end(), var) - a.vars.begin();
        a.vars.erase(a.vars.begin() + i);
        a.domains.erase(a.domains.begin() + i);
        vector<value> keys;

        cout << "collapsing" << endl;

        for (auto &[ v, fa ] : a.rows) {
                //if (v == 17)
                //        BREAKPOINT("DELETE DOTS");
                fa_collapse_level(fa, i);
                //if (v == 17)
                //        BREAKPOINT("CHECK DOTS");
                keys.push_back(v);
        }

        cout << "done" << endl;

        //automata_dot(a, "dot");

        #ifdef REDUCTION_MIN
        sort(keys.begin(), keys.end());
        #else
        sort(keys.begin(), keys.end(), greater<value>());
        #endif

        if (a.vars.size() == 0) {
                return keys.front();
        }

        vector<value> empty;
        cout << "minimizing" << endl;

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

        cout << "done" << endl;

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
                        //if (max_iter == 1) {
                        //        BREAKPOINT("DELETE DOTS");
                        //}
                        cout << "Processing bucket " << *it << " with " << buckets[*it].size() << " functions" << endl;
                        //for (auto a : buckets[*it])
                        //        print_table(compute_table(a));
                        auto h = join_bucket(buckets[*it], pos, domains);
                        cout << "now reduce" << endl;
                        //automata_dot(h, "dot");
                        optimal += reduce_var(h, *it);
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
