#include "be.hpp"

#define REDUCTION_MIN
//#define REDUCTION_MAX

#define JOIN_OP(X, Y) (X + Y)
#define JOIN_INIT (0)

#define OP_FREE_OLD(OP, FREE, X, Y) { auto __TMP = (X); (X) = OP(X, Y); FREE(__TMP); }

static automata join(automata const &a1, automata const &a2, vector<size_t> domains) {

        automata join;
        auto v1 = vector<size_t>(a1.vars);
        auto v2 = vector<size_t>(a2.vars);
        vector<size_t> shared;
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
                                //cout << i + 1 << "/" << n_comb << endl;
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
}

static automata join_bucket(vector<automata> const &bucket, vector<size_t> domains) {

        auto res = copy_automata(bucket.front());

	for (auto it = next(bucket.begin()); it != bucket.end(); ++it) {
	        auto old = res;
	        res = join(old, *it, domains);
	}

        return res;
}

static value reduce_var(automata &a, size_t var) {

        const auto i = find(a.vars.begin(), a.vars.end(), var) - a.vars.begin();
        a.vars.erase(a.vars.begin() + i);
        a.domains.erase(a.domains.begin() + i);
        vector<value> keys;

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

        return 0;
}

void bucket_elimination(vector<vector<automata>> &buckets, vector<size_t> const &order,
                        vector<size_t> const &pos, vector<size_t> const &domains, size_t max_iter) {

        value optimal = 0;

        for (auto it = order.rbegin(); it != order.rend(); ++it) {
                cout << endl << "Processing bucket " << *it << " with " << buckets[*it].size() << " functions" << endl;
                auto h = join_bucket(buckets[*it], domains);
                automata_dot(h, "dot");
                optimal += reduce_var(h, *it);
                if (h.vars.size() > 0) {
                        automata_dot(h, "dot");
                        cout << "Placed in bucket " << push_bucket(h, buckets, pos) << endl;
                }
                if (--max_iter == 0) {
                        return;
                }
        }

        cout << endl << "Optimal solution = " << optimal << endl;
}
