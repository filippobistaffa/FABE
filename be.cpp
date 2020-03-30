#include "be.hpp"

#define REDUCTION_OP(X, Y) (min(X, Y))
#define REDUCTION_INIT (numeric_limits<value>::max())

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
        set_intersection(v1.begin(),v1.end(), v2.begin(),v2.end(), back_inserter(shared));
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

        for (auto r1 : a1.rows) {
                for (auto r2 : a2.rows) {
                        const value join_val = JOIN_OP(r1.first, r2.first);
                        for (size_t i = 0; i < n_comb; ++i) {
                                struct fa *faj;
                                auto fa1 = fa_clone(r1.second);
                                auto fa2 = fa_clone(r2.second);
                                unordered_map<value, struct fa *>::iterator it;
                                auto comb = get_combination(i, sd);
                                char* re = new char[shared.size() + 1];
                                re[shared.size()] = 0;
                                for (size_t v = 0; v < shared.size(); ++v) {
                                        re[v] = ALPHABET[comb[v]];
                                        fa_filter_letter(fa1, pos[v].first, ALPHABET[comb[v]]);
                                        fa_filter_letter(fa2, pos[v].second, ALPHABET[comb[v]]);
                                        if (fa_is_basic(fa1, FA_EMPTY) || fa_is_basic(fa2, FA_EMPTY)) {
                                                goto next_comb;
                                        }
                                        fa_collapse_level(fa1, pos[v].first);
                                        fa_collapse_level(fa2, pos[v].second);
                                        //fa_make_dot(fa1, "dot/a1-val=%.0f-%s=%s.dot", r1.first, vea2str(shared).c_str(), vea2str(comb).c_str());
                                        //fa_make_dot(fa2, "dot/a2-val=%.0f-%s=%s.dot", r2.first, vea2str(shared).c_str(), vea2str(comb).c_str());
                                }
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

automata join_bucket(vector<automata> const &bucket, vector<size_t> domains) {

        auto res = bucket.front();

	for (auto it = next(bucket.begin()); it != bucket.end(); ++it) {
	        auto old = res;
	        res = join(old, *it, domains);
	        free_automata(old);
	}

        return res;
}

/* void add_intersect_minus(vector<row> const &rows, boost::dynamic_bitset<> p, boost::dynamic_bitset<> n, vector<row> &out) {

        row res = {
                vector<size_t>(),
                fa_make_basic(FA_TOTAL),
                REDUCTION_INIT
        };

        for EACH_SET_BIT(p, i) {
                auto tmp = res.fa;
                res.fa = fa_intersect(res.fa, rows[i].fa);
                fa_free(tmp);
                if (fa_is_basic(res.fa, FA_EMPTY)) {
                        return;
                } else {
                        res.v = REDUCTION_OP(res.v, rows[i].v);
                }
        }

        for EACH_SET_BIT(n, i) {
                auto tmp = res.fa;
                res.fa = fa_minus(res.fa, rows[i].fa);
                fa_free(tmp);
                if (fa_is_basic(res.fa, FA_EMPTY)) {
                        return;
                }
        }

        out.push_back(res);
}*/

/*void reduce_var(cost const &c, size_t var) {

        for (auto row : c.rows) {
                fa_collapse_level(row.fa, var);
        }

        if (c.rows.size() > 1) {
                vector<row> tmp;
                const size_t n = c.rows.size();
                const size_t ceil_div2 = 1 + ((n - 1) / 2);
                boost::dynamic_bitset<> p(n);
                for (auto k = n; k >= ceil_div2; --k) {
                        p.reset();
                        for (auto j = 0; j < k; ++j) {
                                p.set(j);
                        }
                        while (true) {
                                boost::dynamic_bitset<> n_p(p);
                                n_p.flip();
                                add_intersect_minus(c.rows, p, n_p, tmp);
                                add_intersect_minus(c.rows, n_p, p, tmp);
                                if (nth_bit(p, 0) == n - k)
                                        break;
                                int i;
                                for (i = k - 1; i >= 0 && nth_bit(p, i) + k - i == n; --i);
                                auto r = nth_bit(p, i);
                                p.reset(r);
                                p.set(r + 1);
                                int j = 2;
                                for (++i; i < k; ++i, ++j) {
                                        p.reset(nth_bit(p, i));
                                        p.set(r + j);
                                }
                        }
                }
        }
}*/
