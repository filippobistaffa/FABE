#include "be.hpp"

#define REDUCTION_OP(X, Y) (min(X, Y))
#define REDUCTION_INIT (numeric_limits<value>::max())

#define JOIN_OP(X, Y) (X + Y)
#define JOIN_INIT (0)

#define OP_FREE_OLD(OP, FREE, X, Y) { auto __TMP = (X); (X) = OP(X, Y); FREE(__TMP); }

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
} */

automata join(automata const &c1, automata const &c2, vector<size_t> domains) {

        automata join;
        auto v1 = vector<size_t>(c1.vars);
        auto v2 = vector<size_t>(c2.vars);
        vector<size_t> shared;
        sort(v1.begin(), v1.end());
        sort(v2.begin(), v2.end());
        set_intersection(v1.begin(),v1.end(), v2.begin(),v2.end(), back_inserter(shared));
        vector<pair<size_t, size_t>> shared_pos;
        vector<size_t> shared_dom;

        for (auto var : shared) {
                shared_pos.push_back(make_pair(find(c1.vars.begin(), c1.vars.end(), var) - c1.vars.begin(),
                                               find(c2.vars.begin(), c2.vars.end(), var) - c2.vars.begin()));
                shared_dom.push_back(domains[var]);
        }

        auto nsv1 = vector<size_t>(c1.vars);
        auto nsv2 = vector<size_t>(c2.vars);

        for (auto it = shared_pos.rbegin(); it != shared_pos.rend(); ++it) {
                nsv1.erase(nsv1.begin() + it->first);
                nsv2.erase(nsv2.begin() + it->second);
        }

        join.vars.insert(join.vars.end(), shared.begin(), shared.end());
        join.vars.insert(join.vars.end(), nsv1.begin(), nsv1.end());
        join.vars.insert(join.vars.end(), nsv2.begin(), nsv2.end());

        for (auto var : join.vars) {
                join.domains.push_back(domains[var]);
        }

        const auto n_comb = accumulate(shared_dom.begin(), shared_dom.end(), 1, std::multiplies<size_t>());

        for (auto r1 : c1.rows) {
                for (auto r2 : c2.rows) {
                        const value join_val = JOIN_OP(r1.first, r2.first);
                        for (size_t i = 0; i < n_comb; ++i) {
                                struct fa *faj;
                                auto fa1 = fa_clone(r1.second);
                                auto fa2 = fa_clone(r2.second);
                                unordered_map<value, struct fa *>::iterator it;
                                auto comb = get_combination(i, shared_dom);
                                char* re = new char[shared.size() + 1];
                                re[shared.size()] = 0;
                                for (size_t v = 0; v < shared.size(); ++v) {
                                        re[v] = ALPHABET[comb[v]];
                                        fa_filter_letter(fa1, shared_pos[v].first, ALPHABET[comb[v]]);
                                        fa_filter_letter(fa2, shared_pos[v].second, ALPHABET[comb[v]]);
                                        if (fa_is_basic(fa1, FA_EMPTY) || fa_is_basic(fa2, FA_EMPTY)) {
                                                goto next_comb;
                                        }
                                        fa_collapse_level(fa1, shared_pos[v].first);
                                        fa_collapse_level(fa2, shared_pos[v].second);
                                        fa_make_dot(fa1, "dot/c1-val=%.0f-%s=%s.dot", r1.first, vec2str(shared).c_str(), vec2str(comb).c_str());
                                        fa_make_dot(fa2, "dot/c2-val=%.0f-%s=%s.dot", r2.first, vec2str(shared).c_str(), vec2str(comb).c_str());
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
