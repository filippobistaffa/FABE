#include "minimise.hpp"

static struct fa *fa_compile_minimise(vector<pair<vector<size_t>, value>>::const_iterator begin,
                           vector<pair<vector<size_t>, value>>::const_iterator end) {

        ostringstream oss;

        for (auto row = begin; row != end; ++row) {
                for (auto i : row->first) {
                        oss << ALPHABET[i];
                }
                oss << "|";
        }

        struct fa *fa;
        fa_compile(oss.str().c_str(), oss.str().length() - 1, &fa);
        fa_minimize(fa);
        return fa;
}

__attribute__((always_inline)) inline
bool are_equal(value a, value b, value tolerance) {

        const value epsilon = numeric_limits<value>::epsilon();
        return fabs(a - b) <= tolerance + epsilon;
};

automata compute_automata(table const &t, value tolerance) {

        automata res = {
                t.vars,
                t.domains,
                unordered_map<value, struct fa *>()
        };

        auto begin = t.rows.begin();
        auto end = begin;

        while (begin != t.rows.end()) {
                while (end != t.rows.end() && are_equal(end->second, begin->second, tolerance)) {
                        end++;
                }
                //cout << "begin " << begin - in.rows.begin() << endl;
                //cout << "end " << end - in.rows.begin() << endl;
                res.rows.insert({ begin->second, fa_compile_minimise(begin, end) });
                begin = end;
        }

        return res;
}
