#include "minimise.hpp"

row automata(vector<row>::const_iterator begin, vector<row>::const_iterator end) {

        ostringstream oss;

        for (auto row = begin; row != end; ++row) {
                for (auto i : row->a) {
                        oss << ALPHABET[i];
                }
                oss << "|";
        }

        row ret;
        ret.v = begin->v;
        fa_compile(oss.str().c_str(), oss.str().length() - 1, &(ret.fa));
        fa_minimize(ret.fa);
        return ret;
}

__attribute__((always_inline)) inline
bool are_equal(value a, value b, value tolerance) {

        const value epsilon = numeric_limits<value>::epsilon();
        return fabs(a - b) <= tolerance + epsilon;
};

cost compress_clusters(cost const &in, value tolerance) {

        cost out = {
                in.vars,
                in.domains,
                vector<row>()
        };

        auto begin = in.rows.begin();
        auto end = begin;

        while (begin != in.rows.end()) {
                while (end != in.rows.end() && are_equal(end->v, begin->v, tolerance)) {
                        end++;
                }
                //cout << "begin " << begin - in.rows.begin() << endl;
                //cout << "end " << end - in.rows.begin() << endl;
                out.rows.push_back(automata(begin, end));
                begin = end;
        }

        return out;
}
