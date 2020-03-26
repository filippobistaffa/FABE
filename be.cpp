#include "be.hpp"

void reduce_var(cost const &c, size_t var) {

        for (auto row : c.rows) {
                fa_collapse_level(row.fa, var);
        }

        if (c.rows.size() > 1) {
                vector<row> tmp;
                const size_t n = c.rows.size();
                const size_t ceil_div2 = 1 + ((n - 1) / 2);
                boost::dynamic_bitset<> p(n);
                for (auto k = 5; k >= ceil_div2; --k) {
                        p.reset();
                        for (auto j = 0; j < k; ++j) {
                                p.set(j);
                        }
                        while (true) {
                                //cout << p << endl;
                                p.flip();
                                //cout << p << endl;
                                p.flip();
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
}
