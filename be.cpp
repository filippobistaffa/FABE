#include "be.hpp"

void reduce_var(cost const &c, size_t var) {

        for (auto row : c.rows) {
                //enumerate_state(fa_state_initial(row.fa), var);
                fa_collapse_level(row.fa, var);
        }
}
