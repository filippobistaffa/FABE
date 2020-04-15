#ifndef BITSET_H_
#define BITSET_H_

#include <boost/dynamic_bitset.hpp>

using namespace std;

#define GET_MACRO(_1,_2,_3,NAME,...) NAME
#define EACH_SET_BIT(...) GET_MACRO(__VA_ARGS__, EACH_SET_BIT_3, EACH_SET_BIT_2)(__VA_ARGS__)
#define EACH_SET_BIT_2(B, I) (auto I = (B).find_first(); I != boost::dynamic_bitset<>::npos; I = (B).find_next(I))
#define EACH_SET_BIT_3(B, I, S) (auto I = (B).find_next(S); I != boost::dynamic_bitset<>::npos; I = (B).find_next(I))

/*__attribute__((always_inline)) inline
size_t nth_bit(boost::dynamic_bitset<> const &b, int i) {

        size_t ret = b.find_first();
        while (--i >= 0) {
                ret = b.find_next(ret);
        }
        return ret;
}*/

#endif /* BITSET_H_ */
