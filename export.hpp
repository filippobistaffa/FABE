#ifndef EXPORT_HPP_
#define EXPORT_HPP_

#include "conversion.hpp"
#include "io.hpp"

struct bin_data {
        char *filename;
        vector<uchar> header;
        vector<vector<uchar>> augmented;
        vector<size_t> n_augmented;
        vector<vector<pair<size_t, size_t>>> intermediate;
        vector<size_t> n_intermediate;
        vector<size_t> evid_offset;
        vector<vector<size_t>> anc;
};

template <typename T>
static inline void buf_push_back(vector<uchar> &buf, T x) {

        uchar *bytes = reinterpret_cast<uchar*>(&x);
        buf.insert(buf.end(), bytes, bytes + sizeof(T));
}

static inline void export_function(automata &h, size_t from, size_t to, struct bin_data &mbe) {

        //cout << "writing table originated from bucket " << from << " in augmented bucket " << to << endl;
        //print_table(compute_table(h));
        //cout << "ancestors of " << from << " : " << vec2str(mbe.anc[from]) << endl;

        for (auto it = mbe.anc[from].rbegin(); it != mbe.anc[from].rend() && *it != to; ++it) {
                //cout << "writing table in intermediate bucket " << *it << endl;
                mbe.intermediate[*it].push_back(make_pair(to, mbe.n_augmented[to]));
                mbe.n_intermediate[*it]++;
        }

        mbe.n_augmented[to]++;
        // function ID
        buf_push_back(mbe.augmented[to], (int)(-from));
        // scope size
        buf_push_back(mbe.augmented[to], (size_t)h.vars.size());
        // scope
        for (auto it = h.vars.rbegin(); it != h.vars.rend(); ++it) {
                buf_push_back(mbe.augmented[to], (int)(*it - mbe.evid_offset[*it]));
        }
        auto [ cpt, n ] = compute_cpt(h);
        // table size
        buf_push_back(mbe.augmented[to], (size_t)n);
        // table values
        uchar *cpt_bytes = reinterpret_cast<uchar*>(cpt);
        mbe.augmented[to].insert(mbe.augmented[to].end(), cpt_bytes, cpt_bytes + sizeof(double) * n);
        delete[] cpt;
}

#include <fstream>      // ofstream

void inline write_binary(bin_data &mbe, double optimal, int ibound) {

        size_t n_vars = mbe.anc.size();
        // open file
        ofstream ofs(mbe.filename, ios::out | ios::binary);
        // header
        buf_push_back(mbe.header, mbe.augmented.size());
        buf_push_back(mbe.header, ibound);
        buf_push_back(mbe.header, optimal);
        ofs.write((char*)mbe.header.data(), mbe.header.size());
        // augmented buckets
        for (size_t var = 0; var < n_vars; ++var) {
                ofs.write((char*)&(mbe.n_augmented[var]), sizeof(size_t));
                ofs.write((char*)mbe.augmented[var].data(), mbe.augmented[var].size());
        }
        // augmented bucket for dummy root
        int root_id = n_vars - 1;
        buf_push_back(mbe.augmented.back(), 1ULL);
        buf_push_back(mbe.augmented.back(), -root_id);
        buf_push_back(mbe.augmented.back(), 0ULL);
        buf_push_back(mbe.augmented.back(), 1ULL);
        buf_push_back(mbe.augmented.back(), optimal);
        ofs.write((char*)mbe.augmented.back().data(), mbe.augmented.back().size());
        //ofs.close();
        //ofs = ofstream("inter.mbe", ios::out | ios::binary);
        // intermediate buckets
        vector<size_t> pfx(mbe.n_augmented.size());
        exclusive_scan(mbe.n_augmented.begin(), mbe.n_augmented.end(), pfx.begin(), 0, plus<>{});
        for (size_t var = 0; var < n_vars; ++var) {
                ofs.write((char*)&(mbe.n_intermediate[var]), sizeof(size_t));
                for (auto p : mbe.intermediate[var]) {
                        size_t idx = pfx[p.first] + p.second;
                        ofs.write((char*)&idx, sizeof(size_t));
                }
        }
        // intermediate bucket for dummy root
        size_t z = 0LL;
        ofs.write((char*)&z, sizeof(size_t));
        // close file
        ofs.close();
}

#endif /* LOG_HPP_ */
