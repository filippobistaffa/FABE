#ifndef EXPORT_HPP_
#define EXPORT_HPP_

struct bin_data {
        char *filename;
        vector<uchar> header;
        vector<vector<uchar>> augmented;
        vector<size_t> n_augmented;
        vector<vector<pair<size_t, size_t>>> intermediate;
        vector<size_t> n_intermediate;
        vector<size_t> evid_offset;
        vector<vector<size_t>> anc;
        double threshold;
};

__attribute__((always_inline)) inline
void alloc_bin_data(bin_data &mbe, size_t n_vars) {

        mbe.augmented = vector<vector<uchar>>(n_vars + 1, vector<uchar>());
        mbe.n_augmented = vector<size_t>(n_vars + 1);
        mbe.intermediate = vector<vector<pair<size_t, size_t>>>(n_vars, vector<pair<size_t, size_t>>());
        mbe.n_intermediate = vector<size_t>(n_vars);
        mbe.evid_offset = vector<size_t>(n_vars);
}

#include <algorithm>    // fill

__attribute__((always_inline)) inline
pair<double *, size_t> compute_cpt(automata const &a, double def) {

        const size_t n = accumulate(a.domains.begin(), a.domains.end(), 1, multiplies<size_t>());
        double *cpt = new double[n];
        fill(cpt, cpt + n, def);

        for (auto const &[ v, fa ] : a.rows) {
                fa_fill_table(fa, &a.domains[0], cpt, v);
        }

        return make_pair(cpt, n);
}

template <typename T>
__attribute__((always_inline)) inline
void buf_push_back(vector<uchar> &buf, T x) {

        uchar *bytes = reinterpret_cast<uchar*>(&x);
        buf.insert(buf.end(), bytes, bytes + sizeof(T));
}

#include "conversion.hpp"
#include "io.hpp"

__attribute__((always_inline)) inline
void export_function(automata &h, int from, int to, struct bin_data &mbe) {

        //cout << "Writing table from bucket " << from << " in augmented bucket " << to << endl;
        //print_table(compute_table(h));
        //cout << "Ancestors of " << from << " : " << vec2str(mbe.anc[from]) << endl;

        for (auto it = mbe.anc[from].rbegin(); it != mbe.anc[from].rend() && *it != to; ++it) {
                //cout << "Writing table in intermediate bucket " << *it << endl;
                mbe.intermediate[*it].push_back(make_pair(to, mbe.n_augmented[to]));
                mbe.n_intermediate[*it]++;
        }

        mbe.n_augmented[to]++;
        // function ID
        buf_push_back(mbe.augmented[to], -from);
        // scope size
        buf_push_back(mbe.augmented[to], (size_t)h.vars.size());
        // scope
        for (auto it = h.vars.rbegin(); it != h.vars.rend(); ++it) {
                buf_push_back(mbe.augmented[to], (int)(*it - mbe.evid_offset[*it]));
        }
        auto [ cpt, n ] = compute_cpt(h, mbe.threshold);
        // table size
        buf_push_back(mbe.augmented[to], (size_t)n);
        // table values
        uchar *cpt_bytes = reinterpret_cast<uchar*>(cpt);
        mbe.augmented[to].insert(mbe.augmented[to].end(), cpt_bytes, cpt_bytes + sizeof(double) * n);
        delete[] cpt;
}

__attribute__((always_inline)) inline
void export_root_function(double value, int from, struct bin_data &mbe) {

        const size_t to = mbe.augmented.size() - 1;
        //cout << "Writing value " << value << " from bucket " << from << " in root augmented bucket" << endl;
        //cout << "Ancestors of " << from << " : " << vec2str(mbe.anc[from]) << endl;

        for (auto it = mbe.anc[from].rbegin(); it != mbe.anc[from].rend() && *it != to; ++it) {
                //cout << "Writing table in intermediate bucket " << *it << endl;
                mbe.intermediate[*it].push_back(make_pair(to, mbe.n_augmented.back()));
                mbe.n_intermediate[*it]++;
        }

        buf_push_back(mbe.augmented.back(), -from);
        buf_push_back(mbe.augmented.back(), 0ULL);
        buf_push_back(mbe.augmented.back(), 1ULL);
        buf_push_back(mbe.augmented.back(), value);
        mbe.n_augmented.back()++;
}

#include <fstream>

__attribute__((always_inline)) inline
void write_binary(bin_data &mbe, double optimal, int ibound, int root) {

        size_t n_vars = mbe.anc.size();
        // open file
        ofstream ofs(mbe.filename, ios::out | ios::binary);
        // header
        buf_push_back(mbe.header, mbe.augmented.size());
        buf_push_back(mbe.header, ibound);
        buf_push_back(mbe.header, optimal);
        ofs.write((char*)mbe.header.data(), mbe.header.size());
        // augmented buckets
        for (size_t var = 0; var < n_vars + 1; ++var) {
                ofs.write((char*)&(mbe.n_augmented[var]), sizeof(size_t));
                ofs.write((char*)mbe.augmented[var].data(), mbe.augmented[var].size());
        }
        // split bin file
        //ofs.close();
        //ofs = ofstream("../fabe_i.bin", ios::out | ios::binary);
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
