#ifndef EXPORT_HPP_
#define EXPORT_HPP_

struct mbe_data {
        char *filename;
        vector<uchar> header;
        vector<vector<uchar>> augmented;
        vector<size_t> n_augmented;
        vector<vector<uchar>> intermediate;
        vector<size_t> n_intermediate;
        vector<size_t> evid_offset;
        size_t cur;
        vector<vector<size_t>> anc;
};

template <typename T>
static inline void buf_push_back(vector<uchar> &buf, T x) {

        uchar *bytes = reinterpret_cast<uchar*>(&x);
        buf.insert(buf.end(), bytes, bytes + sizeof(T));
}

#include <fstream>      // ofstream

void inline write_binary(mbe_data &mbe, double optimal, int ibound) {

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
        // final augmented bucket
        int root_id = n_vars - 1;
        buf_push_back(mbe.augmented.back(), 1ULL);
        buf_push_back(mbe.augmented.back(), -root_id);
        buf_push_back(mbe.augmented.back(), 0ULL);
        buf_push_back(mbe.augmented.back(), 1ULL);
        buf_push_back(mbe.augmented.back(), optimal);
        ofs.write((char*)mbe.augmented.back().data(), mbe.augmented.back().size());
        // close file
        ofs.close();
}

#endif /* LOG_HPP_ */
