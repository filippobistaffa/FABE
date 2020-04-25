/* Ckmeans_1d_dp.cpp -- Performs 1-D k-means by a dynamic programming
 * approach that is guaranteed to be optimal.
 *
 * Copyright (C) 2010-2016 Mingzhou Song and Haizhou Wang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "kmeans.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <cstring>

//#define DEBUG

inline double ssq
  (const size_t j, const size_t i,
   const vector<double> & sum_x, // running sum of xi
   const vector<double> & sum_x_sq, // running sum of xi^2
   const vector<double> & sum_w=vector<double>(0))  // running sum of weights
{
  double sji(0.0);

  if(sum_w.empty()) { // equally weighted version
    if(j >= i) {
      sji = 0.0;
    } else if(j > 0) {
      double muji = (sum_x[i] - sum_x[j-1]) / (i - j + 1);
      sji = sum_x_sq[i] - sum_x_sq[j-1] - (i - j + 1) * muji * muji;
    } else {
      sji = sum_x_sq[i] - sum_x[i] * sum_x[i] / (i+1);
    }
  } else { // unequally weighted version
    if(sum_w[j] >= sum_w[i]) {
      sji = 0.0;
    } else if(j > 0) {
      double muji = (sum_x[i] - sum_x[j-1]) / (sum_w[i] - sum_w[j-1]);
      sji = sum_x_sq[i] - sum_x_sq[j-1] - (sum_w[i] - sum_w[j-1]) * muji * muji;
    } else {
      sji = sum_x_sq[i] - sum_x[i] * sum_x[i] / sum_w[i];
    }
  }

  sji = (sji < 0) ? 0 : sji;
  return sji;
}

inline double sabs
  (const size_t j, const size_t i,
   const vector<double> & sum_x, // running sum of xi
   const vector<double> & sum_w)  // running sum of weights
{
  double sji(0.0);

  if(sum_w.empty()) { // equally weighted version
    if(j >= i) {
      sji = 0.0;
    } else if(j > 0) {
      size_t l = (i + j) >> 1; // l is the index to the median of the cluster

      if(((i-j+1) % 2) == 1) {
        // If i-j+1 is odd, we have
        //   sum (x_l - x_m) over m = j .. l-1
        //   sum (x_m - x_l) over m = l+1 .. i
        sji = - sum_x[l-1] + sum_x[j-1] + sum_x[i] - sum_x[l];
      } else {
        // If i-j+1 is even, we have
        //   sum (x_l - x_m) over m = j .. l
        //   sum (x_m - x_l) over m = l+1 .. i
        sji = - sum_x[l] + sum_x[j-1] + sum_x[i] - sum_x[l];
      }
    } else { // j==0
      size_t l = i >> 1; // l is the index to the median of the cluster

      if(((i+1) % 2) == 1) {
        // If i-j+1 is odd, we have
        //   sum (x_m - x_l) over m = 0 .. l-1
        //   sum (x_l - x_m) over m = l+1 .. i
        sji = - sum_x[l-1] + sum_x[i] - sum_x[l];
      } else {
        // If i-j+1 is even, we have
        //   sum (x_m - x_l) over m = 0 .. l
        //   sum (x_l - x_m) over m = l+1 .. i
        sji = - sum_x[l] + sum_x[i] - sum_x[l];
      }
    }
  } else { // unequally weighted version
    // no exact solutions are known.
  }

  sji = (sji < 0) ? 0 : sji;
  return sji;
}

inline double dissimilarity
  (const enum DISSIMILARITY dis,
   const size_t j, const size_t i,
   const vector<double> & sum_x, // running sum of xi
   const vector<double> & sum_x_sq, // running sum of xi^2
   const vector<double> & sum_w,  // running sum of weights
   const vector<double> & sum_w_sq)  // running sum of square of weights
{
  double d=0;

  switch(dis) {
  case L1:
    d = sabs(j, i, sum_x, sum_w);
    break;
  case L2:
    d = ssq(j, i, sum_x, sum_x_sq, sum_w);
    break;
  }
  return d;
}

inline void reduce_in_place
  (int imin, int imax, int istep, int q,
   const vector<size_t> & js,
   vector<size_t> & js_red,
   const vector< vector<double> > & S,
   const vector< vector<size_t> > & J,
   const vector<double> & sum_x,
   const vector<double> & sum_x_sq,
   const vector<double> & sum_w,
   const vector<double> & sum_w_sq,
   const enum DISSIMILARITY criterion)
{
  int N = (imax - imin) / istep + 1;

  js_red = js;

  if(N >= js.size()) {
    return;
  }

  // Two positions to move candidate j's back and forth
  int left = -1; // points to last favorable position / column
  int right = 0; // points to current position / column

  size_t m = js_red.size();

  while(m > N) { // js_reduced has more than N positions / columns

    int p = (left + 1);

    int i = (imin + p * istep);
    size_t j = (js_red[right]);
    double Sl = (S[q-1][j-1] +
      dissimilarity(criterion, j, i, sum_x, sum_x_sq, sum_w, sum_w_sq));
      // ssq(j, i, sum_x, sum_x_sq, sum_w));

    size_t jplus1 = (js_red[right+1]);
    double Slplus1 = (S[q-1][jplus1-1] +
      dissimilarity(criterion, jplus1, i, sum_x, sum_x_sq, sum_w, sum_w_sq));
      // ssq(jplus1, i, sum_x, sum_x_sq, sum_w));

    if(Sl < Slplus1 && p < N-1) {
      js_red[ ++ left ] = j; // i += istep;
      right ++; // move on to next position / column p+1
    } else if(Sl < Slplus1 && p == N-1) {
      js_red[ ++ right ] = j; // delete position / column p+1
      m --;
    } else { // (Sl >= Slplus1)
      if(p > 0) { // i > imin
        // delete position / column p and
        //   move back to previous position / column p-1:
        js_red[right] = js_red[left --];
        // p --; // i -= istep;
      } else {
        right ++; // delete position / column 0
      }
      m --;
    }
  }

  for(int r=(left+1); r < m; ++r) {
    js_red[r] = js_red[right++];
  }

  js_red.resize(m);
  return;
}

inline void fill_even_positions
  (int imin, int imax, int istep, int q,
   const vector<size_t> & js,
   vector< vector<double> > & S,
   vector< vector<size_t> > & J,
   const vector<double> & sum_x,
   const vector<double> & sum_x_sq,
   const vector<double> & sum_w,
   const vector<double> & sum_w_sq,
   const enum DISSIMILARITY criterion)
{
  // Derive j for even rows (0-based)
  size_t n = (js.size());
  int istepx2 = (istep << 1);
  size_t jl = (js[0]);
  for(int i=(imin), r(0); i<=imax; i+=istepx2) {

    // auto jmin = (i == imin) ? js[0] : J[q][i - istep];

    while(js[r] < jl) {
      // Increase r until it points to a value of at least jmin
      r ++;
    }

    // Initialize S[q][i] and J[q][i]
    S[q][i] = S[q-1][js[r]-1] +
      dissimilarity(criterion, js[r], i, sum_x, sum_x_sq, sum_w, sum_w_sq);
      // ssq(js[r], i, sum_x, sum_x_sq, sum_w);
    J[q][i] = js[r]; // rmin

    // Look for minimum S upto jmax within js
    int jh = (int) ( (i + istep <= imax) ? J[q][i + istep] : js[n-1] );

    int jmax = std::min((int)jh, (int)i);

    double sjimin(
        dissimilarity(criterion, jmax, i, sum_x, sum_x_sq, sum_w, sum_w_sq)
        // ssq(jmax, i, sum_x, sum_x_sq, sum_w)
      );

    for(++ r; r < n && js[r]<=jmax; r++) {

      const size_t & jabs = js[r];

      if(jabs > i) break;

      if(jabs < J[q-1][i]) continue;

      double s =
        dissimilarity(criterion, jabs, i, sum_x, sum_x_sq, sum_w, sum_w_sq);
        // (ssq(jabs, i, sum_x, sum_x_sq, sum_w));
      double Sj = (S[q-1][jabs-1] + s);

      if(Sj <= S[q][i]) {
        S[q][i] = Sj;
        J[q][i] = js[r];
      } else if(S[q-1][jabs-1] + sjimin > S[q][i]) {
        break;
      } /*else if(S[q-1][js[rmin]-1] + s > S[q][i]) {
 break;
      } */
    }
    r --;
    jl = jh;
  }
}

inline void find_min_from_candidates
  (int imin, int imax, int istep, int q,
   const vector<size_t> & js,
   vector< vector<double> > & S,
   vector< vector<size_t> > & J,
   const vector<double> & sum_x,
   const vector<double> & sum_x_sq,
   const vector<double> & sum_w,
   const vector<double> & sum_w_sq,
   const enum DISSIMILARITY criterion)
{
  size_t rmin_prev = (0);

  for(int i=(imin); i<=imax; i+=istep) {

    size_t rmin = (rmin_prev);

    // Initialization of S[q][i] and J[q][i]
    S[q][i] = S[q-1][js[rmin] - 1] +
      dissimilarity(criterion, js[rmin], i, sum_x, sum_x_sq, sum_w, sum_w_sq);
      // ssq(js[rmin], i, sum_x, sum_x_sq, sum_w);
    J[q][i] = js[rmin];

    for(size_t r = (rmin+1); r<js.size(); ++r) {

      const size_t & j_abs = (js[r]);

      if(j_abs < J[q-1][i]) continue;
      if(j_abs > i) break;

      double Sj = (S[q-1][j_abs - 1] +
        dissimilarity(criterion, j_abs, i, sum_x, sum_x_sq, sum_w, sum_w_sq));
        // ssq(j_abs, i, sum_x, sum_x_sq, sum_w));
      if(Sj <= S[q][i]) {
        S[q][i] = Sj;
        J[q][i] = js[r];
        rmin_prev = r;
      }
    }
  }
}

inline void SMAWK
  (int imin, int imax, int istep, int q,
   const vector<size_t> & js,
   vector< vector<double> > & S,
   vector< vector<size_t> > & J,
   const vector<double> & sum_x,
   const vector<double> & sum_x_sq,
   const vector<double> & sum_w,
   const vector<double> & sum_w_sq,
   const enum DISSIMILARITY criterion)
{
#ifdef DEBUG_REDUCE
  std::cout << "i:" << '[' << imin << ',' << imax << ']' << '+' << istep
            << std::endl;
#endif

  if(imax - imin <= 0 * istep) { // base case only one element left

    find_min_from_candidates(
      imin, imax, istep, q, js, S, J, sum_x, sum_x_sq, sum_w,
      sum_w_sq, criterion
    );

  } else {

    // REDUCE

#ifdef DEBUG //_REDUCE
    std::cout << "js:";
    for (size_t l=0; l < js.size(); ++l) {
      std::cout << js[l] << ",";
    }
    std::cout << std::endl;
    std::cout << std::endl;
#endif

    vector<size_t> js_odd;

    reduce_in_place(imin, imax, istep, q, js, js_odd,
                    S, J, sum_x, sum_x_sq, sum_w,
                    sum_w_sq, criterion);

    int istepx2 = (istep << 1);
    int imin_odd = (imin + istep);
    int imax_odd = (imin_odd + (imax - imin_odd) / istepx2 * istepx2);

    // Recursion on odd rows (0-based):
    SMAWK(imin_odd, imax_odd, istepx2,
          q, js_odd, S, J, sum_x, sum_x_sq, sum_w,
          sum_w_sq, criterion);

#ifdef DEBUG // _REDUCE
    std::cout << "js_odd (reduced):";
    for (size_t l=0; l<js_odd.size(); ++l) {
      std::cout << js_odd[l] << ",";
    }
    std::cout << std::endl << std::endl;

    std::cout << "even pos:";
    for (int i=imin; i<imax; i+=istepx2) {
      std::cout << i << ",";
    }
    std::cout << std::endl << std::endl;

#endif

    fill_even_positions(imin, imax, istep, q, js,
                        S, J, sum_x, sum_x_sq, sum_w,
                        sum_w_sq, criterion);
  }
}

inline void fill_row_q_SMAWK
  (int imin, int imax, int q,
   vector< vector<double> > & S,
   vector< vector<size_t> > & J,
   const vector<double> & sum_x,
   const vector<double> & sum_x_sq,
   const vector<double> & sum_w,
   const vector<double> & sum_w_sq,
   const enum DISSIMILARITY criterion)
{
  // Assumption: each cluster must have at least one point.

  vector<size_t> js(imax-q+1);
  int abs = (q);
  std::generate(js.begin(), js.end(), [&] { return abs++; } );

  SMAWK(imin, imax, 1, q, js, S, J, sum_x, sum_x_sq, sum_w, sum_w_sq, criterion);
}

inline void fill_dp_matrix
  (const vector<double> & x, // data
   const vector<double> & w, // weight
   vector< vector< double > > & S,
   vector< vector< size_t > > & J,
   //const enum METHOD method,
   const enum DISSIMILARITY criterion)
  /*
   x: One dimension vector to be clustered, must be sorted (in any order).
   S: K x N matrix. S[q][i] is the sum of squares of the distance from
   each x[i] to its cluster mean when there are exactly x[i] is the
   last point in cluster q
   J: K x N backtrack matrix

   NOTE: All vector indices in this program start at position 0
   */
{
  const int K = (int) S.size();
  const int N = (int) S[0].size();

  vector<double> sum_x(N), sum_x_sq(N);
  vector<double> sum_w(w.size()), sum_w_sq(w.size());

  vector<int> jseq;

  double shift = x[N/2]; // median. used to shift the values of x to
  //  improve numerical stability

  if(w.empty()) { // equally weighted
    sum_x[0] = x[0] - shift;
    sum_x_sq[0] = (x[0] - shift) * (x[0] - shift);
  } else { // unequally weighted
    sum_x[0] = w[0] * (x[0] - shift);
    sum_x_sq[0] = w[0] * (x[0] - shift) * (x[0] - shift);
    sum_w[0] = w[0];
    sum_w_sq[0] = w[0] * w[0];
  }

  S[0][0] = 0;
  J[0][0] = 0;

  for(int i = 1; i < N; ++i) {

    if(w.empty()) { // equally weighted
      sum_x[i] = sum_x[i-1] + x[i] - shift;
      sum_x_sq[i] = sum_x_sq[i-1] + (x[i] - shift) * (x[i] - shift);
    } else { // unequally weighted
      sum_x[i] = sum_x[i-1] + w[i] * (x[i] - shift);
      sum_x_sq[i] = sum_x_sq[i-1] + w[i] * (x[i] - shift) * (x[i] - shift);
      sum_w[i] = sum_w[i-1] + w[i];
      sum_w_sq[i] = sum_w_sq[i-1] + w[i]*w[i];
    }

    // Initialize for q = 0
    S[0][i] = dissimilarity(criterion, 0, i, sum_x, sum_x_sq, sum_w, sum_w_sq); // ssq(0, i, sum_x, sum_x_sq, sum_w);
    J[0][i] = 0;
  }

#ifdef DEBUG
  for(size_t i=0; i<x.size(); ++i) {
    std::cout << x[i] << ',';
  }
  std::cout << std::endl;
  vector<vector<double>> SS(S);
  vector<vector<size_t>> JJ(J);

#endif

  for(int q = 1; q < K; ++q) {
    int imin;
    if(q < K - 1) {
      imin = std::max(1, q);
    } else {
      // No need to compute S[K-1][0] ... S[K-1][N-2]
      imin = N-1;
    }

#ifdef DEBUG
    // std::cout << std::endl << "q=" << q << ":" << std::endl;
#endif
    // fill_row_k_linear_recursive(imin, N-1, 1, q, jseq, S, J, sum_x, sum_x_sq);
    // fill_row_k_linear(imin, N-1, q, S, J, sum_x, sum_x_sq);
    //if(method == LINEAR) {
      fill_row_q_SMAWK(imin, N-1, q, S, J, sum_x, sum_x_sq, sum_w, sum_w_sq, criterion);
    //} else if(method == LOGLINEAR) {
      //fill_row_q_log_linear(imin, N-1, q, q, N-1, S, J, sum_x, sum_x_sq, sum_w, sum_w_sq, criterion);
    //} else { // QUADRATIC
      //fill_row_q(imin, N-1, q, S, J, sum_x, sum_x_sq, sum_w, sum_w_sq, criterion);
    //}

#ifdef DEBUG

    fill_row_q_log_linear(imin, N-1, q, q, N-1, SS, JJ, sum_x, sum_x_sq, sum_w, sum_w_sq, criterion);

    for(int i=imin; i<N; ++i) {
      if(S[q][i] != SS[q][i] || J[q][i] != JJ[q][i]) {
        std::cout << "ERROR: q=" << q << ", i=" << i << std::endl;
        std::cout << "\tS=" << S[q][i] << "\tJ=" << J[q][i] << std::endl;
        std::cout << "Truth\tSS=" << SS[q][i] << "\tJJ=" << JJ[q][i];
        std::cout << std::endl;
        assert(false);

      } else {
        /*
         std::cout << "OK: q=" << q << ", i=" << i << std::endl;
         std::cout << "\tS=" << S[q][i] << "\tJ=" << J[q][i] << std::endl;
         std::cout << "Truth\tSS=" << SS[q][i] << "\tJJ=" << JJ[q][i];
         std::cout << std::endl;
         */
      }

    }
#endif
  }

#ifdef DEBUG
  std::cout << "Linear & log-linear code returned identical dp index matrix."
            << std::endl;
#endif
}

inline void backtrack_L1
  (const vector<double> & x,
   const vector< vector< size_t > > & J,
   size_t* cluster, double* centers, double* withinss,
   size_t* count)
{
  const size_t K = J.size();
  const size_t N = J[0].size();
  size_t cluster_right = N-1;
  size_t cluster_left;

  // Backtrack the clusters from the dynamic programming matrix
  for(int q = ((int)K)-1; q >= 0; --q) {
    cluster_left = J[q][cluster_right];

    for(size_t i = cluster_left; i <= cluster_right; ++i)
      cluster[i] = q;

    centers[q] = x[(cluster_right+cluster_left) >> 1];

    for(size_t i = cluster_left; i <= cluster_right; ++i)
      withinss[q] += std::fabs(x[i] - centers[q]);

    count[q] = cluster_right - cluster_left + 1;

    if(q > 0) {
      cluster_right = cluster_left - 1;
    }
  }
}

inline void backtrack_L2
  (const vector<double> & x,
   const vector< vector< size_t > > & J,
   size_t* cluster, double* centers, double* withinss,
   size_t* count)
{
  const size_t K = J.size();
  const size_t N = J[0].size();
  size_t cluster_right = N-1;
  size_t cluster_left;

  // Backtrack the clusters from the dynamic programming matrix
  for(int q = ((int)K)-1; q >= 0; --q) {
    cluster_left = J[q][cluster_right];

    for(size_t i = cluster_left; i <= cluster_right; ++i)
      cluster[i] = q;

    double sum = 0.0;

    for(size_t i = cluster_left; i <= cluster_right; ++i)
      sum += x[i];

    centers[q] = sum / (cluster_right-cluster_left+1);

    for(size_t i = cluster_left; i <= cluster_right; ++i)
      withinss[q] += (x[i] - centers[q]) * (x[i] - centers[q]);

    count[q] = cluster_right - cluster_left + 1;

    if(q > 0) {
      cluster_right = cluster_left - 1;
    }
  }
}

pair<vector<size_t>, vector<double>> kmeans_1d_dp(const vector<double> &x, size_t k,
                                                  const enum DISSIMILARITY criterion) {

  vector<size_t> cluster(x.size());
  vector<double> center(k);
  vector<double> withinss(k);
  vector<size_t> size(k);
  vector<double> y;
  vector<vector< double>> S(k, vector<double>(x.size()));
  vector<vector< size_t>> J(k, vector<size_t>(x.size()));

  if (criterion == L1) {
    fill_dp_matrix(x, y, S, J, L1);
    backtrack_L1(x, J, &cluster[0], &center[0], &withinss[0], &size[0]);
  } else {
    //EWL2::fill_dp_matrix(x, y, S, J);
    //backtrack_L2(x, J, &cluster[0], &centers[0], &withinss[0], &size[0]);
  }

  return make_pair(cluster, center);
}
