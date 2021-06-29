#ifndef MPS2CI_H 
#define MPS2CI_H 

#include <vector>
#include "global.h"
#include "timer.h"
#include <tuple>
#include <utility>
#include <unordered_map>
#include "CIcoeff_from_mps_TRIE.h"
#include <cmath>        // std::abs
#include "boost_serialization_unordered_map.h"
#include "boost_serialization_tuple.h"
//#include "boost_serialization_rowvector.h"

// SQ(d, leftq_for_d, rightq_for_d+1, Ms_for_d + dotMs_for_d+1) 
typedef std::tuple<int, int, int, int> SQ;

struct key_hash2 : public std::unary_function<SQ, std::size_t>
{
   std::size_t operator()(const SQ& k) const
   {
      return std::get<0>(k) ^ std::get<1>(k) ^ std::get<2>(k) ^ std::get<3>(k);
   }
};

struct key_equal2 : public std::binary_function<SQ, SQ, bool>
{
   bool operator()(const SQ& v0, const SQ& v1) const
   {
      return (
               std::get<0>(v0) == std::get<0>(v1) &&
               std::get<1>(v0) == std::get<1>(v1) &&
               std::get<2>(v0) == std::get<2>(v1) &&
               std::get<3>(v0) == std::get<3>(v1)
             );
   }
};

//  ex_typ:
//  rf, a, b, aa, ab, bb, aaa, aab, abb, bbb,
//  aaaa, aaab, aabb, abbb, bbbb, etc 
typedef std::string ex_typ;

//  occ_string:
//  2 : doubly occupied
//  a : alpha  occupied 
//  b : beta   occupied 
//  0 : unoccupied 
typedef std::string partial_occ;

// ex_idx:
// pair< pair<alpha_hole_idx, alpha_ptcl_idx>, pair<beta_hole_idx, beta_ptcl_idx> >
typedef std::pair< std::pair<std::vector<int>, std::vector<int>>, std::pair<std::vector<int>, std::vector<int>> > partial_ex;

typedef std::unordered_map<SQ,RowVector,key_hash2,key_equal2> partial_wv;

template <typename T>
vector<size_t> sort_indices(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return std::abs(v[i1]) > std::abs(v[i2]);});

  return idx;
}


#endif
