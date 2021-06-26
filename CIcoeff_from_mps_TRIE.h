#ifndef CITRIE_H 
#define CITRIE_H 

#include <vector>
#include "global.h"
#include "timer.h"
#include <tuple>
#include <utility>
#include <unordered_map>

// SQ(d, leftq_for_d, rightq_for_d+1, Ms_for_d + dotMs_for_d+1) 
typedef std::tuple<int, int, int, int> SQ;

struct key_hash : public std::unary_function<SQ, std::size_t>
{
   std::size_t operator()(const SQ& k) const
   {
      return std::get<0>(k) ^ std::get<1>(k) ^ std::get<2>(k) ^ std::get<3>(k);
   }
};

struct key_equal : public std::binary_function<SQ, SQ, bool>
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

typedef std::unordered_map<SQ,RowVector,key_hash,key_equal> SQ_RV;

#endif
