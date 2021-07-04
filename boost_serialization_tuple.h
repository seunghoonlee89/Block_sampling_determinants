#ifndef _SERIALIZE_BOOST_TUPLE_H_
#define _SERIALIZE_BOOST_TUPLE_H_

#include <boost/serialization/serialization.hpp>
#include <tuple>

namespace boost {
namespace serialization {

template <class Archive, typename T0, typename T1, typename T2, typename T3>
inline void serialize(Archive & ar, std::tuple<T0, T1, T2, T3> & t, const unsigned int file_version)
{
   ar & std::get<0>(t) & std::get<1>(t) & std::get<2>(t) & std::get<3>(t);
}

} 
} 

#endif 
