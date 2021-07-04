#ifndef _SERIALIZE_BOOST_ROWVECTOR_H_
#define _SERIALIZE_BOOST_ROWVECTOR_H_

#include <boost/serialization/serialization.hpp>
#include <newmat.h>

namespace boost {
namespace serialization {

template <class Archive>
inline void serialize(Archive & ar, RowVector & t, const unsigned int file_version)
{
   int dim = t.Ncols();
   ar & dim;
   for (int i=0; i<dim; i++){
       ar & t.element(i);
   }
}

} 
} 

#endif 
