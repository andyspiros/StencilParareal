#pragma once

#include <boost/mpl/integral_c.hpp>
#include "Definitions.h"

/**
* @struct DataParameterDynamic
* Type used to access data fields stored in the context 
*/
template<int VParameterIndex>
struct DataParameterDynamic
{
    // define types
    typedef boost::mpl::integral_c<int, VParameterIndex> ParameterIndex;
    const int iOffset;
    const int jOffset;
    const int kOffset;

    __ACC__
    DataParameterDynamic(int i, int j, int k)
        : iOffset(i), jOffset(j), kOffset(k)
    {};
};

