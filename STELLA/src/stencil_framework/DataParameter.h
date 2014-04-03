#pragma once

#include <boost/mpl/integral_c.hpp>
#include "Definitions.h"
#include "Offset.h"
#include "DataParameterDynamic.h"

/**
* @struct DataParameter
* Type used to access data fields stored in the context 
*/
template<
    int VParameterIndex,
    typename TOffset>
struct DataParameter
{
    // define types
    typedef boost::mpl::integral_c<int, VParameterIndex> ParameterIndex;
    typedef TOffset Offset;

    /**
    * @return a data parameter object with offset (0,0,0)
    */
    __ACC__
    static DataParameter<VParameterIndex, TOffset> Center()
    {
        return DataParameter<VParameterIndex, TOffset>();
    }

    /**
    * @return a data parameter object with offset TOffsetToAdd
    */
    template<typename TOffsetToAdd>
    __ACC__
    static DataParameter<VParameterIndex, typename add_offsets<TOffset, TOffsetToAdd>::type> At(TOffsetToAdd)
    {
        return DataParameter<VParameterIndex, typename add_offsets<TOffset, TOffsetToAdd>::type>();
    }

    /**
    * @return a data parameter object with offset TOffsetToAdd
    */
    __ACC__
    static DataParameterDynamic<VParameterIndex> At(int i, int j, int k)
    {
        return DataParameterDynamic<VParameterIndex>(i, j, k);
    }
};

// define a dummy parameter
typedef DataParameter<-1, Offset<0, 0, 0> > DummyParameter;

