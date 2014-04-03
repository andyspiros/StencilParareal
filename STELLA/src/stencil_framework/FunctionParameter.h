#pragma once

#include "Offset.h"

/**
* @struct FunctionParameter
* Type used to run a function using the context object
* (Note that the Center and At methods allow us passing a FunctionParameter as parameter to another stencil function)
*/
template<
    template<typename> class TStencilFunction,
    typename TDomain,
    typename TOffset,
    typename TParameterList>
struct FunctionParameter
{
    // define types
    typedef TOffset Offset;
    typedef TDomain Domain;
    typedef TParameterList ParameterList;
    
    /**
    * @return a function parameter object with offset (0,0,0)
    */
    __ACC__
    static FunctionParameter<
        TStencilFunction, 
        TDomain,
        TOffset, 
        TParameterList> 
    Center()
    {
        return FunctionParameter<
            TStencilFunction, 
            TDomain,
            TOffset, 
            TParameterList
        >();
    }

    /**
    * @return a function parameter object with offset TOffsetToAdd
    */
    template<typename TOffsetToAdd>
    __ACC__
    static FunctionParameter<
        TStencilFunction, 
        TDomain,
        typename add_offsets<TOffset, TOffsetToAdd>::type, 
        TParameterList> 
    At(TOffsetToAdd)
    {
        return FunctionParameter<
            TStencilFunction, 
            TDomain,
            typename add_offsets<TOffset, TOffsetToAdd>::type, 
            TParameterList
        >();
    }
};



  
