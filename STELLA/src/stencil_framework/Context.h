#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/at.hpp>
#include "SharedInfrastructure.h"
#include "FunctionParameter.h"
#include "FunctionParameterEvaluator.h"
#include "DataParameter.h"
#include "DataParameterEvaluator.h"
#include "StencilStage.h"
#include "HasParameter.h"

#include <boost/mpl/map/map40.hpp> // depends on MAX_TUPLE_SIZE
#include <boost/mpl/set/set40.hpp> // depends on MAX_TUPLE_SIZE

/**
* @struct ContextDescriptor
* Structure containing the type information needed to setup a context object
*/
template<
    typename TStencilConfiguration,
    typename TTupleElements>
struct ContextDescriptor
{
    // define stencil configuration information
    typedef TStencilConfiguration StencilConfiguration;
    typedef typename TStencilConfiguration::ValueType ValueType;
    typedef typename TStencilConfiguration::BlockSize BlockSize;
    
    // define tuple elements information
    typedef TTupleElements TupleElements;
    typedef typename TupleElements::ElementIndexes ElementIndexes;
    typedef typename TupleElements::ElementTypes ElementTypes;
    
    // compute a map which maps element indexes to element types
    // used to compute the parameter return type
    typedef typename boost::mpl::fold<
        boost::mpl::range_c<int, 0, boost::mpl::size<ElementIndexes>::value>,
        boost::mpl::map0<>,
        boost::mpl::insert<
            boost::mpl::_1,
            boost::mpl::pair<
                boost::mpl::at<ElementIndexes, boost::mpl::_2>,
                boost::mpl::at<ElementTypes, boost::mpl::_2>
            >
        >
    >::type IndexToElementMap;
};

/**
* @class Context
* The context class is a wrapper around the parameter tuple of a stencil. 
* It offers a bracket operator which evaluates data and function parameters. 
* At the same time the class offers methods to move the iterators / pointers in the parameter tuple. 
* There are specialized implementations for the different back ends of the library (CUDA, OpenMP)
*/
template<
    typename TContextDescriptor,
    typename TContextImpl>
class Context
{
    DISALLOW_COPY_AND_ASSIGN(Context);
protected: 
    __ACC__
    Context() {}
    __ACC__
    ~Context() {}

public:
    typedef typename TContextDescriptor::ValueType ValueType;
    typedef typename TContextDescriptor::BlockSize BlockSize;
    // compute a map which maps element indexes to element types
    // used to compute the parameter return type
    typedef typename TContextDescriptor::IndexToElementMap IndexToElementMap;

    /**
    * Advance the context position by a constant offset
    * The iteration mask allows to apply the method to a subset of the context
    */
    template<typename TIterationMask, int VI, int VJ, int VK>
    __ACC__
    void Advance()
    {
        static_cast<TContextImpl*>(this)->template AdvanceImpl<TIterationMask, VI, VJ, VK>();
    }

    /**
    * Advance the context position by a variable a k offset
    * The iteration mask allows to apply the method to a subset of the context
    * @param kOffset kOffset relative to the current position
    */
    template<typename TIterationMask, int VI, int VJ>
    __ACC__
    void Advance(const int kOffset)
    {
        static_cast<TContextImpl*>(this)->template AdvanceImpl<TIterationMask, VI, VJ>(kOffset);
    }

    /**
    * Restore the context position and move the position by an ijk offset
    * The iteration mask allows to apply the method to a subset of the context
    * @param iOffset iOffset relative to block origin
    * @param jOffset jOffset relative to block origin
    */
    template<typename TIterationMask, int VK>
    __ACC__
    void RestoreAndAdvance(const int iOffset, const int jOffset)
    {
        static_cast<TContextImpl*>(this)->template RestoreAndAdvanceImpl<TIterationMask, VK>(iOffset, jOffset);
    }

    /**
    * Restore the context position and move the position by an ijk offset
    * The iteration mask allows to apply the method to a subset of the context
    * @param iOffset iOffset relative to block origin
    * @param jOffset jOffset relative to block origin
    * @param kOffset kOffset relative to block origin
    */
    template<typename TIterationMask>
    __ACC__
    void RestoreAndAdvance(const int iOffset, const int jOffset, const int kOffset)
    {
        static_cast<TContextImpl*>(this)->template RestoreAndAdvanceImpl<TIterationMask>(iOffset, jOffset, kOffset);
    }

    /**
    * Method sliding the context caches
    */
    template<typename TSlideMask>
    __ACC__
    void SlideCaches() 
    {
        static_cast<TContextImpl*>(this)->template SlideCachesImpl<TSlideMask>();
    }
    
    /**
    * Method filling the cache front elements
    */
    template<
        typename TFillMask,
        typename TDomain>
    __ACC__
    void FillCachesFront() 
    {
        static_cast<TContextImpl*>(this)->template FillCachesFrontImpl<TFillMask, TDomain>();
    }

	/**
    * Method filling all cache elements except for the back elements
    */
    template<
        typename TFillMask,
        typename TDomain>
    __ACC__
    void FillCaches() 
    {
        static_cast<TContextImpl*>(this)->template FillCachesImpl<TFillMask, TDomain>();
    }

	/**
    * Method flushing the cache back elements
    */
    template<
        typename TFlushMask,
        typename TDomain>
    __ACC__
    void FlushCachesBack() 
    {
        static_cast<TContextImpl*>(this)->template FlushCachesBackImpl<TFlushMask, TDomain>();
    }

	/**
    * Method flushing all cache elements except for the back elements
    */
    template<
        typename TFlushMask,
        typename TDomain>
    __ACC__
    void FlushCaches() 
    {
        static_cast<TContextImpl*>(this)->template FlushCachesImpl<TFlushMask, TDomain>();
    }

    /**
    * Overload of the bracket operator accessing a parameter tuple element
    */
    template<
        typename TOffset,
        int VParameterIndex>
    __ACC__
    typename DataParameterEvaluator<
        IndexToElementMap,
        VParameterIndex
    >::ReturnType
    operator[] (DataParameter<VParameterIndex, TOffset>)
    {
        return static_cast<TContextImpl*>(this)->template AccessImpl<
            typename DataParameterEvaluator<
                IndexToElementMap,
                VParameterIndex
            >::ReturnType,
            VParameterIndex,
            TOffset,
            cUseCache
        >();
    }

    /**
    * Overload of the bracket operator which accesses a tuple element
    * where the offset is defined at run time
    */
    template<int VParameterIndex>
    __ACC__
    typename DataParameterEvaluator<
        IndexToElementMap,
        VParameterIndex
    >::ReturnType
    operator[] (const DataParameterDynamic<VParameterIndex> offset)
    {
        return static_cast<TContextImpl*>(this)->template AccessImpl<
            typename DataParameterEvaluator<
                IndexToElementMap,
                VParameterIndex
            >::ReturnType,
            VParameterIndex,
            cUseCache
        >(offset);
    }

    /**
    * Access method bypassing the software cache 
    */
    template<
        typename TOffset,
        int VParameterIndex>
    __ACC__
    typename DataParameterEvaluator<
        IndexToElementMap,
        VParameterIndex
    >::ReturnType
    BypassCache(DataParameter<VParameterIndex, TOffset>)
    {
        return static_cast<TContextImpl*>(this)->template AccessImpl<
            typename DataParameterEvaluator<
                IndexToElementMap,
                VParameterIndex
            >::ReturnType,
            VParameterIndex,
            TOffset,
            cBypassCache
        >();
    }

    /**
    * Overload of the bracket operator which executes a stencil function
    */
    template<
        template<typename> class TStencilFunction,
        typename TDomain,
        typename TOffset,
        typename TIterationMask>
    __ACC__
    typename FunctionParameterEvaluator<
        Context,
        TStencilFunction,
        TDomain,
        TOffset,
        TIterationMask
    >::ReturnType
    operator[] (FunctionParameter<TStencilFunction, TDomain, TOffset, TIterationMask>)
    {
        return FunctionParameterEvaluator<
            Context,
            TStencilFunction,
            TDomain,
            TOffset,
            TIterationMask
        >::Evaluate(*this);
    }

    /**
    * @return the calculation domain size in k dimension
    */
    __ACC__
    int kSize() const { return static_cast<const TContextImpl*>(this)->kSizeImpl(); } 
};



  
