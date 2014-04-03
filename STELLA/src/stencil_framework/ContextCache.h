#pragma once

#include <boost/mpl/assert.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/copy_if.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/if.hpp>
#include "SharedInfrastructure.h"
#include "StencilSweepDescriptor.h"
#include "StencilSweepGroupDescriptor.h"
#include "Cache.h"
#include "ContextCacheFunctors.h"

/**
* @class ContextCache
* Class encapsulating software cache functionality which can be used by context implementations
*/
template<
    typename TContextDescriptor,
    typename TStencilSweepDescriptor> 
class ContextCache;

// specialization for empty cache
template<
    typename TContextDescriptor,
    typename TStencilStages,
    KLoopDirection VKLoopDirection> 
class ContextCache<TContextDescriptor, StencilSweepDescriptor<boost::mpl::void_, TStencilStages, VKLoopDirection> >
{
public:
    // define an empty register cache element map
    typedef boost::mpl::map0<> IndexToCacheElementMap;
    
    // define a void ij cache tuple
    typedef boost::mpl::void_ IJCacheTuple;

    __ACC__
    ContextCache() {};
    __ACC__
    ~ContextCache() {};

    /**
    * Method sliding all caches
    */
    template<typename TSlideMask>
    __ACC__
    void Slide() {}

	/**
    * Method filling the front element of all caches
    * @param context context used to fill the caches
    */
    template< 
        typename TFillMask,
        typename TDomain,
        typename TContext>
    __ACC__
    void FillFront(TContext& context) {}

	/**
    * Method filling all elements except for the front element of all caches
    * @param context context used to fill the caches
    */
    template< 
        typename TFillMask,
        typename TDomain,
        typename TContext>
    __ACC__
    void Fill(TContext& context) {}

	/**
    * Method flushing the back element of all caches
    * @param context context used to flush the caches
    */
    template< 
        typename TFlushMask,
        typename TDomain,
        typename TContext>
    __ACC__
    void FlushBack(TContext& context) {}

	/**
    * Method flushing all elements except for the back element of all caches
    * @param context context used to flush the caches
    */
    template< 
        typename TFlushMask,
        typename TDomain,
        typename TContext>
    __ACC__
    void Flush(TContext& context) {}
};

// specialization for populated cache
template<
    typename TContextDescriptor,
    typename TCaches,
    typename TStencilStages,
    KLoopDirection VKLoopDirection> 
class ContextCache<TContextDescriptor, StencilSweepDescriptor<TCaches, TStencilStages, VKLoopDirection> >
{
public:
    // extract the index to element map from the context descriptor
    typedef typename TContextDescriptor::IndexToElementMap IndexToElementMap;

    // remove caches which are not used by the stencil stages
    typedef typename boost::mpl::copy_if<
        TCaches,
        stencil_stages_have_parameter<TStencilStages, FullDomain, cache_index<boost::mpl::_> >
    >::type Caches;

    // make sure all caches are used by at least one stencil stage
    BOOST_MPL_ASSERT_MSG( 
        (boost::mpl::size<Caches>::value == boost::mpl::size<TCaches>::value),
        CONTEXT_CACHE_ELEMENT_NOT_USED_BY_STAGES,
        (TCaches) 
    );

    // separate k and ij caches
    typedef typename boost::mpl::copy_if<
        Caches,
        is_k_cache<boost::mpl::_>
    >::type KCaches;
    typedef typename boost::mpl::copy_if<
        Caches,
        is_ij_cache<boost::mpl::_>
    >::type IJCaches;
    
    // define the k cache elements
    typedef typename k_cache_elements<
        KCaches,
        IndexToElementMap,
        VKLoopDirection
    >::type KCacheElements;

    // define the ij cache elements
    typedef typename ij_cache_elements<
        IJCaches,
        IndexToElementMap,
        TStencilStages,
        typename TContextDescriptor::BlockSize
    >::type IJCacheElements;

    // define a map which associates a parameter indexes to a caches
    typedef typename create_index_to_cache_element_map<
        KCacheElements,
        IJCacheElements
    >::type IndexToCacheElementMap;

    // define an ij cache tuple type
    typedef acc::Tuple<IJCacheElements> IJCacheTuple;

    __ACC__
    ContextCache() {};
    __ACC__
    ~ContextCache() {};

	/**
    * Method sliding all caches
    */
    template<typename TSlideMask>
    __ACC__
    void Slide() 
	{
		acc::modify_tuple<
            ContextCacheSlideFunctor,
            TSlideMask
        >(kCacheTuple_);
	}

	/**
    * Method filling the front element of all caches
    * @param context context used to fill the caches
    */
    template< 
        typename TFillMask,
        typename TDomain,
        typename TContext>
    __ACC__
    void FillFront(TContext& context) 
	{
		acc::modify_tuple<
            ContextCacheFillFrontFunctor<TContext, TDomain>,
            TFillMask,
            TContext
        >(kCacheTuple_, context);
	}

	/**
    * Method filling all elements except for the back elements of all caches
    * @param context context used to fill the caches
    */
    template< 
        typename TFillMask,
        typename TDomain,
        typename TContext>
    __ACC__
    void Fill(TContext& context) 
	{
		acc::modify_tuple<
            ContextCacheFillFunctor<TContext, TDomain>,
            TFillMask,
            TContext
        >(kCacheTuple_, context);
	}

	/**
    * Method flushing the back element of all caches
    * @param context context used to flush the caches
    */
    template< 
        typename TFlushMask,
        typename TDomain,
        typename TContext>
    __ACC__
    void FlushBack(TContext& context) 
	{
		acc::modify_tuple<
            ContextCacheFlushBackFunctor<TContext, TDomain>,
            TFlushMask,
            TContext
        >(kCacheTuple_, context);
	}

	/**
    * Method flushing all elements except for the back elements of all caches
    * @param context context used to flush the caches
    */
    template< 
        typename TFlushMask,
        typename TDomain,
        typename TContext>
    __ACC__
    void Flush(TContext& context) 
	{
		acc::modify_tuple<
            ContextCacheFlushFunctor<TContext, TDomain>,
            TFlushMask,
            TContext
        >(kCacheTuple_, context);
	}

    /**
    * @return the k cache tuple
    */
    __ACC__
    acc::Tuple<KCacheElements>& kCacheTuple() { return kCacheTuple_; }

private:
    acc::Tuple<KCacheElements> kCacheTuple_;
};

/**
* @struct compute_ij_cache_tuple
* Meta function computing the context ij cache tuple type given a context descriptor and a stencil sweep descriptor
*/
template<
    typename TContextDescriptor,
    typename TStencilSweepDescriptor> 
struct compute_ij_cache_tuple
{
    typedef typename ContextCache<
        TContextDescriptor, 
        TStencilSweepDescriptor
    >::IJCacheTuple type;
};

/**
* @struct maximum_ij_cache_tuple
* Meta function returning the maximum ij cache tuple given two tuples
*/
template<
    typename TIJCacheTuple1,
    typename TIJCacheTuple2> 
struct maximum_ij_cache_tuple
{
    typedef typename boost::mpl::if_c<
        (sizeof(TIJCacheTuple1) > sizeof(TIJCacheTuple2)),
        TIJCacheTuple1,
        TIJCacheTuple2
    >::type type;
};

/**
* @struct stencil_sweep_descriptor_maximum_ij_cache_tuple
* Meta function computing the maximum ij cache tuple given a vector of stencil sweeps
* (in case there are no ij caches void is returned as tuple type)
*/
template<
    typename TContextDescriptor,
    typename TStencilSweepDescriptors> 
struct stencil_sweep_descriptor_maximum_ij_cache_tuple
{
    typedef typename boost::mpl::fold<
        TStencilSweepDescriptors,
        boost::mpl::void_,
        maximum_ij_cache_tuple<
            boost::mpl::_1,
            compute_ij_cache_tuple<
                TContextDescriptor,
                boost::mpl::_2
            >
        >
    >::type type;
};

/**
* @struct stencil_sweep_group_descriptor_maximum_ij_cache_tuple
* Meta function computing the maximum ij cache tuple given a vector of stencil sweep groups
* (in case there are no ij caches void is returned as tuple type)
*/
template<
    typename TContextDescriptor,
    typename TStencilSweepGroupDescriptors> 
struct stencil_sweep_group_descriptor_maximum_ij_cache_tuple
{
    typedef typename boost::mpl::fold<
        TStencilSweepGroupDescriptors,
        boost::mpl::void_,
        maximum_ij_cache_tuple<
            boost::mpl::_1,
            stencil_sweep_descriptor_maximum_ij_cache_tuple<
                TContextDescriptor,
                stencil_sweep_group_descriptor_sweeps<boost::mpl::_2>
            >
        >
    >::type type;
};

