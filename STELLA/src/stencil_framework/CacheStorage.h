#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/copy_if.hpp>
#include <boost/mpl/push_front.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/plus.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/eval_if.hpp>
#include "SharedInfrastructure.h"
#include "Offset.h"
#include "KRange.h"
#include "KWindow.h"
#include "Cache.h"
#include "CacheStorageFunctors.h"

#include <boost/mpl/vector/vector10.hpp>

/**
* @class KCacheStorage
* Class implementing a KCache as sliding window with a given size stored using an array
*/
template<
    typename TValue,
    typename TKWindow,
    typename TKRange, 
    int VParameterIndex, 
    CacheIOPolicy VIOPolicy,
    KLoopDirection VKLoopDirection>
class KCacheStorage
{
    DISALLOW_COPY_AND_ASSIGN(KCacheStorage);
public:
    BOOST_STATIC_ASSERT(is_k_window<TKWindow>::value);
    BOOST_STATIC_ASSERT(is_k_range<TKRange>::value);

    // define the array size
    typedef boost::mpl::integral_c<int, TKWindow::KPlusOffset::value - TKWindow::KMinusOffset::value + 1> Size;

    // define the cache element slide, fill and flush functors
    typedef KCacheStorageSlideFunctor SlideFunctor;
    typedef KCacheStorageFillFunctor<
        DataParameter<VParameterIndex, Offset<0, 0, 0> >, 
        TKWindow, 
        VKLoopDirection
    > FillFunctor;
    typedef KCacheStorageFlushFunctor<
        DataParameter<VParameterIndex, Offset<0, 0, 0> >, 
        TKWindow, 
        VKLoopDirection
    > FlushFunctor;

    // compute a vector holding all cache element k offsets
    typedef typename boost::mpl::eval_if_c<
        VKLoopDirection == cKIncrement,
        boost::mpl::fold<
            boost::mpl::range_c<int, TKWindow::KMinusOffset::value, TKWindow::KPlusOffset::value + 1>,
            boost::mpl::vector0<>,
            boost::mpl::push_back<boost::mpl::_1, boost::mpl::_2>
        >,
        boost::mpl::fold<
            boost::mpl::range_c<int, TKWindow::KMinusOffset::value, TKWindow::KPlusOffset::value + 1>,
            boost::mpl::vector0<>,
            boost::mpl::push_front<boost::mpl::_1, boost::mpl::_2>
        >
    >::type AllKOffsets;
    
    // define array type
    typedef TValue& ReturnType;
    typedef acc::Array<TValue, Size::value> ArrayType;
    
    __ACC__
    KCacheStorage() {}
    __ACC__
    ~KCacheStorage() {}

    /**
    * Method sliding the cache by one position
    */
    __ACC__
    void Slide() 
    {
        // compute the slide array indexes
        typedef typename ArrayType::ElementIndexes ArrayIndexes;
        typedef typename boost::mpl::reverse_copy<
            typename boost::mpl::eval_if<
                boost::mpl::empty<ArrayIndexes>,
                boost::mpl::identity<ArrayIndexes>,
                boost::mpl::pop_back<ArrayIndexes>
            >::type 
        >::type SlideArrayIndexes;

        // execute the slide
        acc::apply_to_all<
            SlideFunctor,
            SlideArrayIndexes,
            ArrayType
        >(cacheArray_);
    }

    /**
    * Method filling the first cache entry from main memory
    * @param context context object used to fill the cache
    */
    template<
        typename TDomain,
        typename TContext>
    __ACC__
    void FillFront(TContext& context)
    {
        // define the fill front offset
        typedef typename boost::mpl::if_c<
            VKLoopDirection == cKIncrement,
            typename TKWindow::KPlusOffset,
            typename TKWindow::KMinusOffset
        >::type FrontKOffset;

        // check if a fill is necessary and return a vector containing the front offset
        // (note if no fill is necessary an empty vector is returned)
        typedef typename boost::mpl::if_<
            in_k_range<TKRange, TDomain, FrontKOffset>,
            boost::mpl::vector1<FrontKOffset>,
            boost::mpl::vector0<>
        >::type FillFrontKOffsets;

        // execute the fill
        acc::apply_to_all<
            FillFunctor,
            FillFrontKOffsets,
            ArrayType,
            TContext
        >(cacheArray_, context);
    }

    /**
    * Method filling all cache elements from main memory
    * (except for the back of the cache which is slided out before the first stage execution)
    * @param context context object used to fill the cache
    */
    template<
        typename TDomain,
        typename TContext>
    __ACC__
    void Fill(TContext& context)
    {
        // filter out all k offsets which are outside the cache k range
        typedef typename boost::mpl::copy_if<
            AllKOffsets, 
            in_k_range<TKRange, TDomain, boost::mpl::_> 
        >::type FillKOffsets;

        // execute the fill
        acc::apply_to_all<
            FillFunctor,
            FillKOffsets,
            ArrayType,
            TContext
        >(cacheArray_, context);
    }

    /**
    * Method flushing the last cache entry back to main memory
    * @param context context object used to flush the cache
    */
    template<
        typename TDomain,
        typename TContext>
    __ACC__
    void FlushBack(TContext& context) 
    {
        // define the flush back offset
        typedef typename boost::mpl::if_c<
            VKLoopDirection == cKIncrement,
            typename TKWindow::KMinusOffset,
            typename TKWindow::KPlusOffset
        >::type BackKOffset;

        // check if a flush is necessary and return a vector containing the back offset
        // (note if no flush is necessary an empty vector is returned)
        typedef typename boost::mpl::if_<
            in_k_range<TKRange, TDomain, BackKOffset>,
            boost::mpl::vector1<BackKOffset>,
            boost::mpl::vector0<>
        >::type FlushBackKOffsets;
            
        // execute the flush
        acc::apply_to_all<
            FlushFunctor,
            FlushBackKOffsets,
            ArrayType,
            TContext
        >(cacheArray_, context);
    }

    /**
    * Method flushing all cache elements to main memory
    * (except for the back of the cache which was flushed after the last stage execution)
    * @param context context object used to flush the cache
    */
    template<
        typename TDomain,
        typename TContext>
    __ACC__
    void Flush(TContext& context)
    {
        // filter out all k offsets which are outside the cache k range
        typedef typename boost::mpl::copy_if<
            AllKOffsets, 
            in_k_range<TKRange, TDomain, boost::mpl::_> 
        >::type FlushKOffsets;

        // execute the flush
        acc::apply_to_all<
            FlushFunctor,
            FlushKOffsets,
            ArrayType,
            TContext
        >(cacheArray_, context);
    }

    /**
    * Method accessing a cache element at a given offset
    * @return cached value
    */
    template<typename TOffset>
    __ACC__
    ReturnType At(TOffset*) 
    { 
        // check the I and J offsets are zero
        BOOST_STATIC_ASSERT(TOffset::I::value == 0 && TOffset::J::value == 0);

        // define the K access offset
        typedef typename TOffset::K KOffset;

        // compute the cache array index for the given k offset
        typedef typename compute_k_cache_storage_array_index<
            KOffset,
            TKWindow,
            VKLoopDirection
        >::type ArrayIndex;

        // return the cached element
        return cacheArray_.At(static_cast<ArrayIndex*>(0)); 
    }

private:
    ArrayType cacheArray_;
};

/**
* @struct create_k_cache_strorage
* Meta function returning the k cache storage for a given k cache
*/
template<
    typename TCache,
    typename TParameter,
    typename TKLoopDirection>
struct create_k_cache_strorage; 

template<
    typename TParameter,
    typename TKLoopDirection,
    int VParameterIndex,
    CacheIOPolicy VIOPolicy,
    typename TKWindow,
    typename TKRange>
struct create_k_cache_strorage<KCache<VParameterIndex, VIOPolicy, TKWindow, TKRange>, TParameter, TKLoopDirection>
{
    typedef KCacheStorage<
        typename value_type<TParameter>::type,
        TKWindow, 
        TKRange,
        VParameterIndex,
        VIOPolicy,
        TKLoopDirection::value
    > type;
};

/**
* @class IJCacheStorage
* Class implementing an IJCache using a block storage storing a single ij level
*/
template<
    typename TValue,
    typename TIJRange,
    typename TBlockSize>
class IJCacheStorage
{
    DISALLOW_COPY_AND_ASSIGN(IJCacheStorage);
public:
    BOOST_STATIC_ASSERT(is_ij_range<TIJRange>::value);
    BOOST_STATIC_ASSERT(is_block_size<TBlockSize>::value);

    // define the return type
    typedef TValue& ReturnType;

    // define a block storage type used to store one ij level
    typedef BlockStorage<
        TValue, 
        BlockStorageFormat<
            TBlockSize, 
            DataFieldIJBoundary<
                TIJRange::IMinusOffset::value, 
                TIJRange::IPlusOffset::value,
                TIJRange::JMinusOffset::value, 
                TIJRange::JPlusOffset::value
            >
        >
    > BlockStorageType;

    __ACC__
    IJCacheStorage() {}
    __ACC__
    ~IJCacheStorage() {}

    /**
    * Method accessing a cache element at a given offset
    * @param iOffset i offset relative to the block origin
    * @param jOffset j offset relative to the block origin
    * @return cached value
    */
    __ACC__
    ReturnType At(const int iOffset, const int jOffset) 
    { 
        return blockStorage_.At(iOffset, jOffset);
    }

private:
    BlockStorageType blockStorage_;
};

/**
* @struct create_ij_cache_strorage
* Meta function returning the ij cache storage for a given ij cache
*/
template<
    typename TCache,
    typename TParameter,
    typename TIJRange,
    typename TBlockSize>
struct create_ij_cache_strorage; 

template<
    typename TParameter,
    typename TIJRange,
    typename TBlockSize,
    int VParameterIndex,
    typename TKRange>
struct create_ij_cache_strorage<IJCache<VParameterIndex, TKRange>, TParameter, TIJRange, TBlockSize>
{
    typedef IJCacheStorage<
        typename value_type<TParameter>::type,
        TIJRange,
        TBlockSize
    > type;
};


