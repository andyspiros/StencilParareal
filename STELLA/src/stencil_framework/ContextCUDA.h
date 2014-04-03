#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/has_key.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/eval_if.hpp>
#include "SharedInfrastructure.h"
#include "LaunchConfiguration.h" 
#include "Context.h"
#include "ContextCUDAFunctors.h"
#include "ContextCache.h"

/**
* @struct ContextCUDASharedData
* Structure storing the relevant data of the context (handed over by the host)
*/
template<
    typename TContextDescriptor,
    template<typename> class TTuple,
    typename TBlockConfiguration, // set to float for host device copy, set to block configuration when used on device
    int VIJCacheStorageSize> // size in bytes used for the ij caching
struct ContextCUDASharedData
{
    // define the descriptor type 
    typedef TContextDescriptor Descriptor;

    // define the strides and data tuple elements
    typedef typename storage_strides_elements<TContextDescriptor>::type StorageStridesElements;
    typedef typename shared_data_elements<TContextDescriptor>::type DataElements;

    TTuple<StorageStridesElements> storageStridesTuple; // tuple with stride information
    TTuple<DataElements> dataTuple; // tuple holding the pointer data handed over by the host
    int kSize; // calculation domain size in k direction
    TBlockConfiguration blockConfiguration; // block size configuration
    double ijCacheStorage[VIJCacheStorageSize / sizeof(double) + 1]; // double aligned storage space for ij caching
};

// dummy context used for the block configuration offset computation
template<
    typename TContextDescriptor,
    template<typename> class TTuple> 
struct ContextCUDASharedData<TContextDescriptor, TTuple, boost::mpl::void_, 0>
{
    typedef TContextDescriptor Descriptor;
    typedef typename storage_strides_elements<TContextDescriptor>::type StorageStridesElements;
    typedef typename shared_data_elements<TContextDescriptor>::type DataElements;

    TTuple<StorageStridesElements> storageStridesTuple; 
    TTuple<DataElements> dataTuple; 
    int kSize; 
    float dummy; // dummy float variable at start position of block configuration
};

/**
* @class ContextCUDA
* Class implementing the context functionality for the CUDA back end.
*/
template<
    typename TSharedData,
    typename TStencilSweepDescriptor>
class ContextCUDA : public Context<typename TSharedData::Descriptor, ContextCUDA<TSharedData, TStencilSweepDescriptor> > // CRTP
{
    DISALLOW_COPY_AND_ASSIGN(ContextCUDA);
public:
    // define the context descriptor
    typedef typename TSharedData::Descriptor ContextDescriptor;

    // use the kSize method
    using Context<ContextDescriptor, ContextCUDA>::kSize;

    // use context type definitions
    typedef typename ContextDescriptor::ElementTypes ElementTypes;
    typedef typename ContextDescriptor::ElementIndexes ElementIndexes;
    typedef typename ContextDescriptor::IndexToElementMap IndexToElementMap;

    // define the context cache type
    typedef ContextCache<ContextDescriptor, TStencilSweepDescriptor> Cache;
    typedef typename Cache::IJCacheTuple IJCacheTuple;
        
    // function computing the iteration mask holding all parameter indexes 
    // used by the stencil stages of the sweep
    template<typename TStencilStages, typename TDomain>
    struct compute_iteration_mask
    {
        // first convert the mask indexes to stride types
        // and then add them to a set of all strides covered by the mask
        typedef typename boost::mpl::fold<
            ElementIndexes, 
            boost::mpl::set0<>,
            boost::mpl::if_<
                boost::mpl::and_<
                    stencil_stages_have_parameter<TStencilStages, TDomain, boost::mpl::_2>,
                    is_iterable<boost::mpl::at<IndexToElementMap, boost::mpl::_2> >
                >,
                boost::mpl::insert<
                    boost::mpl::_1, 
                    create_storage_strides<boost::mpl::at<IndexToElementMap, boost::mpl::_2> > 
                >, 
                boost::mpl::_1
            >
        >::type type;
    };
    
    // define the index tuple element type
    typedef typename TSharedData::StorageStridesElements::ElementIndexes StorageStridesIndexes;
    typedef TupleElements<
        StorageStridesIndexes,
        typename boost::mpl::transform<
            StorageStridesIndexes,
            DataFieldCUDAStorageIndex<boost::mpl::_>
        >::type
    > IndexElements;
    
    // define the register tuple element type
    typedef typename register_data_elements<
        ContextDescriptor,
        TStencilSweepDescriptor
    >::type RegisterDataElements;
   
    // define a default constructor
    __ACC__
    ContextCUDA() 
    { 
        pSharedData_ = NULL; 
        iBoundaryOffset_ = 0;
        jBoundaryOffset_ = 0;
        iThreadOffset_ = 0;
        jThreadOffset_ = 0;
    }
    __ACC__
    ~ContextCUDA() {}

    /**
    * Method initializing the context
    * @param pSharedData pointer to constant context data as storage strides and parameter pointers
    */
    __ACC__
    void Init(TSharedData* pSharedData) 
    { 
        pSharedData_ = pSharedData; 
    }

    /**
    * Advance the context position by a constant offset
    * The iteration mask allows to apply the method to a subset of the context
    */
    template<typename TIterationMask, int VI, int VJ, int VK>
    __ACC__
    void AdvanceImpl()
    {
        ACC_ASSERT(pSharedData_);
        acc::modify_2_tuples<
            AdvanceFunctor<VI, VJ, VK>,
            TIterationMask
        >(pSharedData_->storageStridesTuple, indexTuple_);
    }

    /**
    * Advance the context position by a variable a k offset
    * The iteration mask allows to apply the method to a subset of the context
    * @param kOffset kOffset relative to the current position
    */
    template<typename TIterationMask, int VI, int VJ>
    __ACC__
    void AdvanceImpl(const int kOffset)
    {
        ACC_ASSERT(pSharedData_);
        acc::modify_2_tuples<
            AdvanceInKFunctor<VI,VJ>,
            TIterationMask,
            const int
        >(pSharedData_->storageStridesTuple, indexTuple_, kOffset);
    }

    /**
    * Restore the context position and move the position by an ijk offset
    * The iteration mask allows to apply the method to a subset of the context
    * @param iOffset iOffset relative to block origin
    * @param jOffset jOffset relative to block origin
    */
    template<typename TIterationMask, int VK>
    __ACC__
    void RestoreAndAdvanceImpl(const int iOffset, const int jOffset)
    {
        ACC_ASSERT(pSharedData_);
        IJOffset ijOffset;
        ijOffset.iBlockIndex = pSharedData_->blockConfiguration.iBlockIndex;
        ijOffset.jBlockIndex = pSharedData_->blockConfiguration.jBlockIndex;
        ijOffset.i = iOffset;
        ijOffset.j = jOffset;

        acc::modify_2_tuples<
            RestoreAndAdvanceInIJFunctor<VK>,
            TIterationMask,
            const IJOffset
        >(pSharedData_->storageStridesTuple, indexTuple_, ijOffset);

        // set the boundary offsets
        iBoundaryOffset_ = computeBoundaryOffset(iOffset, pSharedData_->blockConfiguration.iStart, pSharedData_->blockConfiguration.iEnd);
        jBoundaryOffset_ = computeBoundaryOffset(jOffset, pSharedData_->blockConfiguration.jStart, pSharedData_->blockConfiguration.jEnd);

        // set the thread offsets
        iThreadOffset_ = iOffset;
        jThreadOffset_ = jOffset;
    }

    /**
    * Restore the context position and move the position by an ijk offset
    * The iteration mask allows to apply the method to a subset of the context
    * @param iOffset iOffset relative to block origin
    * @param jOffset jOffset relative to block origin
    */
    template<typename TIterationMask>
    __ACC__
    void RestoreAndAdvanceImpl(const int iOffset, const int jOffset, const int kOffset)
    {
        ACC_ASSERT(pSharedData_);
        IJKOffset ijkOffset;
        ijkOffset.iBlockIndex = pSharedData_->blockConfiguration.iBlockIndex;
        ijkOffset.jBlockIndex = pSharedData_->blockConfiguration.jBlockIndex;
        ijkOffset.i = iOffset;
        ijkOffset.j = jOffset;
        ijkOffset.k = kOffset;

        acc::modify_2_tuples<
            RestoreAndAdvanceInIJKFunctor,
            TIterationMask,
            const IJKOffset
        >(pSharedData_->storageStridesTuple, indexTuple_, ijkOffset);

        // set the boundary offsets
        iBoundaryOffset_ = computeBoundaryOffset(iOffset, pSharedData_->blockConfiguration.iStart, pSharedData_->blockConfiguration.iEnd);
        jBoundaryOffset_ = computeBoundaryOffset(jOffset, pSharedData_->blockConfiguration.jStart, pSharedData_->blockConfiguration.jEnd);
        
        // set the thread offsets
        iThreadOffset_ = iOffset;
        jThreadOffset_ = jOffset;
    }

    /**
    * Overload of the bracket operator which accesses a tuple element
    */
    // specialize for compile-time static Offsets
    template<
        typename TReturnType,
        int VParameterIndex,
        typename TOffset,
        CachePolicy VCachePolicy>
    __ACC__
    TReturnType AccessImpl()
    {
        // get the cache register to element map
        typedef typename Cache::IndexToCacheElementMap IndexToCacheElementMap;

        // compute the parameter type depending on the cached elements and the caching policy
        // (either get the parameter type from the cache or the normal element map)
        typedef typename boost::mpl::eval_if_c<
            (
                VCachePolicy == cUseCache &&    
                boost::mpl::has_key<IndexToCacheElementMap, boost::mpl::integral_c<int, VParameterIndex> >::value
            ),
            boost::mpl::at<IndexToCacheElementMap, boost::mpl::integral_c<int, VParameterIndex> >,
            boost::mpl::at<IndexToElementMap, boost::mpl::integral_c<int, VParameterIndex> >
        >::type ParameterType;

        // make sure we do not access a dummy storage
        BOOST_STATIC_ASSERT(!is_dummy_storage<ParameterType>::value);

        // perform the parameter type specific parameter access
        return accessImpl<
            TReturnType,
            VParameterIndex,
            TOffset,
            VCachePolicy
        >(static_cast<ParameterType*>(0));
    }

    // specialize for run-time static Offsets
    template<
        typename TReturnType,
        int VParameterIndex,
        CachePolicy VCachePolicy>
    __ACC__
    TReturnType AccessImpl(const DataParameterDynamic<VParameterIndex> offset)
    {
        //  // compute the parameter type used to choose the right access methods
        //  // (note that this is equivalent to a lookup in the IndexToElementMap in case the cache is empty)
        //  typedef typename Cache::template compute_parameter_type<
        //      VParameterIndex, 
        //      VCachePolicy
        //  >::type ParameterType;

        // get the cache register to element map
        typedef typename Cache::IndexToCacheElementMap IndexToCacheElementMap;

        // compute the parameter type depending on the cached elements and the caching policy
        // (either get the parameter type from the cache or the normal element map)
        typedef typename boost::mpl::eval_if_c<
            (
                VCachePolicy == cUseCache &&    
                boost::mpl::has_key<IndexToCacheElementMap, boost::mpl::integral_c<int, VParameterIndex> >::value
            ),
            boost::mpl::at<IndexToCacheElementMap, boost::mpl::integral_c<int, VParameterIndex> >,
            boost::mpl::at<IndexToElementMap, boost::mpl::integral_c<int, VParameterIndex> >
        >::type ParameterType;
       
        return accessImpl<
            TReturnType,
            VParameterIndex
        >(
            static_cast<ParameterType*>(0),
            offset
        );
    }

    /**
    * Method sliding the context cache
    */
    template<typename TSlideMask>
    __ACC__
    void SlideCachesImpl() 
    {
        contextCache_.template Slide<TSlideMask>(); 
    }
    
    /**
    * Method filling the cache front elements
    */
    template<
        typename TFillMask,
        typename TDomain>
    __ACC__
    void FillCachesFrontImpl() 
    {
        contextCache_.template FillFront<TFillMask, TDomain>(*this); 
    }

	/**
    * Method filling all cache elements except for the back elements
    */
    template<
        typename TFillMask,
        typename TDomain>
    __ACC__
    void FillCachesImpl() 
    {
        contextCache_.template Fill<TFillMask, TDomain>(*this); 
    }

	/**
    * Method flushing the cache back elements
    */
    template<
        typename TFlushMask,
        typename TDomain>
    __ACC__
    void FlushCachesBackImpl() 
    {
         contextCache_.template FlushBack<TFlushMask, TDomain>(*this); 
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
        contextCache_.template Flush<TFlushMask, TDomain>(*this); 
    }

    // CUDA specific implementation

    /**
    * @return true if the thread is inside the ij range
    */
    template<typename TIJRange>
    __ACC__
    bool IsInRange() const 
    {
        BOOST_STATIC_ASSERT(is_ij_range<TIJRange>::value);
        return
            iBoundaryOffset_ >= TIJRange::IMinusOffset::value && iBoundaryOffset_ <= TIJRange::IPlusOffset::value &&
            jBoundaryOffset_ >= TIJRange::JMinusOffset::value && jBoundaryOffset_ <= TIJRange::JPlusOffset::value &&
            (TIJRange::CornerPolicy::value == cComplete || iBoundaryOffset_ == 0 || jBoundaryOffset_ == 0); 
    }

    /**
    * @return calculation domain size in k dimension
    */
    __ACC__
    int kSizeImpl() const 
    { 
        ACC_ASSERT(pSharedData_);    
        return pSharedData_->kSize; 
    }

private:
    // access data fields with static offset
    template<
        typename TReturnType,
        int VParameterIndex,
        typename TOffset,
        CachePolicy VCachePolicy,
        typename TStoragePointer>
    __ACC__
    TReturnType accessImpl(DataFieldCUDAStoragePointer<TStoragePointer>*)
    {
        // define the associated strides type 
        typedef typename DataFieldCUDAStoragePointer<TStoragePointer>::StorageStridesType StorageStridesType;

        // check the shared data pointer is valid
        ACC_ASSERT(pSharedData_);

        // depending on the caching policy either access the data field using At or BypassCache
        // (note that the BypassCache method is not reading read only data fields via texture cache
        // which makes sense if we fill want to fill a register or shared memory cache without polluting the texture cache)
        if(VCachePolicy == cUseCache)
        {
            return pSharedData_->dataTuple(static_cast<boost::mpl::integral_c<int, VParameterIndex>*>(0)).At(
                indexTuple_(static_cast<StorageStridesType*>(0)), 
                pSharedData_->storageStridesTuple(static_cast<StorageStridesType*>(0)), 
                TOffset::I::value, 
                TOffset::J::value, 
                TOffset::K::value
            );
        }
        else
        {
            return pSharedData_->dataTuple(static_cast<boost::mpl::integral_c<int, VParameterIndex>*>(0)).BypassCache(
                indexTuple_(static_cast<StorageStridesType*>(0)), 
                pSharedData_->storageStridesTuple(static_cast<StorageStridesType*>(0)), 
                TOffset::I::value, 
                TOffset::J::value, 
                TOffset::K::value
            );
        }
    }

    // data fields with dynamic offset
    template<
        typename TReturnType,
        int VParameterIndex,
        typename TStoragePointer>
    __ACC__
    TReturnType accessImpl(DataFieldCUDAStoragePointer<TStoragePointer>*, const DataParameterDynamic<VParameterIndex> offset)
    {
        // define the associated strides type 
        typedef typename DataFieldCUDAStoragePointer<TStoragePointer>::StorageStridesType StorageStridesType;

        // use the at function
        ACC_ASSERT(pSharedData_);
        return pSharedData_->dataTuple(static_cast<boost::mpl::integral_c<int, VParameterIndex>*>(0)).At(
            indexTuple_(static_cast<StorageStridesType*>(0)), 
            pSharedData_->storageStridesTuple(static_cast<StorageStridesType*>(0)), 
            offset.iOffset,
            offset.jOffset,
            offset.kOffset
        );
    }

    // access k cache values with static offset
    template<
        typename TReturnType,
        int VParameterIndex,
        typename TOffset,
        CachePolicy VCachePolicy,
        typename TValue,
        typename TKWindow, 
        typename TKRange, 
        CacheIOPolicy VIOPolicy,
        KLoopDirection VKLoopDirection>
    __ACC__
    TReturnType accessImpl(KCacheStorage<TValue, TKWindow, TKRange, VParameterIndex, VIOPolicy, VKLoopDirection>*)
    {
        return contextCache_.kCacheTuple()(static_cast<boost::mpl::integral_c<int, VParameterIndex>*>(0)).At(static_cast<TOffset*>(0));
    }
    // TODO: implement access to k cache values with dynamical offset

    // access ij cache values with static offset
    template<
        typename TReturnType,
        int VParameterIndex,
        typename TOffset,
        CachePolicy VCachePolicy,
        typename TValue,
        typename TIJRange,
        typename TBlockSize>
    __ACC__
    TReturnType accessImpl(IJCacheStorage<TValue, TIJRange, TBlockSize>*)
    {
        // ij caches do not support accesses with a k offset
        BOOST_STATIC_ASSERT(TOffset::K::value == 0);
        
        // check the shared data pointer is valid
        ACC_ASSERT(pSharedData_);

        // convert the ij cache storage pointer into an ij cache tuple pointer
        IJCacheTuple* pIJCacheTuple = reinterpret_cast<IJCacheTuple*>(pSharedData_->ijCacheStorage);
        return (*pIJCacheTuple)(static_cast<boost::mpl::integral_c<int, VParameterIndex>*>(0)).At( 
            iThreadOffset_ + TOffset::I::value,
            jThreadOffset_ + TOffset::J::value
        );
    }
    // TODO: implement access to ij cache values with dynamical offset

    // access read only scalars
    template<
        typename TReturnType,
        int VParameterIndex,
        typename TOffset,
        CachePolicy VCachePolicy,
        typename TValue>
    __ACC__
    TReturnType accessImpl(ScalarStorage<TValue, cReadOnly>*)
    {
        ACC_ASSERT(pSharedData_);
        return pSharedData_->dataTuple(static_cast<boost::mpl::integral_c<int, VParameterIndex>*>(0)).value();
    }
    // TODO: implement access to scalar values with dynamical offset


    // access read write scalars
    template<
        typename TReturnType,
        int VParameterIndex,
        typename TOffset,
        CachePolicy VCachePolicy,
        typename TValue>
    __ACC__
    TReturnType accessImpl(ScalarStorage<TValue, cReadWrite>*)
    {
        return registerDataTuple_(static_cast<boost::mpl::integral_c<int, VParameterIndex>*>(0)).value();
    }
    // TODO: implement access to scalar values with dynamical offset

    // method computing the boundary offset given block start and end
    __ACC__
    int computeBoundaryOffset(const int threadOffset, const int blockStart, const int blockEnd) const
    {
        int result = 0; // return 0 inside the block
        if(threadOffset < blockStart)
        {
            result = threadOffset - blockStart; // negative values at the minus boundary
        }
        else if(threadOffset >= blockEnd) 
        {
            result = threadOffset - blockEnd + 1; // positive values at the plus boundary
        }
        return result;
    }

    // offsets with respect to the thread boundaries, 0 for threads inside the block
    int iBoundaryOffset_; 
    int jBoundaryOffset_;

    // ij thread offsets with respect to the block origin 
    int iThreadOffset_;
    int jThreadOffset_;
    
    // parameter and cache data
    acc::Tuple<IndexElements> indexTuple_;
    acc::Tuple<RegisterDataElements> registerDataTuple_;
    Cache contextCache_;
    TSharedData* pSharedData_;
};

