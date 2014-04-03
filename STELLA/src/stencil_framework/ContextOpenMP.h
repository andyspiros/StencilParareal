#pragma once

#include <boost/mpl/fold.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/and.hpp>
#include "SharedInfrastructure.h"
#include "LaunchConfiguration.h"
#include "Context.h"
#include "ContextOpenMPFunctors.h"
#include "DataFieldOpenMPStorageIterator.h"

/**
* @class ContextOpenMP
* Class implementing the context functionality for the OpenMP back end
*/
template<typename TContextDescriptor>
class ContextOpenMP : public Context<TContextDescriptor, ContextOpenMP<TContextDescriptor> > // CRTP
{
public:
    // compute the tuple type
    typedef typename TContextDescriptor::TupleElements TupleElements;
    typedef Tuple<TupleElements> DataTupleType;
    typedef typename TContextDescriptor::IndexToElementMap IndexToElementMap;

    // function computing the iteration mask holding all parameter indexes 
    // used by the stencil stages of the sweep
    template<typename TStencilStages, typename TDomain>
    struct compute_iteration_mask
    {
        typedef typename boost::mpl::fold<
            typename TupleElements::ElementIndexes,
            boost::mpl::set0<>,
            boost::mpl::if_<
                boost::mpl::and_<
                    stencil_stages_have_parameter<TStencilStages, TDomain, boost::mpl::_2>,
                    is_iterable<boost::mpl::at<IndexToElementMap, boost::mpl::_2> >
                >,
                boost::mpl::insert<boost::mpl::_1, boost::mpl::_2>,
                boost::mpl::_1
            >
        >::type type;
    };

    ContextOpenMP() 
    { 
        kSize_ = 0; 
        iBlockStart_ = 0;
        iBlockEnd_ = 0;
        jBlockStart_ = 0;
        jBlockEnd_ = 0;
    }
    ~ContextOpenMP() {}
    
    ContextOpenMP(const ContextOpenMP& other) { *this = other; }
    ContextOpenMP& operator= (const ContextOpenMP& other)
    {
        kSize_ = other.kSize_; 
        iBlockStart_ = other.iBlockStart_;
        iBlockEnd_ = other.iBlockEnd_;
        jBlockStart_ = other.jBlockStart_;
        jBlockEnd_ = other.jBlockEnd_;
        dataTuple_ = other.dataTuple_;
        return *this;
    }

    /**
    * Method initializing the context by setting the kSize
    * @param kSize calculation domain size in k dimension
    */
    void Init(const int kSize) { kSize_ = kSize; }

    /**
    * Advance the context position by a constant offset
    * The iteration mask allows to apply the method to a subset of the context
    */
    template<typename TIterationMask, int VI, int VJ, int VK>
    void AdvanceImpl()
    {
        modify_tuple<
            AdvanceFunctor<VI,VJ,VK>,
            TIterationMask
        >(dataTuple_);
    }

    /**
    * Advance the context position by a variable a k offset
    * The iteration mask allows to apply the method to a subset of the context
    * @param kOffset kOffset relative to the current position
    */
    template<typename TIterationMask, int VI, int VJ>
    void AdvanceImpl(const int kOffset)
    {
        modify_tuple<
            AdvanceInKFunctor<VI,VJ>,
            TIterationMask,
            const int
        >(dataTuple_, kOffset);
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
        IJOffset ijOffset;
        ijOffset.i = iOffset;
        ijOffset.j = jOffset;

        modify_tuple<
            RestoreAndAdvanceInIJFunctor<VK>,
            TIterationMask,
            const IJOffset
        >(dataTuple_, ijOffset);
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
        IJKOffset ijkOffset;
        ijkOffset.i = iOffset;
        ijkOffset.j = jOffset;
        ijkOffset.k = kOffset; 

        modify_tuple<
            RestoreAndAdvanceInIJKFunctor,
            TIterationMask,
            const IJKOffset
        >(dataTuple_, ijkOffset);
    }

    /**
    * Method sliding the context cache
    */
    template<typename TSlideMask>
    __ACC__
    void SlideCachesImpl() 
    {
        // empty method as the OpenMP back end supports no caching
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
        // empty method as the OpenMP back end supports no caching
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
        // empty method as the OpenMP back end supports no caching
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
        // empty method as the OpenMP back end supports no caching
    }

	/**
    * Method flushing all cache elements except for the back elements
    */
    template<
        typename TFlushMask,
        typename TDomain>
    __ACC__
    void FlushCachesImpl() 
    {
        // empty method as the OpenMP back end supports no caching
    }
 
    /**
    * Overload of the bracket operator which accesses a tuple element.
    * Static (compile time) offset
    */
    template<
        typename TReturnType,
        int VParameterIndex,
        typename TOffset,
        CachePolicy VCachePolicy>
    TReturnType AccessImpl()
    {
        return accessImpl<
            TReturnType,
            VParameterIndex,
            TOffset
        >(
            static_cast<        
                typename boost::mpl::at<
                    IndexToElementMap, 
                    boost::mpl::integral_c<int, VParameterIndex> 
                >::type*
            >(0)
        );
    }

    /**
    * Overload of the bracket operator which accesses a tuple element.
    * Dynamic (run time) offset
    */
    template<
        typename TReturnType,
        int VParameterIndex,
        CachePolicy VCachePolicy>
    TReturnType AccessImpl(const DataParameterDynamic<VParameterIndex> offset)
    {
        return accessImpl<
            TReturnType,
            VParameterIndex
        >(
            static_cast<        
                typename boost::mpl::at<
                    IndexToElementMap, 
                    boost::mpl::integral_c<int, VParameterIndex> 
                >::type*
            >(0),
            offset
        );
    }

    // OpenMP specific interface

    /**
    * Method setting the context block configuration
    * @param blockConfiguration block size parameters
    */
    void SetBlockConfiguration(const BlockConfiguration& blockConfiguration) 
    { 
        iBlockStart_ = blockConfiguration.iStart;
        iBlockEnd_ = blockConfiguration.iEnd;
        jBlockStart_ = blockConfiguration.jStart;
        jBlockEnd_ = blockConfiguration.jEnd;
    }

    /**
    * Advance the memorized position of the context by an i offset
    * The parameter mask allows to apply the method to a subset of the context
    * @param iOffset offset in i direction
    */
    template<typename TIterationMask>
    void AdvanceMemorizedPositionInI(const int iOffset)
    {
        modify_tuple<
            AdvanceMemorizedPositionInIFunctor,
            TIterationMask,
            const int
        >(dataTuple_, iOffset);
    }

    /**
    * Advance the memorized offset of the context by an j offset
    * The parameter mask allows to apply the method to a subset of the context
    * @param jOffset offset in j direction
    */
    template<typename TIterationMask>
    void AdvanceMemorizedPositionInJ(const int jOffset)
    {
        modify_tuple<
            AdvanceMemorizedPositionInJFunctor,
            TIterationMask,
            const int
        >(dataTuple_, jOffset);
    }

    /**
    * @return calculation domain size in k dimension
    */
    int kSizeImpl() const { return kSize_; }

    /**
    * @return i block start offset relative to the block origin
    */
    int iBlockStart() const { return iBlockStart_; }

    /**
    * @return i block end offset relative to the block origin
    */
    int iBlockEnd() const { return iBlockEnd_; }
    
    /**
    * @return j block start offset relative to the block origin
    */
    int jBlockStart() const { return jBlockStart_; }
    
    /**
    * @return j block end offset relative to the block origin
    */
    int jBlockEnd() const { return jBlockEnd_; }

    /**
    * @return tuple containing context data
    */
    DataTupleType& dataTuple() { return dataTuple_; }

private:
    // access methods to read data fields and scalars

    // static offset
    template<
        typename TReturnType,
        int VParameterIndex,
        typename TOffset,
        typename TValue,
        typename TStorageStrides, 
        AccessPolicy VAccessPolicy>
    TReturnType accessImpl(DataFieldOpenMPStorageIterator<TValue, TStorageStrides, VAccessPolicy>*)
    {
        return dataTuple_(static_cast<boost::mpl::integral_c<int, VParameterIndex>*>(0)).At(TOffset::I::value, TOffset::J::value, TOffset::K::value);    
    }

    // dynamic offset
    template<
        typename TReturnType,
        int VParameterIndex,
        typename TValue,
        typename TStorageStrides, 
        AccessPolicy VAccessPolicy>
    TReturnType accessImpl(
            DataFieldOpenMPStorageIterator<TValue, TStorageStrides, VAccessPolicy>*,
            const DataParameterDynamic<VParameterIndex> offset)
    {
        return dataTuple_(
                static_cast< boost::mpl::integral_c<int, VParameterIndex>*>(0)).At(
                    offset.iOffset, offset.jOffset, offset.kOffset
        );
    }

    template<
        typename TReturnType,
        int VParameterIndex,
        typename TOffset,
        typename TValue, 
        AccessPolicy VAccessPolicy>
    TReturnType accessImpl(ScalarStorage<TValue, VAccessPolicy>*)
    {
        return dataTuple_(static_cast<boost::mpl::integral_c<int, VParameterIndex>*>(0)).value();
    }


    // TODO: no need to specialize for dynamic offsets used to access scalars?

    // k column size
    int kSize_;
    // block start / stop indexes
    int iBlockStart_;
    int iBlockEnd_;
    int jBlockStart_;
    int jBlockEnd_;
    // parameters
    DataTupleType dataTuple_;
};
  
