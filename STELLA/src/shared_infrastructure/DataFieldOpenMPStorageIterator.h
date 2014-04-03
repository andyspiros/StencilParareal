#pragma once

#include <cassert>
#include <boost/static_assert.hpp>
#include <boost/type_traits/add_const.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/bool.hpp>
#include "Enums.h"
#include "ParameterTraits.h"
#include "PositionChecker.h"
#include "IJKSize.h"
#include "DataFieldStorageStrides.h"

/**
* @class DataFieldStorageIterator
* The DataFieldStorageIterator stores a pointer to data field storage.
* It can be dereference and moved in 3 dimensions. 
* In case the underlying data field does not cover a dimension an advance is just ignored.
*/
template<
    typename TValue,
    typename TStorageFormat,
    AccessPolicy VAccessPolicy>
class DataFieldOpenMPStorageIterator
{
public:    
    typedef TValue ValueType;
    typedef TValue* PointerType;
    typedef TValue& ReferenceType;
    typedef typename boost::add_const<ValueType>::type& ConstReferenceType;
    typedef typename boost::mpl::if_c<
        VAccessPolicy == cReadWrite, 
        ReferenceType,
        ConstReferenceType
    >::type ReturnType;

    DataFieldOpenMPStorageIterator()
    {
        pCurrentPosition_ = NULL;
        pMemorizedPosition_ = NULL;
    }
    ~DataFieldOpenMPStorageIterator() {}
    
    DataFieldOpenMPStorageIterator(const DataFieldOpenMPStorageIterator& other) { *this = other; }
    DataFieldOpenMPStorageIterator& operator= (const DataFieldOpenMPStorageIterator& other)
    {
        // copy the members
#ifndef NDEBUG
        positionChecker_ = other.positionChecker_;
        memorizedPositionChecker_ = other.memorizedPositionChecker_;
#endif
        pCurrentPosition_ = other.pCurrentPosition_;
        pMemorizedPosition_ = other.pMemorizedPosition_;
        storageStrides_ = other.storageStrides_;
        // by convention, always return *this
        return *this;
    }

    /**
    * Init the iterator to the position pStorage and set store position to pStorage.
    * @param pStorage pointer to the data field storage
    * @param size Size of the storage, includes calculation domain and boundary
    */
    void Init(PointerType pStorage, const IJKSize& size)
    {
#ifndef NDEBUG
        IJKIndex floor, ceiling;
        floor.Init(0,0,0);
        ceiling.Init(size.iSize() - 1, size.jSize() - 1, size.kSize() - 1);
        positionChecker_.Init(floor, ceiling, 
            has_storage_in<TStorageFormat, cDimI>::value, 
            has_storage_in<TStorageFormat, cDimJ>::value, 
            has_storage_in<TStorageFormat, cDimK>::value
        );
        memorizedPositionChecker_ = positionChecker_;
#endif
        // set the current position and copy the stride information
        pCurrentPosition_ = pStorage;
        pMemorizedPosition_ = pCurrentPosition_;
        storageStrides_.Init(size);
    }
    
    /**
    * Set new valid range for access range check in debug mode
    * @param floor smallest allowed index
    * @param ceiling largest allowed index
    */
    void SetValidRange(const IJKIndex& floor, const IJKIndex& ceiling)
    {
#ifndef NDEBUG
        positionChecker_.SetValidRange(floor, ceiling);
        memorizedPositionChecker_.SetValidRange(floor, ceiling);
#endif
    }
    
    /**
    * Advance the iterator in all three dimensions
    * @param i offset in i direction
    * @param j offset in j direction
    * @param k offset in k direction
    */
    void Advance(const int i, const int j, const int k)
    {
        assert(pCurrentPosition_);
#ifndef NDEBUG
        positionChecker_.Advance(i, j, k);
#endif
        pCurrentPosition_ += storageStrides_.ComputeStride(i,j,k);
    }
   
    /**
    * Dereference the iterator at an offset in all three dimensions
    * @param i offset in i direction
    * @param j offset in j direction
    * @param k offset in k direction
    * @return reference to the value at the offset position
    */
    ReturnType At(const int i, const int j, const int k) const
    {
        assert(pCurrentPosition_);
#ifndef NDEBUG
        positionChecker_.CheckPositionAtOffset(i, j, k);
#endif
        return *(pCurrentPosition_ + storageStrides_.ComputeStride(i,j,k));
    }
    
    /**
    * Memorize the position of the iterator handed over as parameter.
    * @param iterator containing position to memorize
    */
    void MemorizePosition(const DataFieldOpenMPStorageIterator& iterator)
    {
        assert(iterator.pCurrentPosition_ != NULL);
#ifndef NDEBUG
        memorizedPositionChecker_ = iterator.positionChecker_;
        // reset the bounding box of the memorized position
        memorizedPositionChecker_.ResetBoundingBox();
#endif
        pMemorizedPosition_ = iterator.pCurrentPosition_;
    }

    /**
    * Advance memorized position in all three dimensions.
    * @param i offset in i direction
    * @param j offset in j direction
    * @param k offset in k direction
    */
    void AdvanceMemorizedPosition(const int i, const int j, const int k)
    {
        assert(pMemorizedPosition_);
#ifndef NDEBUG
        memorizedPositionChecker_.Advance(i, j, k);
#endif
        pMemorizedPosition_ += storageStrides_.ComputeStride(i,j,k);
    }

    /**
    * Set current position to the memorized position
    */
    void RestoreMemorizedPosition()
    {
        assert(pMemorizedPosition_);
#ifndef NDEBUG
        // merge the bounding box with the restored position
        memorizedPositionChecker_.MergeBoundingBox(positionChecker_);
        // setup the new position checker
        positionChecker_ = memorizedPositionChecker_;
#endif
        pCurrentPosition_ = pMemorizedPosition_;
    }

#ifndef NDEBUG
    /**
    * @return the position checker
    */
    PositionChecker& positionChecker() const { return positionChecker_; }

    /**
    * @return the memorized position checker
    */
    PositionChecker& memorizedPositionChecker() const { return memorizedPositionChecker_; }
#endif
       
private:
    PointerType pCurrentPosition_;
    PointerType pMemorizedPosition_;
    DataFieldStorageStrides<typename TStorageFormat::StorageOrder> storageStrides_; 
#ifndef NDEBUG
    mutable PositionChecker positionChecker_;
    mutable PositionChecker memorizedPositionChecker_;
#endif
};

// return type specialization
template<
    typename TValue,
    typename TStorageFormat,
    AccessPolicy VAccessPolicy>
struct return_type<DataFieldOpenMPStorageIterator<TValue, TStorageFormat, VAccessPolicy> >
{
    typedef typename DataFieldOpenMPStorageIterator<TValue, TStorageFormat, VAccessPolicy>::ReturnType type;
};

// value type specialization
template<
    typename TValue,
    typename TStorageFormat,
    AccessPolicy VAccessPolicy>
struct value_type<DataFieldOpenMPStorageIterator<TValue, TStorageFormat, VAccessPolicy> >
{
    typedef typename DataFieldOpenMPStorageIterator<TValue, TStorageFormat, VAccessPolicy>::ReturnType type;
};

// is iterable specialization
template<
    typename TValue,
    typename TStorageFormat,
    AccessPolicy VAccessPolicy>
struct is_iterable<DataFieldOpenMPStorageIterator<TValue, TStorageFormat, VAccessPolicy> > : boost::mpl::true_ {};
