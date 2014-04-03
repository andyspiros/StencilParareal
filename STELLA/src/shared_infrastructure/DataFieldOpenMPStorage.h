#pragma once

#include <cassert>
#include <vector>
#include <cstring>
#include <boost/static_assert.hpp>
#include "Definitions.h"
#include "Enums.h"
#include "DataFieldStorageFormat.h"
#include "DataFieldStorageStrides.h"
#include "DataFieldOpenMPStorageIterator.h"
#include "DataFieldStorage.h"

/**
* @class DataFieldOpenMPStorage
* DataField storage class used by the OpenMP back end. It provides iterator based data access.
*/
template<
    typename TValue,
    typename TStorageFormat>
class DataFieldOpenMPStorage : public DataFieldStorage<TValue, TStorageFormat, DataFieldOpenMPStorage<TValue, TStorageFormat> > // CRTP
{
    BOOST_STATIC_ASSERT(is_data_field_storage_format<TStorageFormat>::value);
public:
    using DataFieldStorage<TValue, TStorageFormat, DataFieldOpenMPStorage>::Assign;
    using DataFieldStorage<TValue, TStorageFormat, DataFieldOpenMPStorage>::paddedSize;
    using DataFieldStorage<TValue, TStorageFormat, DataFieldOpenMPStorage>::originOffset;
    using DataFieldStorage<TValue, TStorageFormat, DataFieldOpenMPStorage>::size;

    typedef DataFieldOpenMPStorageIterator<TValue, TStorageFormat, cReadWrite> StorageIteratorType;
    typedef DataFieldOpenMPStorageIterator<TValue, TStorageFormat, cReadOnly> ConstStorageIteratorType;
    typedef typename DataFieldStorage<TValue, TStorageFormat, DataFieldOpenMPStorage>::ValueType ValueType;
    typedef typename DataFieldStorage<TValue, TStorageFormat, DataFieldOpenMPStorage>::PointerType PointerType;
    typedef typename DataFieldStorage<TValue, TStorageFormat, DataFieldOpenMPStorage>::ExternalStorageType ExternalStorageType;

    __CPU__
    DataFieldOpenMPStorage() 
    { 
        pStorageBase_ = NULL; 
        storageSize_ = 0; 
    }
    __CPU__
    ~DataFieldOpenMPStorage() 
    {
        clearStorageMemory();
    }

    __CPU__
    DataFieldOpenMPStorage(const DataFieldOpenMPStorage& other) { *this = other; }
    __CPU__
    DataFieldOpenMPStorage& operator= (const DataFieldOpenMPStorage& other)
    {
        // run the base assign
        Assign(other);
        floor_ = other.floor_;
        ceiling_ = other.ceiling_;
        storageSize_ = other.storageSize_;

        // finish initialization in case the other field was initialized
        if(storageSize_ != 0)
        {  
            // setup the storage
            initStorageMemory();

            // copy the data to the storage
            memcpy(static_cast<void*>(pStorageBase_), static_cast<void*>(other.pStorageBase_), sizeof(ValueType) * storageSize_);
           
            // setup the iterators
            InitializeStorageIteratorToOrigin(originIterator_);
            InitializeStorageIteratorToOrigin(constOriginIterator_);
        }
        else
        {
            clearStorageMemory();
        }

        // by convention
        return *this;
    }

    /**
    * Init the OpenMP storage
    */
    __CPU__
    void InitImpl()
    {
        // init the storage size
        initStorageSize();

        // allocate the storage
        initStorageMemory();

        // initialize the iterators
        initIterators();
    }

    /**
    * Init the OpenMP storage
    * @param externalStorage external storage
    */
    __CPU__
    void InitImpl(const ExternalStorageType& externalStorage)
    {        
        // init the storage size
        initStorageSize();

        // setup the storage memory
        initStorageMemory(externalStorage);
        
        // initialize the iterators
        initIterators();
    }

    /**
    * Set OpenMP storage to external storage
    * @param externalStorage external storage
    */
    __CPU__
    void SetExternalStorageImpl(const ExternalStorageType& externalStorage)
    {
        // setup the storage memory without size calculation
        initStorageMemory(externalStorage);

        // initialize the iterators
        initIterators();
    }

    /**
    * @return pointer to the base of the storage
    */
    __CPU__
    PointerType pStorageBaseImpl() const { return pStorageBase_; }

    /**
    * Set new valid range for access range check in debug mode
    * @param floor smallest allowed index
    * @param ceiling largest allowed index
    */
    __CPU__
    void SetValidRange(const IJKIndex& floor, const IJKIndex& ceiling)
    {
        // assert if valid range greater than size of data field storage
        assert(floor.iIndex() >= 0 && floor.jIndex() >= 0 && floor.kIndex() >= 0);
        assert(floor.iIndex() <= ceiling.iIndex() && floor.jIndex() <= ceiling.jIndex() && floor.kIndex() <= ceiling.kIndex());
        assert(ceiling.iIndex() < size().iSize() && ceiling.jIndex() < size().jSize() && ceiling.kIndex() < size().kSize());

        floor_ = floor;
        ceiling_ = ceiling;
        originIterator_.SetValidRange(floor_, ceiling_);
        constOriginIterator_.SetValidRange(floor_, ceiling_);
    }

    /**
    * @return an iterator pointing to the origin.
    */
    __CPU__
    const StorageIteratorType& originIterator() const { return originIterator_; }

    /**
    * @return an iterator pointing to the origin.
    */
    __CPU__
    const ConstStorageIteratorType& constOriginIterator() const { return constOriginIterator_; }

    /**
    * @param iterator iterator which is initialized to the origin position of the data field
    */
    template<AccessPolicy VAccessPolicy>
    __CPU__
    void InitializeStorageIteratorToOrigin(DataFieldOpenMPStorageIterator<TValue, TStorageFormat, VAccessPolicy>& iterator) const
    {
        iterator.Init(pStorageBase_, paddedSize());
        iterator.Advance(originOffset().iIndex(), originOffset().jIndex(), originOffset().kIndex());
        iterator.MemorizePosition(iterator); // set restore point
        iterator.SetValidRange(floor_, ceiling_);
    }
    
    /**
    * Memorize the origin iterator position
    */
    __CPU__
    void MemorizePositionOfOrigin(StorageIteratorType& iterator) const
    {
        iterator.MemorizePosition(originIterator());
    }

    /**
    * Memorize the const origin iterator position
    */
    __CPU__
    void MemorizePositionOfOrigin(ConstStorageIteratorType& iterator) const
    {
        iterator.MemorizePosition(constOriginIterator());
    }

private:
    // method freeing private memory of the storage
    __CPU__
    void clearStorageMemory()
    {
        storageMemory_.clear();
        pStorageBase_ = NULL;
    }

    // setup the size information return the origin stride
    __CPU__
    void initStorageSize()
    {
        // allocate the memory (first calculate the size using the storage size functions)
        storageSize_ = paddedSize().iSize() * paddedSize().jSize() * paddedSize().kSize();
        assert(storageSize_ != 0);
    }

    // compute the strides of the origin
    __CPU__
    int computeOriginStride()
    {
        // init the strides
        DataFieldStorageStrides<typename TStorageFormat::StorageOrder> storageStrides;
        storageStrides.Init(paddedSize());

        // compute origin stride
        return storageStrides.ComputeStride(originOffset().iIndex(), originOffset().jIndex(), originOffset().kIndex());
    }

    // method setting up storage memory
    __CPU__
    void initStorageMemory()
    {
        // verify the storageSize_ variable is up to date
        assert(storageSize_ == paddedSize().iSize() * paddedSize().jSize() * paddedSize().kSize());

        // allocate the memory
        storageMemory_.resize(storageSize_ + storage_alignment<TStorageFormat>::value);

        // setup the pointers
        int originStride = computeOriginStride();
        pStorageBase_ = aligned_pointer<TStorageFormat>(&storageMemory_[0] + originStride) - originStride;
    }

    // method setting up storage memory using external storage
    __CPU__
    void initStorageMemory(const ExternalStorageType& externalStorage)
    {
        // verify the external storage size
        assert(paddedSize().iSize() == externalStorage.paddedSize().iSize());
        assert(paddedSize().jSize() == externalStorage.paddedSize().jSize());
        assert(paddedSize().kSize() == externalStorage.paddedSize().kSize());

        // verify the storageSize_ variable is up to date
        assert(storageSize_ == paddedSize().iSize() * paddedSize().jSize() * paddedSize().kSize());

        // free all internal memory of the storage
        clearStorageMemory();

        // set the storage pointer to the external memory
        pStorageBase_ = externalStorage.pStoragePointer();

        // make sure the external memory is aligned
        int originStride = computeOriginStride();
        assert(pStorageBase_ == aligned_pointer<TStorageFormat>(pStorageBase_ + originStride) - originStride);           
    }

    // method initializing iterators
    __CPU__
    void initIterators()
    {
        // computed the valid range
        floor_.Init(0, 0, 0);
        ceiling_.Init(size().iSize() - 1, size().jSize() - 1, size().kSize() - 1);

        // initialize the iterator and move it to the origin
        InitializeStorageIteratorToOrigin(originIterator_);
        InitializeStorageIteratorToOrigin(constOriginIterator_);
    }

    IJKIndex floor_;
    IJKIndex ceiling_;

    PointerType pStorageBase_;
    std::vector<ValueType> storageMemory_; // vector used to allocate memory if no external memory is provided
    int storageSize_;

    StorageIteratorType originIterator_;
    ConstStorageIteratorType constOriginIterator_;
};

  
