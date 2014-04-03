#pragma once

#include <cassert>
#include <algorithm>
#include <boost/static_assert.hpp>
#include <boost/type_traits/add_const.hpp>
#include <boost/mpl/bool.hpp>
#include "KBoundary.h"
#include "DataField.h"
#include "DataFieldStorageFormat.h"
#include "DataFieldOpenMPStorage.h"
#include "DataFieldOpenMPStorageIterator.h"

/**
* @class DataFieldOpenMP
* OpenMP implementation of the data field class
*/
template<
    typename TValue,
    typename TStorageFormat>
class DataFieldOpenMP : public DataField<DataFieldOpenMPStorage<TValue, TStorageFormat>, DataFieldOpenMP<TValue, TStorageFormat> > // CRTP interface class
{
    BOOST_STATIC_ASSERT(is_data_field_storage_format<TStorageFormat>::value);
public:
    using DataField<DataFieldOpenMPStorage<TValue, TStorageFormat>, DataFieldOpenMP>::Assign;
    using DataField<DataFieldOpenMPStorage<TValue, TStorageFormat>, DataFieldOpenMP>::isInitialized;
    using DataField<DataFieldOpenMPStorage<TValue, TStorageFormat>, DataFieldOpenMP>::calculationDomain;

    typedef typename DataField<DataFieldOpenMPStorage<TValue, TStorageFormat>, DataFieldOpenMP>::StorageType StorageType;
    typedef typename DataField<DataFieldOpenMPStorage<TValue, TStorageFormat>, DataFieldOpenMP>::ValueType ValueType;
    typedef typename DataField<DataFieldOpenMPStorage<TValue, TStorageFormat>, DataFieldOpenMP>::PointerType PointerType;
    typedef typename DataField<DataFieldOpenMPStorage<TValue, TStorageFormat>, DataFieldOpenMP>::ReferenceType ReferenceType;
    typedef typename DataField<DataFieldOpenMPStorage<TValue, TStorageFormat>, DataFieldOpenMP>::ConstReferenceType ConstReferenceType;
    typedef typename DataField<DataFieldOpenMPStorage<TValue, TStorageFormat>, DataFieldOpenMP>::ExternalStorageType ExternalStorageType;

    DataFieldOpenMP() { pStorage_ = new StorageType(); }
    ~DataFieldOpenMP() { delete pStorage_; }

    DataFieldOpenMP(const DataFieldOpenMP& other) 
    { 
        // allocate storage here
        pStorage_ = new StorageType();
        *this = other; 
    }
    DataFieldOpenMP& operator= (const DataFieldOpenMP& other) 
    {
        // assign the base values
        Assign(other);

        // copy the storage using its copy operator
        *pStorage_ = *(other.pStorage_);
        
        // by convention
        return *this;
    }   

    /** 
    * Implement the OpenMP data field initialization
    * @param kBoundary k boundary
    */
    void InitImpl(const KBoundary& kBoundary)
    {   
        // make sure the implementation is called via DataField
        assert(isInitialized());

        // initialize data field storage
        pStorage_->Init(calculationDomain(), kBoundary);
    }

    /** 
    * Implement the OpenMP data field initialization
    * @param externalStorage external storage    
    * @param kBoundary k boundary
    */
    void InitImpl(const ExternalStorageType& externalStorage, const KBoundary& kBoundary)
    {   
        // make sure the implementation is called via DataField
        assert(isInitialized());

        // initialize data field storage
        pStorage_->Init(externalStorage, calculationDomain(), kBoundary);
    }

    /**
    * Implement the OpenMP set external storage operation
    * @param externalStorage external storage
    */
    __CPU__
    void SetExternalStorageImpl(const ExternalStorageType& externalStorage)
    {
        // set the external storage on the storage
        pStorage_->SetExternalStorage(externalStorage);
    }

    /** 
    * Implement the OpenMP swap with operation 
    */
    void SwapWithImpl(DataFieldOpenMP& other)
    {
        // swap the storage pointers
        std::swap(pStorage_, other.pStorage_);
    }
    
    /** 
    * Implement the OpenMP access operator
    * @return value at i j k offset
    */
    ReferenceType AccessImpl(const int i, const int j, const int k)
    {
        return pStorage_->originIterator().At(i, j, k);
    }

    /** 
    * Implement the OpenMP access operator
    * @return value at i j k offset
    */
    ConstReferenceType AccessImpl(const int i, const int j, const int k) const
    {
        return pStorage_->originIterator().At(i, j, k);
    }

    /**
    * @return reference to data field storage buffer
    */
    const StorageType& storageImpl() const { return *pStorage_; }
    
private:
    StorageType* pStorage_;
};

// add data field trait
template<
    typename TValue,
    typename TStorageFormat>
struct is_data_field<DataFieldOpenMP<TValue, TStorageFormat> > : boost::mpl::true_ {};

/**
* @struct create_storage_iterator
* Meta function creating a iterator given a data field
*/
template<
    typename TDataField,
    ParameterIntend VParameterIntend>
struct create_storage_iterator;

// for in out fields create a read write iterator
template<
    typename TValue,
    typename TStorageFormat>
struct create_storage_iterator<DataFieldOpenMP<TValue, TStorageFormat>, cInOut>
{
    typedef DataFieldOpenMPStorageIterator<TValue, TStorageFormat, cReadWrite> type;
};

// for in fields create a read only iterator
template<
    typename TValue,
    typename TStorageFormat>
struct create_storage_iterator<DataFieldOpenMP<TValue, TStorageFormat>, cIn>
{
    typedef DataFieldOpenMPStorageIterator<TValue, TStorageFormat, cReadOnly> type;
};
