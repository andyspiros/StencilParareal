#pragma once

#include <cassert>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/contains.hpp>
#include "Definitions.h"
#include "Enums.h"
#include "IJKSize.h"
#include "IJKIndex.h"
#include "KBoundary.h"
#include "ExternalStorage.h"
#include "DataFieldStorageFormat.h"

/**
* @class DataFieldStorage
* Universal storage class providing data field memory
*/
template<
    typename TValue,
    typename TStorageFormat,
    typename TDataFieldStorageImpl>
class DataFieldStorage
{
    DISALLOW_COPY_AND_ASSIGN(DataFieldStorage);
protected:
    __CPU__
    DataFieldStorage() {}
    __CPU__
    ~DataFieldStorage() {}

    // protect assign method used by implementation for assignment operator
    __CPU__
    void Assign(const DataFieldStorage& other)
    {
        size_ = other.size_;
        allocatedSize_ = other.allocatedSize_;
        paddedSize_ = other.paddedSize_;
        originOffset_ = other.originOffset_;
    }

public:
    // define the data field storage interface
    typedef TValue ValueType;
    typedef TValue* PointerType;
    typedef TStorageFormat StorageFormat;
    typedef ExternalStorage<ValueType> ExternalStorageType;   

    /**
    * Init the data field storage
    * @param calculationDomain size of the calculation domain
    * @param kBoundary boundary in k dimension
    */
    __CPU__
    void Init(const IJKSize& calculationDomain, const KBoundary& kBoundary)
    {   
        // init the size information
        initDataFieldStorage(calculationDomain, kBoundary);

        // call the init implementation
        static_cast<TDataFieldStorageImpl*>(this)->InitImpl();
    }

    /**
    * Init the data field storage using an external storage
    * @param externalStorage external storage
    * @param calculationDomain size of the calculation domain
    * @param kBoundary boundary in k dimension
    */
    __CPU__
    void Init(const ExternalStorageType& externalStorage, const IJKSize& calculationDomain, const KBoundary& kBoundary)
    {   
        // init the size information
        initDataFieldStorage(calculationDomain, kBoundary);

        // call the init implementation
        static_cast<TDataFieldStorageImpl*>(this)->InitImpl(externalStorage);
    }

    /**
    * Set data field storage to external storage
    * @param externalStorage external storage
    */
    __CPU__
    void SetExternalStorage(const ExternalStorageType& externalStorage)
    {
        // call the set external storage implementation
        static_cast<TDataFieldStorageImpl*>(this)->SetExternalStorageImpl(externalStorage);
    }

    /**
    * @return pointer to the beginning of the data field storage
    */
    __CPU__
    PointerType pStorageBase() const 
    { 
        return static_cast<const TDataFieldStorageImpl*>(this)->pStorageBaseImpl(); 
    }

    /**
    * @return origin offset in i, j, k direction
    */
    __CPU__
    const IJKIndex& originOffset() const { return originOffset_; }

    /**
    * @return total size
    */
    __CPU__
    const IJKSize& size() const { return size_; }

    /**
    * @return allocated size (depending on the rank)
    */
    __CPU__
    const IJKSize& allocatedSize() const { return allocatedSize_; }

    /**
    * @return padded size (depending on the rank)
    */
    __CPU__
    const IJKSize& paddedSize() const { return paddedSize_; }
  
    /**
    * @return rank (~dimensionality) of the field
    */
    __CPU__
    int rank() const { return storage_rank<TStorageFormat>::value; }   

private:
    // setup the size information
    __CPU__
    void initDataFieldStorage(const IJKSize& calculationDomain, const KBoundary& kBoundary)
    {
        // assert that sizes are not negative
        assert(calculationDomain.iSize() > 0 && calculationDomain.jSize() > 0 && calculationDomain.kSize() > 0);
        assert(calculationDomain.kSize() - kBoundary.kMinusOffset() + kBoundary.kPlusOffset() > 0);

        // define the ij-boundary
        typedef typename TStorageFormat::IJBoundary IJBoundary;

        // compute size and origin index of the storage
        size_.Init(
            calculationDomain.iSize() - IJBoundary::IMinusOffset::value + IJBoundary::IPlusOffset::value,
            calculationDomain.jSize() - IJBoundary::JMinusOffset::value + IJBoundary::JPlusOffset::value,
            calculationDomain.kSize() - kBoundary.kMinusOffset() + kBoundary.kPlusOffset()
        );
        originOffset_.Init(
            -IJBoundary::IMinusOffset::value,
            -IJBoundary::JMinusOffset::value,
            -kBoundary.kMinusOffset()
        );

        // compute the allocated size
        allocatedSize_.Init(
            allocation_size<TStorageFormat, cDimI>(size_.iSize()),
            allocation_size<TStorageFormat, cDimJ>(size_.jSize()),
            allocation_size<TStorageFormat, cDimK>(size_.kSize())
        );

        // compute the padded size
        paddedSize_.Init(
            padded_size<TStorageFormat, cDimI>(allocatedSize_.iSize()),
            padded_size<TStorageFormat, cDimJ>(allocatedSize_.jSize()),
            padded_size<TStorageFormat, cDimK>(allocatedSize_.kSize())
        );
    }

    IJKSize size_;
    IJKSize allocatedSize_;
    IJKSize paddedSize_;
    IJKIndex originOffset_;
};

  
