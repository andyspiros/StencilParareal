#pragma once

#include <boost/static_assert.hpp>
#include <boost/type_traits/add_const.hpp>
#include "KBoundary.h"
#include "DataFieldStorageFormat.h"
#include "DataField.h"
#include "DataFieldCUDAStorage.h"
#include "DataFieldCUDAStoragePointer.h"

/**
* @class DataFieldCUDA
* CUDA implementation of the data field class
*/
template<
    typename TValue,
    typename TStorageFormat>
class DataFieldCUDA : public DataField<DataFieldCUDAStorage<TValue, TStorageFormat>, DataFieldCUDA<TValue, TStorageFormat> > // CRTP 
{
    BOOST_STATIC_ASSERT(is_data_field_storage_format<TStorageFormat>::value);
public:
    using DataField<DataFieldCUDAStorage<TValue, TStorageFormat>, DataFieldCUDA>::Assign;
    using DataField<DataFieldCUDAStorage<TValue, TStorageFormat>, DataFieldCUDA>::isInitialized;
    using DataField<DataFieldCUDAStorage<TValue, TStorageFormat>, DataFieldCUDA>::calculationDomain;
    using DataField<DataFieldCUDAStorage<TValue, TStorageFormat>, DataFieldCUDA>::boundary;
    
    typedef typename DataField<DataFieldCUDAStorage<TValue, TStorageFormat>, DataFieldCUDA>::PointerType PointerType;
    typedef typename DataField<DataFieldCUDAStorage<TValue, TStorageFormat>, DataFieldCUDA>::ReferenceType ReferenceType;
    typedef typename DataField<DataFieldCUDAStorage<TValue, TStorageFormat>, DataFieldCUDA>::ConstReferenceType ConstReferenceType;
    typedef typename DataField<DataFieldCUDAStorage<TValue, TStorageFormat>, DataFieldCUDA>::ExternalStorageType ExternalStorageType;
    typedef typename DataField<DataFieldCUDAStorage<TValue, TStorageFormat>, DataFieldCUDA>::StorageType DeviceStorageType;
    typedef DataFieldOpenMPStorage<TValue, TStorageFormat> HostStorageType;

    __CPU__
    DataFieldCUDA() 
    { 
        pDeviceStorage_ = new DeviceStorageType(); 
        pHostStorage_ = NULL;
    }
    __CPU__
    ~DataFieldCUDA() 
    { 
        assert(pDeviceStorage_);
        delete pDeviceStorage_; 
    }

    __CPU__
    DataFieldCUDA(const DataFieldCUDA& other) 
    { 
        // allocate storage here
        pDeviceStorage_ = new DeviceStorageType(); 
        pHostStorage_ = NULL;
        *this = other; 
    }
    __CPU__
    DataFieldCUDA& operator= (const DataFieldCUDA& other) 
    {
        // make sure the data is synchronized back to the device storage
        // it could be an external storage with a lifetime that exceeds the data field lifetime
        SynchronizeDeviceStorage();

        // assign the base values
        Assign(other);

        // always set the device storage even in case its content is not valid
        // the dimensions have to be copied nevertheless
        *pDeviceStorage_ = *(other.pDeviceStorage_);
       
        // either create a copy of the host storage or set it to null
        if(other.pHostStorage_)
        {
            pHostStorage_ = new HostStorageType();
            *pHostStorage_ = *(other.pHostStorage_);
        }
        else
        {   
            pHostStorage_ = NULL;
        }

        // by convention
        return *this;
    } 

    /**
    * Implement the CUDA data field initialization
    * @param kBoundary k boundary
    */
    __CPU__
    void InitImpl(const KBoundary& kBoundary)
    {
        // make sure the implementation is called via DataField
        assert(isInitialized());

        // make sure the data is synchronized back to the device storage
        // it could be an external storage with a lifetime that exceeds the data field lifetime
        SynchronizeDeviceStorage();

        // initialize the device storage
        pDeviceStorage_->Init(calculationDomain(), kBoundary);
    }

    /**
    * Implement the CUDA data field initialization
    * @param externalStorage external storage    
    * @param kBoundary k boundary
    */
    __CPU__
    void InitImpl(const ExternalStorageType& externalStorage, const KBoundary& kBoundary)
    {
        // make sure the implementation is called via DataField
        assert(isInitialized());

        // make sure the data is synchronized back to the device storage
        // it could be an external storage with a lifetime that exceeds the data field lifetime
        SynchronizeDeviceStorage();

        // initialize the device storage
        pDeviceStorage_->Init(externalStorage, calculationDomain(), kBoundary);
    }
    
    /**
    * Implement the CUDA set external storage operation
    * @param externalStorage external storage
    */
    __CPU__
    void SetExternalStorageImpl(const ExternalStorageType& externalStorage)
    {
        // make sure the data is synchronized back to the device storage
        // it could be an external storage with a lifetime that exceeds the data field lifetime
        SynchronizeDeviceStorage();

        // initialize the device storage with external memory
        pDeviceStorage_->SetExternalStorage(externalStorage);
    }

    /**
    * Implement the CUDA swap with operation 
    */
    __CPU__
    void SwapWithImpl(DataFieldCUDA& other)
    {
        // swap the storage pointers
        std::swap(pDeviceStorage_, other.pDeviceStorage_);
        std::swap(pHostStorage_, other.pHostStorage_);
    }

    /** 
    * Implement the CUDA access operator
    * @return value at i j k offset
    */
    __CPU__
    ReferenceType AccessImpl(const int i, const int j, const int k)
    {
        SynchronizeHostStorage();
        return pHostStorage_->originIterator().At(i, j, k);
    }

    /** 
    * Implement the CUDA access operator
    * @return value at i j k offset
    */
    __CPU__
    ConstReferenceType AccessImpl(const int i, const int j, const int k) const
    {
        SynchronizeHostStorage();
        return pHostStorage_->originIterator().At(i, j, k);
    }

    /**
    * @return reference to the device storage
    */
    __CPU__
    const DeviceStorageType& storageImpl() const { return *pDeviceStorage_; }

    /**
    * Copy the device to the host storage if necessary
    */
    __CPU__
    void SynchronizeHostStorage() const
    {
        // create a host copy and update it with the device storage
        if(!pHostStorage_)
        {
            initHostStorage();
            pDeviceStorage_->CopyDeviceToHostStorage(*pHostStorage_);
        }
    }

    /**
    * Copy the host to the device storage if necessary
    */
    __CPU__
    void SynchronizeDeviceStorage() const
    {
        // update the device storage with the host and free the host copy
        if(pHostStorage_)
        {
            pDeviceStorage_->CopyHostToDeviceStorage(*pHostStorage_);
            clearHostStorage();
        }
    }

    /**
    * @return reference to the device storage
    */
    __CPU__
    const DeviceStorageType& deviceStorage() const { return *pDeviceStorage_; }

    /**
    * @return true if device buffer is up to date
    */
    __CPU__
    bool isDeviceUpToDate() const { return !pHostStorage_; }

private:
    // method initializing host storage
    void initHostStorage() const
    {
        // make sure the implementation is called via DataField
        assert(isInitialized());

        // clean up an existing host storage if necessary
        clearHostStorage();

        // allocate the new host storage
        pHostStorage_ = new HostStorageType();

        // initialize the new host storage
        KBoundary kBoundary;
        kBoundary.Init(boundary().kMinusOffset(), boundary().kPlusOffset());
        pHostStorage_->Init(calculationDomain(), kBoundary);
    }

    // method cleaning up host storage if necessary
    void clearHostStorage() const
    {
        if(pHostStorage_) 
        {
            delete pHostStorage_;
            pHostStorage_ = NULL;
        }    
    }

    DeviceStorageType* pDeviceStorage_;
    mutable HostStorageType* pHostStorage_; // only allocated if the valid copy is on the host!
};

// add data field trait
template<
    typename TValue,
    typename TStorageFormat>
struct is_data_field<DataFieldCUDA<TValue, TStorageFormat> > : boost::mpl::true_ {};
  
/**
* @struct create_storage_pointer
* Meta function creating a storage pointer given a data field
*/
template<
    typename TDataField,
    typename TBlockSize,
    ParameterIntend VParameterIntend,
    int VParameterIndex>
struct create_storage_pointer;

template<
    typename TValue,
    typename TStorageFormat,
    typename TBlockSize,
    int VParameterIndex>
struct create_storage_pointer<DataFieldCUDA<TValue, TStorageFormat>, TBlockSize, cIn, VParameterIndex>
{
    typedef DataFieldCUDAStoragePointer<ScalarPointer<TValue, TStorageFormat, TBlockSize, cReadOnly> > type;
};

template<
    typename TValue,
    typename TStorageFormat,
    typename TBlockSize,
    int VParameterIndex>
struct create_storage_pointer<DataFieldCUDA<TValue, TStorageFormat>, TBlockSize, cInOut, VParameterIndex>
{
    typedef DataFieldCUDAStoragePointer<ScalarPointer<TValue, TStorageFormat, TBlockSize, cReadWrite> > type;
};

