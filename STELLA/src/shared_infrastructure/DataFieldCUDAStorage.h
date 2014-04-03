#pragma once

#include <cassert>
#include "CUDADefinitions.h"
#include "BlockSize.h"
#include "DataFieldStorageFormat.h"
#include "DataFieldOpenMPStorage.h"
#include "DataFieldCUDAStoragePointer.h"

/**
* @class DataFieldCUDAStorage
* Storage class which abstracts the device storage of a DataField from the host perspective.
*/
template<
    typename TValue,
    typename TStorageFormat>
class DataFieldCUDAStorage : public DataFieldStorage<TValue, TStorageFormat, DataFieldCUDAStorage<TValue, TStorageFormat> > // CRTP
{
    BOOST_STATIC_ASSERT(is_data_field_storage_format<TStorageFormat>::value);
public:
    using DataFieldStorage<TValue, TStorageFormat, DataFieldCUDAStorage>::Assign;
    using DataFieldStorage<TValue, TStorageFormat, DataFieldCUDAStorage>::paddedSize;
    using DataFieldStorage<TValue, TStorageFormat, DataFieldCUDAStorage>::originOffset;

    typedef typename DataFieldStorage<TValue, TStorageFormat, DataFieldCUDAStorage>::ValueType ValueType;
    typedef typename DataFieldStorage<TValue, TStorageFormat, DataFieldCUDAStorage>::PointerType PointerType;
    typedef typename DataFieldStorage<TValue, TStorageFormat, DataFieldCUDAStorage>::ExternalStorageType ExternalStorageType;
 
    __CPU__
    DataFieldCUDAStorage() 
    { 
        pStorageMemory_ = NULL;
        pStorageBase_ = NULL;
        pStorageOrigin_ = NULL;
        storageSize_ = 0;
        initializedExternally_ = false;
    }
    __CPU__
    ~DataFieldCUDAStorage() 
    {
        clearStorageMemory();
    }

    __CPU__
    DataFieldCUDAStorage(const DataFieldCUDAStorage& other) 
    { 
        // make sure pStorageMemory is initialized to NULL 
        // otherwise the assignment operator will try the free the memory
        pStorageMemory_ = NULL; 
        *this = other; 
    }
    __CPU__
    DataFieldCUDAStorage& operator= (const DataFieldCUDAStorage& other)
    {
        // run the base assign
        Assign(other);
        storageSize_ = other.storageSize_;
        storageStrides_ = other.storageStrides_;
        initializedExternally_ = false;

        // finish initialization in case the other field was initialized
        if(storageSize_ != 0)
        {
            // provide storage
            initStorageMemory();

            // copy the data from the other storage
            cudaMemcpy(pStorageBase_, other.pStorageBase_, sizeof(ValueType) * storageSize_, cudaMemcpyDeviceToDevice);
            assertNoCUDAError("DataFieldCUDAStorage");
        }
        else
        {
            clearStorageMemory();
        }

        // by convention
        return *this;
    }

    /**
    * Init the CUDA storage
    */
    __CPU__
    void InitImpl() 
    { 
        // init storage size
        initStorageSizeAndStrides();

        // init the storage memory
        initStorageMemory(); 
    }

    /**
    * Init the CUDA storage
    * @param externalStorage external storage
    */
    __CPU__
    void InitImpl(const ExternalStorageType& externalStorage)
    {
        // init the storage size
        initStorageSizeAndStrides();

        // init the storage memory using an external storage
        initStorageMemory(externalStorage);
    }

    /**
    * Set CUDA storage to external storage
    * @param externalStorage external storage
    */
    __CPU__
    void SetExternalStorageImpl(const ExternalStorageType& externalStorage)
    {
        // init the storage memory using an external storage 
        // without previously recomputing the storage size
        initStorageMemory(externalStorage);
    }

    /**
    * @return pointer to the base of the storage
    */
    __CPU__
    PointerType pStorageBaseImpl() const { return pStorageBase_; }

    /**
    * Method copying the content of the device storage to the host storage
    */
    __CPU__
    void CopyDeviceToHostStorage(const DataFieldOpenMPStorage<TValue, TStorageFormat>& hostStorage) const
    {
        // get the pointer to the base of the host storage
        PointerType pHostStorageBase = &hostStorage.originIterator().At(
            -originOffset().iIndex(),
            -originOffset().jIndex(),
            -originOffset().kIndex()
        );

        // copy device storage to host storage
        cudaMemcpy(pHostStorageBase, pStorageBase_, sizeof(ValueType) * storageSize_, cudaMemcpyDeviceToHost);
        assertNoCUDAError("DataFieldCUDAStorage");
    }

    /**
    * Method copying the content of the host storage to the device storage
    */
    __CPU__
    void CopyHostToDeviceStorage(const DataFieldOpenMPStorage<TValue, TStorageFormat>& hostStorage) const
    {
        // get the pointer to the base of the host storage
        PointerType pHostStorageBase = &hostStorage.originIterator().At(
            -originOffset().iIndex(),
            -originOffset().jIndex(),
            -originOffset().kIndex()
        );

        // copy device storage to host storage
        cudaMemcpy(pStorageBase_, pHostStorageBase, sizeof(ValueType) * storageSize_, cudaMemcpyHostToDevice);
        assertNoCUDAError("DataFieldCUDAStorage");
    }

    /**
    * Method initializing a storage pointer to the origin position
    * @param storagePointer pointer to initialize
    */
    template<typename TStoragePointer>
    __CPU__
    void InitializeStoragePointerToOrigin(DataFieldCUDAStoragePointer<TStoragePointer>& storagePointer) const 
    {
        storagePointer.Init(pStorageOrigin_);
    }

    /**
    * Method initializing CUDA storage strides. 
    * In case the strides where already initialized for a different size the application terminates with an error.
    * (the method is used to setup the strides for CUDA kernels where there is only one stride for all fields of the same type.
    * @param storageStrides strides to initialize
    */
    template<typename TStorageStrides>
    __CPU__
    void TryToInitializeStorageStrides(DataFieldCUDAStorageStrides<TStorageStrides>& storageStrides) const 
    {
        // if the strides are uninitialized -> initialize them
        if( storageStrides.strides().ComputeStride(1, 0, 0) == -1 ||
            storageStrides.strides().ComputeStride(0, 1, 0) == -1 ||
            storageStrides.strides().ComputeStride(0, 0, 1) == -1 )
        {
            storageStrides.Init(storageStrides_);
        }
        // otherwise check the strides are identical
        else if(storageStrides.strides().ComputeStride(1, 0, 0) != storageStrides_.ComputeStride(1, 0, 0) ||
                storageStrides.strides().ComputeStride(0, 1, 0) != storageStrides_.ComputeStride(0, 1, 0) ||
                storageStrides.strides().ComputeStride(0, 0, 1) != storageStrides_.ComputeStride(0, 0, 1) )
        {
            std::cerr << "DataFieldCUDAStorage -> failed to initialize storage strides, already set to different value" << std::endl;
            assert(false);
            exit(-1);
        }
    }

    /**
    * @return device pointer to the origin of the storage
    */
    __CPU__
    PointerType pStorageOrigin() const { return pStorageOrigin_; };

private:
    // method freeing private memory of the storage
    __CPU__
    void clearStorageMemory()
    {
        if(pStorageMemory_) 
        {
            if(!initializedExternally_)
            {
                cudaFree(pStorageMemory_);
                assertNoCUDAError("DataFieldCUDAStorage");
            }
            pStorageMemory_ = NULL;
            pStorageBase_ = NULL;
            pStorageOrigin_ = NULL;
        }
    }

    // compute the strides of the origin
    __CPU__
    int computeOriginStride()
    {
        return storageStrides_.ComputeStride(originOffset().iIndex(), originOffset().jIndex(), originOffset().kIndex());
    }

    // setup the size information return the origin stride
    __CPU__
    void initStorageSizeAndStrides()
    {
        // allocate the memory (first calculate the size using the storage size functions)
        storageSize_ = paddedSize().iSize() * paddedSize().jSize() * paddedSize().kSize();
        assert(storageSize_ != 0);

        // init the strides
        storageStrides_.Init(paddedSize());
    }
    
    // method setting up storage memory
    __CPU__
    void initStorageMemory()
    {
        // free memory if it was allocated 
        clearStorageMemory();
        
        // verify the storageSize_ variable is up to date
        assert(storageSize_ == paddedSize().iSize() * paddedSize().jSize() * paddedSize().kSize());

        // allocate the memory
        cudaMalloc(reinterpret_cast<void**>(&pStorageMemory_), sizeof(ValueType) * (storageSize_ + storage_alignment<TStorageFormat>::value));
        assertNoCUDAError("DataFieldCUDAStorage");

        // setup pointers
        int originStride = computeOriginStride();
        pStorageOrigin_ = aligned_pointer<TStorageFormat>(pStorageMemory_ + originStride);
        pStorageBase_ = pStorageOrigin_ - originStride;
    }

    // method initializing the storage memory using external storage
    __CPU__
    void initStorageMemory(const ExternalStorageType& externalStorage)
    {
        // verify the external storage size
        assert(paddedSize().iSize() == externalStorage.paddedSize().iSize());
        assert(paddedSize().jSize() == externalStorage.paddedSize().jSize());
        assert(paddedSize().kSize() == externalStorage.paddedSize().kSize());

        // verify the storageSize_ variable is up to date
        assert(storageSize_ == paddedSize().iSize() * paddedSize().jSize() * paddedSize().kSize());

        // clear storage memory if external memory is available
        clearStorageMemory();

        // setup the pointers
        int originStride = computeOriginStride();
        pStorageBase_ = externalStorage.pStoragePointer();
        pStorageOrigin_ = pStorageBase_ + originStride;

        // make sure the external memory is aligned
        assert(pStorageOrigin_ == aligned_pointer<TStorageFormat>(pStorageOrigin_));
        initializedExternally_ = true;
    }

    bool initializedExternally_;
    int storageSize_;
    PointerType pStorageMemory_; // pointer to the start of the memory block (unaligned)
    PointerType pStorageBase_; // pointer to the actual start of the array
    PointerType pStorageOrigin_; // pointer to the origin

    DataFieldStorageStrides<typename TStorageFormat::StorageOrder> storageStrides_; // stride member variable
};

  
