#pragma once

#include <boost/mpl/if.hpp>
#include <boost/mpl/bool.hpp>
#include "Enums.h"
#include "ParameterTraits.h"
#include "BlockStorage.h"
#include "DataFieldCUDAStorageIndex.h"
#include "DataFieldCUDAStorageStrides.h"

/**
* @class ScalarPointer
* Used by DataFieldCUDAStoragePointer class to implement a pointer to a scalar element
*/
template<
    typename TValue,
    typename TStorageFormat,
    typename TBlockSize,
    AccessPolicy VAccessPolicy>
class ScalarPointer
{
public:
    // define the value types
    typedef TValue ValueType;
    typedef ValueType& ReferenceType;
    typedef typename boost::mpl::if_c<
        VAccessPolicy == cReadWrite, 
        ReferenceType,
        ValueType
    >::type ReturnType;
    typedef TValue* PointerType;

    // define the storage stride type
    typedef DataFieldCUDAStorageStrides<ScalarStorageStrides<TStorageFormat, TBlockSize> > StorageStridesType;
    
    __ACC_CPU__
    ScalarPointer() { pValue_ = NULL; }
    __ACC_CPU__
    ~ScalarPointer() {}

    __ACC_CPU__
    ScalarPointer(const ScalarPointer& other) { *this = other; }
    __ACC_CPU__
    ScalarPointer& operator= (const ScalarPointer& other) 
    { 
        pValue_ = other.pValue_;
        // by convention
        return *this;
    }

    /**
    * Initialize the pointer
    */
    __CPU__
    void Init(const PointerType pValue) { pValue_ = pValue; }
    
    /**
    * @return true if pointers match
    */
    __CPU__
    bool PointsTo(const PointerType pValue) const { return pValue_ == pValue; }

    /** 
    * Access operator
    * @param index index of the access
    * @return value at the index position
    */
    __ACC__
    ReturnType operator[] (const int index) const 
    {
        return accessImpl(index, static_cast<boost::mpl::integral_c<AccessPolicy, VAccessPolicy>*>(0));
    }

    /** 
    * Access method which bypasses the texture cache
    * (used in order to implement software managed caching which is not polluting the hardware caches)
    * @param index index of the access
    * @return value at the index position
    */
    __ACC__
    ReturnType BypassCache(const int index) const 
    {
        // use the read write access policy which doesn't use the texture memory
        return accessImpl(index, static_cast<boost::mpl::integral_c<AccessPolicy, cReadWrite>*>(0));
    }


private:
    // provide separate access functions for read write and read only pointers
    __ACC__
    ReturnType accessImpl(const int index, boost::mpl::integral_c<AccessPolicy, cReadWrite>*) const
    {
        return pValue_[index];
    }
    __ACC__
    ReturnType accessImpl(const int index, boost::mpl::integral_c<AccessPolicy, cReadOnly>*) const
    {
#if __CUDA_ARCH__ >= 350
        // on Kepler use ldg to read directly via read only cache
        return __ldg(pValue_ + index);
#else
        return pValue_[index];
#endif
    }

    PointerType pValue_;
};

/**
* @class BlockStoragePointer
* Used by DataFieldCUDAStoragePointer class to implement a pointer to a block storage element
*/
template<
    typename TValue,
    typename TStorageFormat,
    typename TBlockStorageFormat,
    AccessPolicy VAccessPolicy>
class BlockStoragePointer
{
public:
    // define the value types
    typedef TValue ValueType;
    typedef ValueType& ReferenceType;
    typedef typename boost::add_const<ValueType>::type& ConstReferenceType;
    typedef typename boost::mpl::if_c<
        VAccessPolicy == cReadWrite, 
        ReferenceType,
        ValueType
    >::type ReturnType;
    typedef BlockStorage<ValueType, TBlockStorageFormat>* PointerType;
    typedef ValueType* BlockElementPointerType;

    // define the storage stride type
    typedef DataFieldCUDAStorageStrides<BlockStorageStrides<TStorageFormat, TBlockStorageFormat> > StorageStridesType;

    __ACC_CPU__
    BlockStoragePointer() { pValue_ = NULL; }
    __ACC_CPU__
    ~BlockStoragePointer() {}

    __ACC_CPU__
    BlockStoragePointer(const BlockStoragePointer& other) { *this = other; }
    __ACC_CPU__
    BlockStoragePointer& operator= (const BlockStoragePointer& other) 
    { 
        pValue_ = other.pValue_;
        // by convention
        return *this;
    }

    __CPU__
    void Init(const PointerType pBlockStorage) 
    { 
        pValue_ = reinterpret_cast<BlockElementPointerType>(pBlockStorage) + TBlockStorageFormat::OriginOffset::value; // add the origin offset to the pointer
    }

    /**
    * @return true if pointers match
    */
    __CPU__
    bool PointsTo(const PointerType pBlockStorage) const 
    { 
        return pValue_ == reinterpret_cast<BlockElementPointerType>(pBlockStorage) + TBlockStorageFormat::OriginOffset::value; 
    }

    /** 
    * Access operator
    * @param index index of the access
    * @return value at the index position
    */
    __ACC__
    ReturnType operator[] (const int index) const 
    { 
        return accessImpl(index, static_cast<boost::mpl::integral_c<AccessPolicy, VAccessPolicy>*>(0));
    }

    /** 
    * Access method bypassing the texture cache
    * (used in order to implement software managed caching which is not polluting the hardware caches)
    * @param index index of the access
    * @return value at the index position
    */
    __ACC__
    ReturnType BypassCache(const int index) const 
    {
        // use the read write access policy which bypasses the texture memory
        return accessImpl(index, static_cast<boost::mpl::integral_c<AccessPolicy, cReadWrite>*>(0));
    }

private:
    // provide separate access functions for read write and read only pointers
    __ACC__
    ReturnType accessImpl(const int index, boost::mpl::integral_c<AccessPolicy, cReadWrite>*) const
    {
        return pValue_[index];
    }
    __ACC__
    ReturnType accessImpl(const int index, boost::mpl::integral_c<AccessPolicy, cReadOnly>*) const
    {
#if __CUDA_ARCH__ >= 350
        // on Kepler use ldg to read directly via read only cache
        return __ldg(pValue_ + index);
#else
        return pValue_[index];
#endif
    }

    BlockElementPointerType pValue_;
};

/**
* @class DataFieldCUDAStoragePointer
* Class pointing to the device storage of a CUDA data field.
* The class provides the functionality of an iterator when it is combined with an index and a stride class.
*/
template<typename TStoragePointer>
class DataFieldCUDAStoragePointer
{
public:    
    // define value types
    typedef typename TStoragePointer::ValueType ValueType;
    typedef typename TStoragePointer::PointerType PointerType;
    typedef typename TStoragePointer::ReturnType ReturnType;

    // define strides and index types
    typedef typename TStoragePointer::StorageStridesType StorageStridesType;
    typedef DataFieldCUDAStorageIndex<StorageStridesType> StorageIndexType;
    
    __ACC_CPU__
    DataFieldCUDAStoragePointer() {}
    __ACC_CPU__
    ~DataFieldCUDAStoragePointer() {}

    __ACC_CPU__
    DataFieldCUDAStoragePointer(const DataFieldCUDAStoragePointer& other) { *this = other; }
    __ACC_CPU__
    DataFieldCUDAStoragePointer& operator= (const DataFieldCUDAStoragePointer& other) 
    { 
        pointer_ = other.pointer_;
        // by convention
        return *this;
    }
    
    /**
    * Method initializing the pointer to the origin of the device storage.
    * @param pOrigin pointer to the origin of the device storage
    */
    __CPU__
    void Init(const PointerType pOrigin)
    {
        pointer_.Init(pOrigin);
    }

    /**
    * Method used to verify if a storage pointer points to the right origin or if an update is necessary
    * @param pOrigin origin pointer to compare with
    * @return true if the storage pointer points the pOrigin
    */
    __CPU__
    bool PointsToOrigin(const PointerType pOrigin) const 
    {
        return pointer_.PointsTo(pOrigin);
    }

    /**
    * Access method dereferencing the pointer at the center position given the storage index
    * @param index data field storage index pointing to the access position
    * @return reference to the value at the offset position
    */
    __ACC__
    ReturnType Center(const StorageIndexType& index) const
    {
        return pointer_[index.index()];
    }

    /**
    * Access method dereferencing the pointer at an offset given the storage strides and index
    * @param index data field storage index pointing to the access position
    * @param storageStrides storageStrides of the underlying data field
    * @param i offset in i direction
    * @param j offset in j direction
    * @param k offset in k direction
    * @return reference to the value at the offset position
    */
    __ACC__
    ReturnType At(const StorageIndexType& index, const StorageStridesType& storageStrides, const int i, const int j, const int k) const
    {
        return pointer_[index.index() + storageStrides.ComputeStride(i, j, k)];
    }

    /**
    * Access method dereferencing the pointer at an offset given the storage strides and index
    * (note that the method bypasses the texture memory of the GPU otherwise its equivalent to the At method)
    * @param index data field storage index pointing to the access position
    * @param storageStrides storageStrides of the underlying data field
    * @param i offset in i direction
    * @param j offset in j direction
    * @param k offset in k direction
    * @return reference to the value at the offset position
    */
    __ACC__
    ReturnType BypassCache(const StorageIndexType& index, const StorageStridesType& storageStrides, const int i, const int j, const int k) const
    {
        return pointer_.BypassCache(index.index() + storageStrides.ComputeStride(i, j, k));
    }

private:
    TStoragePointer pointer_;
};

// return type specialization
template<typename TStoragePointer>
struct return_type<DataFieldCUDAStoragePointer<TStoragePointer> >
{
    typedef typename DataFieldCUDAStoragePointer<TStoragePointer>::ReturnType type;
};

// value type specialization
template<typename TStoragePointer>
struct value_type<DataFieldCUDAStoragePointer<TStoragePointer> >
{
    typedef typename DataFieldCUDAStoragePointer<TStoragePointer>::ValueType type;
};

// is iterable specialization
template<typename TStoragePointer>
struct is_iterable<DataFieldCUDAStoragePointer<TStoragePointer> > : boost::mpl::true_ {};

/**
* @struct create_storage_strides
* Meta function returning the storage strides associated to a CUDA storage pointer
*/
template<typename T>
struct create_storage_strides; 

template<typename TStoragePointer>
struct create_storage_strides<DataFieldCUDAStoragePointer<TStoragePointer> > 
{
    typedef typename DataFieldCUDAStoragePointer<TStoragePointer>::StorageStridesType type; 
};


