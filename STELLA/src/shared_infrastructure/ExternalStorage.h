#pragma once

#include "Definitions.h"
#include "IJKSize.h"

/**
* @class ExternalStorage 
* Class holding a pointer and the size of external storage which can be used to initialize a data field
*/
template<typename TValue>
class ExternalStorage
{
public:
    typedef TValue* PointerType;

    __CPU__
    ExternalStorage() 
    { 
        pStoragePointer_ = NULL;
        paddedSize_.Init(0, 0, 0);
    }
    __CPU__
    ~ExternalStorage() {}

    __CPU__
    ExternalStorage(const ExternalStorage& other) { *this = other; }
    __CPU__
    ExternalStorage& operator= (const ExternalStorage& other)
    {
        pStoragePointer_ = other.pStoragePointer_;
        paddedSize_ = other.paddedSize_;
        // by convention
        return *this;
    }

    /**
    * Method initializing the external storage
    * @param pStoragePointer pointer to the beginning of the storage
    * @param paddedSize total storage size
    */
    __CPU__
    void Init(PointerType pStoragePointer, const IJKSize& paddedSize)
    {
        pStoragePointer_ = pStoragePointer;
        paddedSize_ = paddedSize;
    }

    /**
    * @return pointer to the beginning of the storage
    */
    __CPU__
    PointerType pStoragePointer() const { return pStoragePointer_; }
    
    /**
    * @return total storage size
    */
    __CPU__
    const IJKSize& paddedSize() const { return paddedSize_; }

private:
    PointerType pStoragePointer_;
    IJKSize paddedSize_;
};

