#pragma once

#include <cassert>
#include "Definitions.h"

/**
* @class IJKSize
* Container for i, j, k Size
*/
class IJKSize
{
public:
    __ACC_CPU__
    IJKSize()
    {
        iSize_ = 0;
        jSize_ = 0;
        kSize_ = 0;
    }
    __ACC_CPU__
    ~IJKSize() {}

    __ACC_CPU__
    IJKSize(const IJKSize& other)
    {
        *this = other;
    }
    __ACC_CPU__
    IJKSize& operator= (const IJKSize& other)
    {
        iSize_ = other.iSize_;
        jSize_ = other.jSize_;
        kSize_ = other.kSize_; 
        // by convention
        return *this;
    }
    __ACC_CPU__
    bool operator== (const IJKSize& other) const
    {
        return (
            (iSize() == other.iSize()) && 
            (jSize() == other.jSize()) && 
            (kSize() == other.kSize())
        );
    }

    /**
    * Init the container
    * @param iSize size in i dimension
    * @param jSize size in j dimension
    * @param kSize size in k dimension
    */
    __ACC_CPU__
    void Init(const int iSize, const int jSize, const int kSize)
    {
        // assert that sizes are not negative
        ACC_ASSERT(iSize >= 0 && jSize >= 0 && kSize >= 0);

        // store the size per dimension
        iSize_ = iSize;
        jSize_ = jSize;
        kSize_ = kSize;
    }

    /**
    * Size in i dimension
    * @return size
    */
    __ACC_CPU__
    int iSize() const { return iSize_; }

    /**
    * Size in j dimension
    * @return size
    */
    __ACC_CPU__
    int jSize() const { return jSize_; }
    
    /**
    * Size in k dimension
    * @return size
    */
    __ACC_CPU__
    int kSize() const { return kSize_; }

    /**
    * @return true if the size is empty
    */
    __ACC_CPU__
    bool empty() const { return iSize_ <= 0 || jSize_ <=0 || kSize_ <=0; }

private:
    int iSize_;
    int jSize_;
    int kSize_;
};
