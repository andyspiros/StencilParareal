#pragma once

#include <cassert>
#include <algorithm>
#include <string>
#include <sstream>
#include "Definitions.h"

/**
* @class IJKIndex
* Container for a i, j, k Index
*/
class IJKIndex
{
public:
    __ACC_CPU__
    IJKIndex()
    {
        iIndex_ = 0;
        jIndex_ = 0;
        kIndex_ = 0;
    }
    __ACC_CPU__
    ~IJKIndex() {}

    __ACC_CPU__
    IJKIndex(const IJKIndex& other)
    {
        *this = other;
    }
    __ACC_CPU__
    IJKIndex& operator= (const IJKIndex& other)
    {
        iIndex_ = other.iIndex_;
        jIndex_ = other.jIndex_;
        kIndex_ = other.kIndex_; 
        // by convention
        return *this;
    }
    __ACC_CPU__
    bool operator== (const IJKIndex& other) const
    {
        return (
            (iIndex() == other.iIndex()) &&
            (jIndex() == other.jIndex()) &&
            (kIndex() == other.kIndex())
        );
    }

    /**
    * Init the container
    * @param iIndex Index in i dimension
    * @param jIndex Index in j dimension
    * @param kIndex Index in k dimension
    */
    __ACC_CPU__
    void Init(const int iIndex, const int jIndex, const int kIndex)
    {
        // store the Index per dimension
        iIndex_ = iIndex;
        jIndex_ = jIndex;
        kIndex_ = kIndex;
    }

    /**
    * Add offset to index
    * @param i Offset in i dimension
    * @param j Offset in j dimension
    * @param k Offset in k dimension
    */
    __ACC_CPU__
    void AddOffset(const int i, const int j, const int k)
    {
        // store the Index per dimension
        iIndex_ += i;
        jIndex_ += j;
        kIndex_ += k;
    }

    /**
    * CommonMin sets indexes to minimum of this and other indexes
    */
    __CPU__
    void CommonMin(const IJKIndex& other)
    {
        iIndex_ = std::min(iIndex_, other.iIndex_);
        jIndex_ = std::min(jIndex_, other.jIndex_);
        kIndex_ = std::min(kIndex_, other.kIndex_);
    }

    /**
    * CommonMax sets indexes to maximum of this and other indexes
    */
    __CPU__
    void CommonMax(const IJKIndex& other)
    {
        iIndex_ = std::max(iIndex_, other.iIndex_);
        jIndex_ = std::max(jIndex_, other.jIndex_);
        kIndex_ = std::max(kIndex_, other.kIndex_);
    }

    /**
    * Index in i dimension
    * @return Index
    */
    __ACC_CPU__
    int iIndex() const { return iIndex_; }

    /**
    * Index in j dimension
    * @return Index
    */
    __ACC_CPU__
    int jIndex() const { return jIndex_; }
    
    /**
    * Index in k dimension
    * @return Index
    */
    __ACC_CPU__
    int kIndex() const { return kIndex_; }

    /**
    * @return the index as a string
    */
    __CPU__
    std::string ToString() const
    {
        std::stringstream stream;
        stream << iIndex_ << "\t" 
            << jIndex_ << "\t" 
            << kIndex_;
        return stream.str();
    }

private:
    int iIndex_;
    int jIndex_;
    int kIndex_;
};

  
