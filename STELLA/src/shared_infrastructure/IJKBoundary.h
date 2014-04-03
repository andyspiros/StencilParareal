#pragma once

#include "Definitions.h"

/**
* @class IJKBoundary
* Container for boundary offsets in all directions of 3D domain
* Positive numbers indicate count arrow direction, negative numbers the opposite.
* Therefore i-/j-/k-minus are typically negative.
* 
*            ^ j-plus
*            |
* i-minus <-----> i-plus
*            |
*            v j-minus
*/
class IJKBoundary
{
public:
    __ACC_CPU__
    IJKBoundary()
    {
        iMinusOffset_ = 0;
        iPlusOffset_  = 0;
        jMinusOffset_ = 0;
        jPlusOffset_  = 0;
        kMinusOffset_ = 0;
        kPlusOffset_  = 0;
    }
    __ACC_CPU__
    ~IJKBoundary() {}

    __ACC_CPU__
    IJKBoundary(const IJKBoundary& other) { *this = other; }
    __ACC_CPU__
    IJKBoundary& operator= (const IJKBoundary& other)
    {
        iMinusOffset_ = other.iMinusOffset_;
        iPlusOffset_  = other.iPlusOffset_;
        jMinusOffset_ = other.jMinusOffset_;
        jPlusOffset_  = other.jPlusOffset_;
        kMinusOffset_ = other.kMinusOffset_;
        kPlusOffset_  = other.kPlusOffset_;
        // by convention
        return *this;
    }
    __ACC_CPU__
    bool operator== (const IJKBoundary& other) const
    {
        return (
            (iMinusOffset() == other.iMinusOffset()) &&
            (iPlusOffset()  == other.iPlusOffset()) &&
            (jMinusOffset() == other.jMinusOffset()) &&
            (jPlusOffset()  == other.jPlusOffset()) &&
            (kMinusOffset() == other.kMinusOffset()) &&
            (kPlusOffset()  == other.kPlusOffset())
        );
    }

    /**
    * Init the container
    * @param iMinusOffset Offset of boundary in i-minus direction of domain
    * @param iPlusOffset Offset of boundary in i-plus direction of domain
    * @param jMinusOffset Offset of boundary in j-minus direction of domain
    * @param jPlusOffset Offset of boundary in j-plus direction of domain
    * @param kMinusOffset Offset of boundary in k-minus direction of domain
    * @param kPlusOffset Offset of boundary in k-plus direction of domain
    */
    __ACC_CPU__
    void Init(const int iMinusOffset, const int iPlusOffset, const int jMinusOffset, const int jPlusOffset, const int kMinusOffset, const int kPlusOffset)
    {
        // store the size per dimension
        iMinusOffset_ = iMinusOffset;
        iPlusOffset_  = iPlusOffset;
        jMinusOffset_ = jMinusOffset;
        jPlusOffset_  = jPlusOffset;
        kMinusOffset_ = kMinusOffset;
        kPlusOffset_  = kPlusOffset;
    }

    /**
    * @return Offset in i-minus direction
    */
    __ACC_CPU__
    int iMinusOffset() const { return iMinusOffset_; }

    /**
    * @return Offset in i-plus direction
    */
    __ACC_CPU__
    int iPlusOffset() const { return iPlusOffset_; }

    /**
    * @return Offset in j-minus direction
    */
    __ACC_CPU__
    int jMinusOffset() const { return jMinusOffset_; }

    /**
    * @return Offset in j-plus direction
    */
    __ACC_CPU__
    int jPlusOffset() const { return jPlusOffset_; }

    /**
    * @return Offset in k-minus direction
    */
    __ACC_CPU__
    int kMinusOffset() const { return kMinusOffset_; }

    /**
    * @return Offset in k-plus direction
    */
    __ACC_CPU__
    int kPlusOffset() const { return kPlusOffset_; }

private:
    int iMinusOffset_;
    int iPlusOffset_;
    int jMinusOffset_;
    int jPlusOffset_;
    int kMinusOffset_;
    int kPlusOffset_;
};

  
