#pragma once

#include "Definitions.h"

/**
* @class IJBoundary
* Container for boundary offsets in i and j direction. 
* Positive numbers indicate count arrow direction, negative numbers the opposite.
* Therefore i-/j-minus are typically negative.
* 
*            ^ j-plus
*            |
* i-minus <-----> i-plus
*            |
*            v j-minus
*/
class IJBoundary
{
public:
    __ACC_CPU__
    IJBoundary()
    {
        iMinusOffset_ = 0;
        iPlusOffset_  = 0;
        jMinusOffset_ = 0;
        jPlusOffset_  = 0;
    }
    __ACC_CPU__
    ~IJBoundary() {}

    __ACC_CPU__
    IJBoundary(const IJBoundary& other) { *this = other; }
    __ACC_CPU__
    IJBoundary& operator= (const IJBoundary& other)
    {
        iMinusOffset_ = other.iMinusOffset_;
        iPlusOffset_  = other.iPlusOffset_;
        jMinusOffset_ = other.jMinusOffset_;
        jPlusOffset_  = other.jPlusOffset_;
        // by convention
        return *this;
    }
    __ACC_CPU__
    bool operator== (const IJBoundary& other) const
    {
        return (
            (iMinusOffset() == other.iMinusOffset()) &&
            (iPlusOffset()  == other.iPlusOffset()) &&
            (jMinusOffset() == other.jMinusOffset()) &&
            (jPlusOffset()  == other.jPlusOffset()) 
        );
    }

    /**
    * Init the container
    * @param iMinusOffset Offset of boundary in i-minus direction of domain
    * @param iPlusOffset Offset of boundary in i-plus direction of domain
    * @param jMinusOffset Offset of boundary in j-minus direction of domain
    * @param jPlusOffset Offset of boundary in j-plus direction of domain
    */
    __ACC_CPU__
    void Init(const int iMinusOffset, const int iPlusOffset, const int jMinusOffset, const int jPlusOffset)
    {
        // store the size per dimension
        iMinusOffset_ = iMinusOffset;
        iPlusOffset_  = iPlusOffset;
        jMinusOffset_ = jMinusOffset;
        jPlusOffset_  = jPlusOffset;
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

private:
    int iMinusOffset_;
    int iPlusOffset_;
    int jMinusOffset_;
    int jPlusOffset_;
};

  
