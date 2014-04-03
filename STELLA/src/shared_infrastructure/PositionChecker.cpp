#include <cassert>
#include <iostream>
#include <algorithm>
#include <limits>
#include "PositionChecker.h"

PositionChecker::PositionChecker() :
    checkI_(false), checkJ_(false), checkK_(false)
{
    currentPosition_.Init(0, 0, 0);
    ResetBoundingBox();
}

PositionChecker::~PositionChecker() {}

PositionChecker::PositionChecker(const PositionChecker& other)
{
    currentPosition_.Init(0, 0, 0);
    *this = other;
}

PositionChecker& PositionChecker::operator= (const PositionChecker& other)
{
    SetValidRange(other.validRangeFloor_, other.validRangeCeiling_);
    currentPosition_ = other.currentPosition_;
    checkI_  = other.checkI_;
    checkJ_  = other.checkJ_;
    checkK_  = other.checkK_;
    boundingBoxFloor_ = other.boundingBoxFloor_;
    boundingBoxCeiling_ = other.boundingBoxCeiling_;

    // by convention
    return *this;
}

void PositionChecker::Init(const IJKIndex& floor, const IJKIndex& ceiling, const bool checkI, const bool checkJ, const bool checkK)
{
    SetValidRange(floor, ceiling);
    ResetBoundingBox();
    currentPosition_.Init(0, 0, 0);
    checkI_ = checkI;
    checkJ_ = checkJ;
    checkK_ = checkK;
}

const IJKIndex& PositionChecker::currentPosition() const
{
    return currentPosition_;
}

void PositionChecker::SetValidRange(const IJKIndex& floor, const IJKIndex& ceiling)
{
    validRangeFloor_ = floor;
    validRangeCeiling_ = ceiling;
}

const IJKIndex& PositionChecker::validRangeFloor() const
{
    return validRangeFloor_;
}

const IJKIndex& PositionChecker::validRangeCeiling() const
{
    return validRangeCeiling_;
}


void PositionChecker::ResetBoundingBox()
{
    boundingBoxFloor_.Init(
        std::numeric_limits<int>::max(),
        std::numeric_limits<int>::max(),
        std::numeric_limits<int>::max()
    );
    boundingBoxCeiling_.Init(
        std::numeric_limits<int>::min(),
        std::numeric_limits<int>::min(),
        std::numeric_limits<int>::min()
    );
}

void PositionChecker::MergeBoundingBox(const PositionChecker& other)
{
    boundingBoxFloor_.CommonMin(other.boundingBoxFloor());
    boundingBoxCeiling_.CommonMax(other.boundingBoxCeiling());
}

const IJKIndex& PositionChecker::boundingBoxFloor() const 
{
    return boundingBoxFloor_;
}

const IJKIndex& PositionChecker::boundingBoxCeiling() const 
{
    return boundingBoxCeiling_;
}

void PositionChecker::SetPosition(const PositionChecker& other)
{
    assert(validRangeFloor_ == other.validRangeFloor_);
    assert(validRangeCeiling_ == other.validRangeCeiling_);
    currentPosition_ = other.currentPosition_;
}

void PositionChecker::Advance(const int i, const int j, const int k)
{
    currentPosition_.AddOffset(i, j, k);
}


void PositionChecker::Advance(Dimension dim, const int d)
{
    switch(dim)
    {
    case cDimI:
        Advance(d, 0, 0);
        break;
    case cDimJ:
        Advance(0, d, 0);
        break;
    case cDimK:
        Advance(0, 0, d);
        break;
    }
}

void PositionChecker::CheckPositionAtOffset(const int i, const int j, const int k)
{
    // initialize in range flags to true
    bool isIPositionOutOfRange = false;
    bool isJPositionOutOfRange = false;
    bool isKPositionOutOfRange = false;

    // check the offset
    if(checkI_) 
    {
        isIPositionOutOfRange = 
            currentPosition_.iIndex() + i < validRangeFloor_.iIndex() || 
            currentPosition_.iIndex() + i > validRangeCeiling_.iIndex();
    }
    if(checkJ_)
    {
        isJPositionOutOfRange = 
            currentPosition_.jIndex() + j < validRangeFloor_.jIndex() || 
            currentPosition_.jIndex() + j > validRangeCeiling_.jIndex();
    }
    if(checkK_) 
    {
        isKPositionOutOfRange = 
            currentPosition_.kIndex() + k < validRangeFloor_.kIndex() || 
            currentPosition_.kIndex() + k > validRangeCeiling_.kIndex();
    }
    
    // print error message and assert in case position is out of range
    if(isIPositionOutOfRange || isJPositionOutOfRange || isKPositionOutOfRange)
    {
#ifdef __OPENMP_BACKEND__
        #pragma omp critical
#endif
        {
            if(isIPositionOutOfRange)
            {
                std::cerr 
                    << "PositionChecker -> i position " << currentPosition_.iIndex() + i
                    << " is out of range [" << validRangeFloor_.iIndex() << ", " << validRangeCeiling_.iIndex() + 1 << ")" << std::endl;
            }
            if(isJPositionOutOfRange)
            {
                std::cerr 
                    << "PositionChecker -> j position " << currentPosition_.jIndex() + j
                    << " is out of range [" << validRangeFloor_.jIndex() << ", " << validRangeCeiling_.jIndex() + 1 << ")" << std::endl;
            }
            if(isKPositionOutOfRange)
            {
                std::cerr 
                    << "PositionChecker -> k position " << currentPosition_.kIndex() + k 
                    << " is out of range [" << validRangeFloor_.kIndex() << ", " << validRangeCeiling_.kIndex() + 1 << ")" << std::endl;
            }
        }

        // assert in case the position is out of range
        assert(false);
    }
    
    // update the bounding box
    IJKIndex accessIndex;
    accessIndex.Init(currentPosition_.iIndex() + i, currentPosition_.jIndex() + j, currentPosition_.kIndex() + k);
    boundingBoxFloor_.CommonMin(accessIndex);
    boundingBoxCeiling_.CommonMax(accessIndex);
}


