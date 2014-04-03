#pragma once

#include "Enums.h"
#include "IJKIndex.h"

/**
* @class IteratorPositionChecker
* Class tracking the position of an iterator, used in debug build to check if there are access violations.
*/
class PositionChecker
{
public:
    PositionChecker();
    ~PositionChecker();

    PositionChecker(const PositionChecker& other);
    PositionChecker& operator= (const PositionChecker& other);

    /**
    * Initialize current position and set new valid range for access range check in debug mode
    * @param floor smallest allowed index
    * @param ceiling largest allowed index
    */
    void Init(const IJKIndex& floor, const IJKIndex& ceiling, const bool checkI, const bool checkJ, const bool checkK);
    
    /**
    * @return the current index
    */
    const IJKIndex& currentPosition() const;

    /**
    * Set new valid range for access range check in debug mode
    * @param floor smallest allowed index
    * @param ceiling largest allowed index
    */
    void SetValidRange(const IJKIndex& floor, const IJKIndex& ceiling);
   
    /**
    * Method resetting the bounding box which is the minimal box 
    * encapsulating all data accesses so far
    */
    void ResetBoundingBox();

    /**
    * Merge the bounding boxes of two position checkers
    * @param other position checker
    */
    void MergeBoundingBox(const PositionChecker& other);

    /**
    * Set current to current of other position checker
    * @param other position checker
    */
    void SetPosition(const PositionChecker& other);

    /**
    * Advance current
    */
    void Advance(const int i, const int j, const int k);

    /**
    * Advance in one dimension
    */
    void Advance(Dimension dim, const int d);

    /**
    * Check access index
    */
    void CheckPositionAtOffset(const int i, const int j, const int k);
    
    /**
    * @return valid range floor
    */
    const IJKIndex& validRangeFloor() const;

    /**
    * @return valid range ceiling
    */
    const IJKIndex& validRangeCeiling() const;
    
    /**
    * @return bounding box floor
    */
    const IJKIndex& boundingBoxFloor() const;

    /**
    * @return bounding box floor
    */
    const IJKIndex& boundingBoxCeiling() const;

private:
    IJKIndex    currentPosition_;
    IJKIndex    validRangeFloor_;
    IJKIndex    validRangeCeiling_;
    IJKIndex    boundingBoxFloor_;
    IJKIndex    boundingBoxCeiling_;
    bool        checkI_;
    bool        checkJ_;
    bool        checkK_;
};
  
