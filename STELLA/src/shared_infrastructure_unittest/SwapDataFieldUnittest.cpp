#include "gtest/gtest.h"
#include "SharedInfrastructure.h"

// test SwapDataField and Iterators
class SwapDataFieldUnittest : public ::testing::Test
{
protected:
    std::string name_;
    // define arrays
    SwapDataField<IJKRealField> testField_;

    virtual void SetUp()
    {
        name_ = "testField";
        const int iCalcDomainSize = 53;
        const int jCalcDomainSize = 59;
        const int kCalcDomainSize = 61;
        const int kMinusOffset = -11;
        const int kPlusOffset = 13;

        // initialization for further tests
        IJKSize calcDomainSize;
        calcDomainSize.Init(iCalcDomainSize, jCalcDomainSize, kCalcDomainSize);
        KBoundary boundary;
        boundary.Init(kMinusOffset, kPlusOffset);
        testField_.Init(name_, calcDomainSize, boundary);
    };
};

TEST_F(SwapDataFieldUnittest, CopySwap)
{
    const int testValue0 = 77;
    const int testValue1 = 99;
    testField_.out()(0,0,0) = testValue0;
    testField_.Swap();
    testField_.out()(0,0,0) = testValue1;

    // check the name property
    ASSERT_EQ(name_ + "_in", testField_.in().name());
    ASSERT_EQ(name_ + "_out", testField_.out().name());
    
    // test if values have been set correctly
    ASSERT_EQ(testValue0, testField_.in()(0,0,0));
    ASSERT_EQ(testValue1, testField_.out()(0,0,0));
    
    // copy fields
    SwapDataField<IJKRealField> copyField = testField_;

    // compare data
    ASSERT_EQ(testField_.in()(0,0,0), copyField.in()(0,0,0));
    ASSERT_EQ(testField_.out()(0,0,0), copyField.out()(0,0,0));

    // swap copy and compare
    copyField.Swap();
    ASSERT_EQ(testField_.in()(0,0,0), copyField.out()(0,0,0));
    ASSERT_EQ(testField_.out()(0,0,0), copyField.in()(0,0,0));

    // check if it is really a copy
    copyField.out()(0,0,0) = copyField.in()(0,0,0);
    ASSERT_NE(testField_.in()(0,0,0), copyField.out()(0,0,0));
    ASSERT_EQ(testField_.out()(0,0,0), copyField.out()(0,0,0));
    ASSERT_EQ(testField_.out()(0,0,0), copyField.in()(0,0,0));
}

TEST_F(SwapDataFieldUnittest, MemorizeSwapPosition)
{
    // initially points to nowhere
    ASSERT_EQ(NULL, testField_.pMemorizedSwapPositionStorage());

    testField_.MemorizeSwapPosition();

    // check the memorized position storage points to the in field 
    ASSERT_EQ(&testField_.in().storage(), testField_.pMemorizedSwapPositionStorage());
    
    testField_.Swap();

    // check after swapping the memorized storage points to the out field
    ASSERT_EQ(&testField_.out().storage(), testField_.pMemorizedSwapPositionStorage());
}  


