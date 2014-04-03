#include "gtest/gtest.h"
#include "SharedInfrastructure.h"

// test JokerDataField and Iterators
class JokerDataFieldUnittest : public ::testing::Test
{
protected:
    // define arrays
    SwapDataField<IJKRealField> testField_;
    JokerDataField<IJKRealField> testJokerDataField_;
    std::string name_;    

    virtual void SetUp()
    {
        name_ = "joker";
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
        testJokerDataField_.Init(name_, testField_.in());
    };
};

TEST_F(JokerDataFieldUnittest, SetField)
{
    testField_.out()(0,0,0) = 77;
    testField_.Swap();
    testField_.out()(0,0,0) = 99;

    // test if values have been set correctly
    ASSERT_EQ(name_, testJokerDataField_.name());

    testJokerDataField_.set_dataField(testField_.in());
    ASSERT_EQ(77, testJokerDataField_.dataField()(0,0,0));

    // double buffer swaps the actual data storage
    testField_.Swap();
    ASSERT_EQ(99, testJokerDataField_.dataField()(0,0,0));

    // assign a new field
    testJokerDataField_.set_dataField(testField_.out());
    ASSERT_EQ(77, testJokerDataField_.dataField()(0,0,0));

    // test field assign
    testJokerDataField_.dataField()(0,0,0) = 99;
    ASSERT_EQ(99, testJokerDataField_.dataField()(0,0,0));
}

  
