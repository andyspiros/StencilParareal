#include "gtest/gtest.h"
#include "SharedInfrastructure.h"
#include "DataFieldOpenMP.h"

// test DataField and Iterators
class DataFieldUnittest : public ::testing::Test
{
protected:
    std::string name;
    static const int iCalcDomainSize;
    static const int jCalcDomainSize;
    static const int kCalcDomainSize;
    static const int kMinusOffset;
    static const int kPlusOffset;
    
    typedef DataFieldOpenMP<Real, DataFieldStorageFormat<DataFieldIJBoundary<-2, 1, -3, 4>, StorageOrder::JIK, DataFieldAlignment<cDimK, 1> > > TestField;

    // define arrays
    TestField test3DField_;
    TestField swap3DField_;

    DataFieldUnittest() : name("testField") {};

    virtual void SetUp()
    {
        // initialization for further tests
        IJKSize calcDomainSize;
        calcDomainSize.Init(iCalcDomainSize, jCalcDomainSize, kCalcDomainSize);
        KBoundary boundary;
        boundary.Init(kMinusOffset, kPlusOffset);
        test3DField_.Init(name, calcDomainSize, boundary);
        swap3DField_.Init(name, calcDomainSize, boundary);
    };

};
const int DataFieldUnittest::iCalcDomainSize = 53;
const int DataFieldUnittest::jCalcDomainSize = 59;
const int DataFieldUnittest::kCalcDomainSize = 61;
const int DataFieldUnittest::kMinusOffset = -11;
const int DataFieldUnittest::kPlusOffset = 13;

TEST_F(DataFieldUnittest, AccessOperator)
{
    // test storage
    test3DField_(0,0,0) = 33;
    test3DField_(1,2,3) = 77;

    ASSERT_EQ(77, test3DField_(1,2,3));
    ASSERT_EQ(33, test3DField_(0,0,0));
}

TEST_F(DataFieldUnittest, SwapWith)
{
    // set cell 
    swap3DField_(0,0,0) = 77;
    test3DField_(0,0,0) = 33;

    ASSERT_EQ(77, swap3DField_(0,0,0));
    ASSERT_EQ(33, test3DField_(0,0,0));

    test3DField_.SwapWith(swap3DField_);

    ASSERT_EQ(33, swap3DField_(0,0,0));
    ASSERT_EQ(77, test3DField_(0,0,0));
}

TEST_F(DataFieldUnittest, Getters)
{
    // test initialization of 3D field
    IJKSize calcDomSize = test3DField_.calculationDomain();
    ASSERT_EQ(iCalcDomainSize, calcDomSize.iSize());
    ASSERT_EQ(jCalcDomainSize, calcDomSize.jSize());
    ASSERT_EQ(kCalcDomainSize, calcDomSize.kSize());

    IJKBoundary boundary = test3DField_.boundary();
    ASSERT_EQ(-2, boundary.iMinusOffset());
    ASSERT_EQ(1, boundary.iPlusOffset());
    ASSERT_EQ(-3, boundary.jMinusOffset());
    ASSERT_EQ(4, boundary.jPlusOffset());
    ASSERT_EQ(kMinusOffset, boundary.kMinusOffset());
    ASSERT_EQ(kPlusOffset, boundary.kPlusOffset());
    ASSERT_EQ(name, test3DField_.name());
    ASSERT_EQ(3, test3DField_.storage().rank());

    IJKIndex originOffset = test3DField_.storage().originOffset();
    ASSERT_EQ(-2, -originOffset.iIndex());
    ASSERT_EQ(-3, -originOffset.jIndex());
    ASSERT_EQ(kMinusOffset, -originOffset.kIndex());
}

TEST_F(DataFieldUnittest, Copy)
{
    test3DField_(0,0,0) = 77;
    // copy field
    TestField copyField;
    ASSERT_FALSE(copyField.isInitialized());
    copyField = test3DField_;
    ASSERT_TRUE(copyField.isInitialized());

    // compare data
    ASSERT_EQ(77, copyField(0,0,0));
    copyField(0,0,0) = 33;
    // check original field is not modified
    ASSERT_EQ(33, copyField(0,0,0));
    ASSERT_EQ(77, test3DField_(0,0,0));

    // compare size etc.
    IJKSize calcDomSize = copyField.calculationDomain();
    ASSERT_EQ(iCalcDomainSize, calcDomSize.iSize());
    ASSERT_EQ(jCalcDomainSize, calcDomSize.jSize());
    ASSERT_EQ(kCalcDomainSize, calcDomSize.kSize());

    IJKBoundary boundary = copyField.boundary();
    ASSERT_EQ(-2, boundary.iMinusOffset());
    ASSERT_EQ(1, boundary.iPlusOffset());
    ASSERT_EQ(-3, boundary.jMinusOffset());
    ASSERT_EQ(4, boundary.jPlusOffset());
    ASSERT_EQ(kMinusOffset, boundary.kMinusOffset());
    ASSERT_EQ(kPlusOffset, boundary.kPlusOffset());
    ASSERT_EQ(name, test3DField_.name());
}
