#include <vector>
#include "gtest/gtest.h"
#include "SharedInfrastructure.h"

// test DataField and Iterators
class DataFieldOpenMPUnittest : public ::testing::Test
{
protected:
    static const int iCalcDomainSize;
    static const int jCalcDomainSize;
    static const int kCalcDomainSize;
    static const int kMinusOffset;
    static const int kPlusOffset;

    IJKRealField test3DField_;
    IJKSize calcDomainSize_;
    KBoundary boundary_;
    std::string name_;

    virtual void SetUp()
    {
        name_ = "testField";

        // initialization for further tests
        calcDomainSize_.Init(iCalcDomainSize, jCalcDomainSize, kCalcDomainSize);
        boundary_.Init(kMinusOffset, kPlusOffset);
        
        test3DField_.Init(name_, calcDomainSize_, boundary_);
    };

};
const int DataFieldOpenMPUnittest::iCalcDomainSize = 53;
const int DataFieldOpenMPUnittest::jCalcDomainSize = 59;
const int DataFieldOpenMPUnittest::kCalcDomainSize = 61;
const int DataFieldOpenMPUnittest::kMinusOffset = -11;
const int DataFieldOpenMPUnittest::kPlusOffset = 13;

TEST_F(DataFieldOpenMPUnittest, Storage)
{
    // test that the storage points to same data as the field
    test3DField_(0,0,0) = 77;
    test3DField_(3,4,5) = 55;

    ASSERT_EQ(77, test3DField_(0,0,0));
    ASSERT_EQ(77, test3DField_.storage().originIterator().At(0,0,0));
    ASSERT_EQ(55, test3DField_(3,4,5));
    ASSERT_EQ(55, test3DField_.storage().originIterator().At(3,4,5));
}

TEST_F(DataFieldOpenMPUnittest, InitWithExternalStorage)
{
    // allocate external memory
    IJKSize paddedSize = test3DField_.storage().paddedSize();
    int size = paddedSize.iSize() * paddedSize.jSize() * paddedSize.kSize();

    std::vector<Real> memory;
    memory.resize(size, 33.0);

    // define external storage
    ExternalStorage<Real> externalStorage;
    externalStorage.Init(&memory[0], paddedSize);

    IJKRealField dataField;
    dataField.Init("dataField", externalStorage, calcDomainSize_, boundary_);

    // check the pointer
    ASSERT_EQ(&memory[0], dataField.storage().pStorageBase()); 
    ASSERT_EQ(33.0, dataField(0,0,0));
    dataField(0,0,0) = 77.0;
    ASSERT_EQ(77.0, memory[&dataField(0,0,0) - dataField.storage().pStorageBase()]);
        
    // allocate new memory for the field
    dataField.Init("dataField", calcDomainSize_, boundary_);
    ASSERT_NE(&memory[0], dataField.storage().pStorageBase()); 

    // check the vector is still valid
    ASSERT_EQ(33.0, memory[0]);
}

TEST_F(DataFieldOpenMPUnittest, SetExternalStorage)
{
    // allocate external memory
    IJKSize paddedSize = test3DField_.storage().paddedSize();
    int size = paddedSize.iSize() * paddedSize.jSize() * paddedSize.kSize();

    std::vector<Real> memory1, memory2;
    memory1.resize(size, 33.0);
    memory2.resize(size, 77.0);
    
    // define external storage
    ExternalStorage<Real> externalStorage1, externalStorage2;
    externalStorage1.Init(&memory1[0], paddedSize);
    externalStorage2.Init(&memory2[0], paddedSize);

    // normal data field initialization
    IJKRealField dataField;
    dataField.Init("dataField", calcDomainSize_, boundary_);

    // test setting the first storage
    dataField.SetExternalStorage(externalStorage1);
    ASSERT_EQ(33.0, dataField(0,0,0));
    dataField(0,0,0) = 44.0; // set the origin

    // test setting the second storage
    dataField.SetExternalStorage(externalStorage2);
    ASSERT_EQ(77.0, dataField(0,0,0));
    dataField(0,0,0) = 88.0; // set the origin

    // go back to the first storage and make sure origin was set
    dataField.SetExternalStorage(externalStorage1);
    ASSERT_EQ(44.0, dataField(0,0,0));
    ASSERT_EQ(33.0, dataField(1,1,1));

    // go back to the second storage and make sure origin was set
    dataField.SetExternalStorage(externalStorage2);
    ASSERT_EQ(88.0, dataField(0,0,0));
    ASSERT_EQ(77.0, dataField(1,1,1));
}