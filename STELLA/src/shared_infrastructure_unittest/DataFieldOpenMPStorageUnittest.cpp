#include <vector>
#include "gtest/gtest.h"
#include "DataFieldOpenMPStorage.h"

typedef DataFieldAlignment<cDimK, 1> Alignment;
typedef DataFieldIJBoundary<-1, 2, -3, 4> Boundary;

// test DataFieldStorage and Iterators
class DataFieldOpenMPStorageUnittest : public ::testing::Test
{
protected:
    // test data
    typedef double DataType3D;
    typedef float DataType2D;
    typedef int DataType1D;

    typedef DataFieldStorageFormat<Boundary, StorageOrder::JIK, Alignment> IJKStorageFormat;
    typedef DataFieldStorageFormat<Boundary, StorageOrder::JI, Alignment> IJStorageFormat;
    typedef DataFieldStorageFormat<Boundary, StorageOrder::I, Alignment> IStorageFormat;

    DataFieldOpenMPStorage<DataType3D, IJKStorageFormat> storage3D_;
    DataFieldOpenMPStorage<DataType2D, IJStorageFormat> storage2D_;
    DataFieldOpenMPStorage<DataType1D, IStorageFormat> storage1D_;
 
    IJKSize calculationDomain_;
    KBoundary kBoundary_;

    virtual void SetUp()
    {
        calculationDomain_.Init(12, 14, 60);
        kBoundary_.Init(-1, 2);

        storage3D_.Init(calculationDomain_, kBoundary_);
        storage2D_.Init(calculationDomain_, kBoundary_);
        storage1D_.Init(calculationDomain_, kBoundary_);
    };

    template<typename TStorage>
    void randomStorageInit(TStorage& storage)
    {
        IJKSize allocSize = storage.allocatedSize();
        typename TStorage::StorageIteratorType iter = storage.originIterator();
        iter.Advance(
            -storage.originOffset().iIndex(), 
            -storage.originOffset().jIndex(), 
            -storage.originOffset().kIndex()
        );
        
        for(int i = 0; i < allocSize.iSize(); ++i)
        {
            for(int j = 0; j < allocSize.jSize(); ++j)
            {
                for(int k = 0; k < allocSize.kSize(); ++k)
                {
                    iter.At(i,j,k) = static_cast<typename TStorage::ValueType>(rand());
                }
            }
        }
    }
};

TEST_F(DataFieldOpenMPStorageUnittest, Copy)
{
    // initialize data field with random numbers
    randomStorageInit(storage3D_);
    randomStorageInit(storage2D_);
    randomStorageInit(storage1D_);
    
    // copy fields
    DataFieldOpenMPStorage<DataType3D, IJKStorageFormat> storage3DCopyConstruct(storage3D_);
    DataFieldOpenMPStorage<DataType2D, IJStorageFormat> storage2DCopyConstruct(storage2D_);
    DataFieldOpenMPStorage<DataType1D, IStorageFormat> storage1DCopyConstruct(storage1D_);

    DataFieldOpenMPStorage<DataType3D, IJKStorageFormat> storage3DAssign;
    DataFieldOpenMPStorage<DataType2D, IJStorageFormat> storage2DAssign;
    DataFieldOpenMPStorage<DataType1D, IStorageFormat> storage1DAssign;

    storage3DAssign = storage3D_;
    storage2DAssign = storage2D_;
    storage1DAssign = storage1D_;

    // assert copies are made
    ASSERT_NE(storage3D_.pStorageBase(), storage3DCopyConstruct.pStorageBase());
    ASSERT_NE(storage3D_.pStorageBase(), storage3DAssign.pStorageBase());
    ASSERT_NE(storage2D_.pStorageBase(), storage2DCopyConstruct.pStorageBase());
    ASSERT_NE(storage2D_.pStorageBase(), storage2DAssign.pStorageBase());
    ASSERT_NE(storage1D_.pStorageBase(), storage1DCopyConstruct.pStorageBase());
    ASSERT_NE(storage1D_.pStorageBase(), storage1DAssign.pStorageBase());

    // check content is copied
    for(int i = -1; i < 12-2; ++i) 
    {
        for(int j = -3; j < 14+4; ++j) 
        {
            for(int k = -1; k < 60+2; ++k) 
            {
                ASSERT_EQ(storage3D_.constOriginIterator().At(i,j,k), storage3DCopyConstruct.constOriginIterator().At(i,j,k));
                ASSERT_EQ(storage3D_.constOriginIterator().At(i,j,k), storage3DAssign.constOriginIterator().At(i,j,k));
                ASSERT_EQ(storage2D_.constOriginIterator().At(i,j,1), storage2DCopyConstruct.constOriginIterator().At(i,j,k));
                ASSERT_EQ(storage2D_.constOriginIterator().At(i,j,1), storage2DAssign.constOriginIterator().At(i,j,k));
                ASSERT_EQ(storage1D_.constOriginIterator().At(i,1,1), storage1DCopyConstruct.constOriginIterator().At(i,j,k));
                ASSERT_EQ(storage1D_.constOriginIterator().At(i,1,1), storage1DAssign.constOriginIterator().At(i,j,k));
            }
        }
    }
}

TEST_F(DataFieldOpenMPStorageUnittest, ExternalStorage)
{
    // setup external memory
    IJKSize paddedSize = storage3D_.paddedSize();
    int size = paddedSize.iSize() * paddedSize.jSize() * paddedSize.kSize();
    std::vector<DataType3D> memory1, memory2;
    memory1.resize(size, 11.0);
    memory2.resize(size, 22.0);

    // define external storage
    ExternalStorage<DataType3D> externalStorage1, externalStorage2;
    externalStorage1.Init(&memory1[0], storage3D_.paddedSize());
    externalStorage2.Init(&memory2[0], storage3D_.paddedSize());
    
    // storage with external storage
    DataFieldOpenMPStorage<DataType3D, IJKStorageFormat> storage;
    storage.Init(externalStorage1, calculationDomain_, kBoundary_);
    
    // check the memory is ok
    ASSERT_EQ(&memory1[0], storage.pStorageBase());
    ASSERT_EQ(11.0, *storage.pStorageBase());

    // set the second external storage
    storage.SetExternalStorage(externalStorage2);

    // check the memory is ok
    ASSERT_EQ(&memory2[0], storage.pStorageBase());
    ASSERT_EQ(22.0, *storage.pStorageBase());
}

TEST_F(DataFieldOpenMPStorageUnittest, Access)
{
    // init storage
    randomStorageInit(storage3D_);

    // get 2 iterators
    DataFieldOpenMPStorage<DataType3D, IJKStorageFormat>::StorageIteratorType iter1 = storage3D_.originIterator();
    DataFieldOpenMPStorage<DataType3D, IJKStorageFormat>::StorageIteratorType iter2 = storage3D_.originIterator();
    
    // set the base value to 99
    *storage3D_.pStorageBase() = 99.0;

    // check iterator advance
    iter1.Advance(-1, -3, -1); // move to storage base
    ASSERT_EQ(99.0, iter1.At(0, 0, 0));

    // check read access
    ASSERT_EQ(99.0, iter2.At(-1, -3, -1)); // access storage base

    // check write access
    iter1.At(0, 0, 0) += 3.0;
    ASSERT_EQ(99.0 + 3.0, iter2.At(-1, -3, -1)); 
}

TEST_F(DataFieldOpenMPStorageUnittest, MemorizedPosition)
{
    DataFieldOpenMPStorage<DataType3D, IJKStorageFormat>::StorageIteratorType iter = storage3D_.originIterator();
    const DataType3D* ptr = &iter.At(0,0,0);
    
    iter.Advance(1,2,3);
    ASSERT_NE(ptr, &iter.At(0,0,0)); // current position different from origin

    iter.RestoreMemorizedPosition(); // reset current position to origin
    ASSERT_EQ(ptr, &iter.At(0,0,0));

    iter.AdvanceMemorizedPosition(1,2,3); // move the memorized position
    iter.RestoreMemorizedPosition(); 
    ASSERT_NE(ptr, &iter.At(0,0,0)); 
    ASSERT_EQ(ptr, &iter.At(-1,-2,-3)); 
}
