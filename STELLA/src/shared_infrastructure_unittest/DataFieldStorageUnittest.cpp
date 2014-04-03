#include "gtest/gtest.h"
#include "SharedInfrastructure.h"

// test DataFieldStorage and Iterators
class DataFieldStorageUnittest : public ::testing::Test
{
protected:
    IJKSize calculationDomain_;
    KBoundary kBoundary_;

    virtual void SetUp()
    {
        calculationDomain_.Init(12, 14, 60);
        kBoundary_.Init(-1, 2);
    };
};

TEST_F(DataFieldStorageUnittest, Getters)
{
    typedef DataFieldStorageFormat<DataFieldIJBoundary<-2, 1, -3, 4>, StorageOrder::JIK, DataFieldAlignment<cDimI, 1> > IJKStorage;
    typedef DataFieldStorageFormat<DataFieldIJBoundary<-2, 1, -3, 4>, StorageOrder::JI, DataFieldAlignment<cDimI, 1> > IJStorage;
    typedef DataFieldStorageFormat<DataFieldIJBoundary<-2, 1, -3, 4>, StorageOrder::I, DataFieldAlignment<cDimI, 1> > IStorage;

    // define storage
#ifdef __CUDA_BACKEND__
    DataFieldCUDAStorage<Real, IJKStorage> storage3D;
    DataFieldCUDAStorage<Real, IJStorage>  storage2D;
    DataFieldCUDAStorage<Real, IStorage>  storage1D;
#else
    DataFieldOpenMPStorage<Real, IJKStorage> storage3D;
    DataFieldOpenMPStorage<Real, IJStorage>  storage2D;
    DataFieldOpenMPStorage<Real, IStorage>  storage1D;
#endif

    // setup sizes
    IJKSize totalSize;
    totalSize.Init(
        calculationDomain_.iSize() + 2 + 1,
        calculationDomain_.jSize() + 3 + 4,
        calculationDomain_.kSize() + 1 + 2
    );

    // initialization for further tests
    storage3D.Init(calculationDomain_, kBoundary_);
    storage2D.Init(calculationDomain_, kBoundary_);
    storage1D.Init(calculationDomain_, kBoundary_);

    // test initialization of 3D field
    ASSERT_TRUE(totalSize == storage3D.size());
    ASSERT_TRUE(totalSize == storage3D.allocatedSize());
    ASSERT_TRUE(totalSize == storage3D.paddedSize());
    ASSERT_EQ(3, storage3D.rank());
    ASSERT_EQ(2, storage3D.originOffset().iIndex());
    ASSERT_EQ(3, storage3D.originOffset().jIndex());
    ASSERT_EQ(1, storage3D.originOffset().kIndex());

    // test initialization of 2D field
    ASSERT_TRUE(totalSize == storage2D.size());
    ASSERT_EQ(totalSize.iSize(), storage2D.allocatedSize().iSize());
    ASSERT_EQ(totalSize.jSize(), storage2D.allocatedSize().jSize());
    ASSERT_EQ(1, storage2D.allocatedSize().kSize());
    ASSERT_EQ(2, storage2D.rank());

    // test initialization of 1D field
    ASSERT_TRUE(totalSize == storage1D.size());
    ASSERT_EQ(totalSize.iSize(), storage1D.allocatedSize().iSize());
    ASSERT_EQ(1, storage1D.allocatedSize().jSize());
    ASSERT_EQ(1, storage1D.allocatedSize().kSize());
    ASSERT_EQ(1, storage1D.rank());   
}

TEST_F(DataFieldStorageUnittest, Alignment)
{
    typedef DataFieldStorageFormat<DataFieldIJBoundary<-2, 1, -3, 4>, StorageOrder::KJI, DataFieldAlignment<cDimI, 5> > AlignedI;
    typedef DataFieldStorageFormat<DataFieldIJBoundary<-2, 1, -3, 4>, StorageOrder::KIJ, DataFieldAlignment<cDimJ, 17> > AlignedJ;
    typedef DataFieldStorageFormat<DataFieldIJBoundary<-2, 1, -3, 4>, StorageOrder::IJK, DataFieldAlignment<cDimK, 11> > AlignedK;

    // define storage
#ifdef __CUDA_BACKEND__
    DataFieldCUDAStorage<Real, AlignedI> storageI;
    DataFieldCUDAStorage<Real, AlignedJ>  storageJ;
    DataFieldCUDAStorage<Real, AlignedK>  storageK;
#else
    DataFieldOpenMPStorage<Real, AlignedI> storageI;
    DataFieldOpenMPStorage<Real, AlignedJ>  storageJ;
    DataFieldOpenMPStorage<Real, AlignedK>  storageK;
#endif

    storageI.Init(calculationDomain_, kBoundary_);
    storageJ.Init(calculationDomain_, kBoundary_);
    storageK.Init(calculationDomain_, kBoundary_);
     
    ASSERT_EQ((size_t)0, (size_t)(storageI.pStorageBase() + 2) % (5 * sizeof(double)));
    ASSERT_EQ(storageI.paddedSize().iSize(), 15); /* 12 + 2 + 1 = 15 ok */

    ASSERT_EQ((size_t)0, (size_t)(storageJ.pStorageBase() + 3) % (17 * sizeof(double)));
    ASSERT_EQ(storageJ.paddedSize().jSize(), 34); /* 13 + 3 + 4 = 20 -> 2*17 = 34 */

    ASSERT_EQ((size_t)0, (size_t)(storageK.pStorageBase() + 1) % (11 * sizeof(double)));
    ASSERT_EQ(storageK.paddedSize().kSize(), 66); /* 60 + 1 + 2 = 63 -> 6*11 = 66 */
}
