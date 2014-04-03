#include "gtest/gtest.h"
#include "Definitions.h"
#include "SharedInfrastructure.h"

// test DataField and Iterators
class DataFieldCUDAUnittest : public ::testing::Test
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
const int DataFieldCUDAUnittest::iCalcDomainSize = 64;
const int DataFieldCUDAUnittest::jCalcDomainSize = 24;
const int DataFieldCUDAUnittest::kCalcDomainSize = 61;
const int DataFieldCUDAUnittest::kMinusOffset = -11;
const int DataFieldCUDAUnittest::kPlusOffset = 13;

template<
    typename TStorageStrides,
    typename TStoragePointer>
__KERNEL__
void addValueToField(const TStorageStrides storageStrides, TStoragePointer storagePointer, const int kSize, Real increment)
{
    DataFieldCUDAStorageIndex<TStorageStrides> index;
    index.Init(storageStrides, blockIdx.x, blockIdx.y, threadIdx.x, threadIdx.y, 0);
    for(int k = 0; k < kSize; ++k)
    {
        storagePointer.Center(index) += increment;
        index.Advance(storageStrides, 0, 0, 1);
    }
}

TEST_F(DataFieldCUDAUnittest, Access)
{
    // test that the storage points to same data as the field
    test3DField_(0, 0, 0) = 77;
    test3DField_(3, 4, 5) = 55;
    test3DField_(iCalcDomainSize-1, jCalcDomainSize-1, kCalcDomainSize-1) = 77;

    ASSERT_EQ(77,test3DField_(0, 0, 0));
    ASSERT_EQ(55,test3DField_(3, 4, 5));
    ASSERT_EQ(77,test3DField_(iCalcDomainSize-1, jCalcDomainSize-1, kCalcDomainSize-1));
    
    // set elements outside calculation domain (kernel shouldn't modify them)
    test3DField_(-1, 0, 0) = 55;
    test3DField_(0, -1, 0) = 55;
    test3DField_(0, 0, -1) = 55;
    test3DField_(iCalcDomainSize, jCalcDomainSize-1, kCalcDomainSize-1) = 55;
    test3DField_(iCalcDomainSize-1, jCalcDomainSize, kCalcDomainSize-1) = 55;
    test3DField_(iCalcDomainSize-1, jCalcDomainSize-1, kCalcDomainSize) = 55;

    // add a value and test again
    test3DField_.SynchronizeDeviceStorage();
    
    dim3 threadsPerBlock(32, 2); 
    dim3 numberOfBlocks(iCalcDomainSize/32, jCalcDomainSize/2);
    ASSERT_TRUE(iCalcDomainSize % 32 == 0 && jCalcDomainSize % 2 == 0);

    typedef DataFieldCUDAStoragePointer<ScalarPointer<Real, typename IJKRealField::StorageFormat, BlockSize<32,2>, cReadWrite> > StoragePointerType;
    typedef StoragePointerType::StorageStridesType StorageStridesType;
    StoragePointerType originPointer;
    StorageStridesType storageStrides;
    
    ASSERT_FALSE(originPointer.PointsToOrigin(test3DField_.deviceStorage().pStorageOrigin()));
    test3DField_.deviceStorage().InitializeStoragePointerToOrigin(originPointer);
    ASSERT_TRUE(originPointer.PointsToOrigin(test3DField_.deviceStorage().pStorageOrigin()));
    test3DField_.deviceStorage().TryToInitializeStorageStrides(storageStrides);

    addValueToField<<<numberOfBlocks,threadsPerBlock>>>(storageStrides, originPointer, kCalcDomainSize, 3.0);

    // synchronization should be done by device 
    ASSERT_EQ(77+3, test3DField_(0, 0, 0));
    ASSERT_EQ(55+3, test3DField_(3, 4, 5));
    ASSERT_EQ(77+3, test3DField_(iCalcDomainSize-1, jCalcDomainSize-1, kCalcDomainSize-1));
    ASSERT_EQ(55, test3DField_(-1, 0, 0));
    ASSERT_EQ(55, test3DField_(0, -1, 0));
    ASSERT_EQ(55, test3DField_(0, 0, -1));
    ASSERT_EQ(55, test3DField_(iCalcDomainSize, jCalcDomainSize-1, kCalcDomainSize-1));
    ASSERT_EQ(55, test3DField_(iCalcDomainSize-1, jCalcDomainSize, kCalcDomainSize-1));
    ASSERT_EQ(55, test3DField_(iCalcDomainSize-1, jCalcDomainSize-1, kCalcDomainSize));
}

TEST_F(DataFieldCUDAUnittest, Copy)
{
    IJKRealField copy;

    // test that the storage points to same data as the field
    test3DField_(0, 0, 0) = 77;
    test3DField_(3, 4, 5) = 55;
    test3DField_(iCalcDomainSize-1, jCalcDomainSize-1, kCalcDomainSize-1) = 77;
    copy = test3DField_;

    // check host copy
    ASSERT_EQ(77, copy(0, 0, 0));
    ASSERT_EQ(55, copy(3, 4, 5));
    ASSERT_EQ(77, copy(iCalcDomainSize-1, jCalcDomainSize-1, kCalcDomainSize-1));
    
    // test it is a real copy
    copy(0, 0, 0) = 0;
    ASSERT_EQ(77, test3DField_(0, 0, 0));
    
    // add a value and test again
    test3DField_.SynchronizeDeviceStorage();
    
    dim3 threadsPerBlock(32, 2); 
    dim3 numberOfBlocks(iCalcDomainSize/32, jCalcDomainSize/2);
    ASSERT_TRUE(iCalcDomainSize % 32 == 0 && jCalcDomainSize % 2 == 0);

    typedef DataFieldCUDAStoragePointer<ScalarPointer<Real, typename IJKRealField::StorageFormat, BlockSize<32,2>, cReadWrite> > StoragePointerType;
    typedef StoragePointerType::StorageStridesType StorageStridesType;
    StoragePointerType originPointer;
    StorageStridesType storageStrides;
    test3DField_.deviceStorage().InitializeStoragePointerToOrigin(originPointer);
    test3DField_.deviceStorage().TryToInitializeStorageStrides(storageStrides);
    addValueToField<<<numberOfBlocks,threadsPerBlock>>>(storageStrides, originPointer, kCalcDomainSize, 3.0);
   
    // check device copy
    copy = test3DField_;
    ASSERT_EQ(77+3, copy(0, 0, 0));
    ASSERT_EQ(55+3, copy(3, 4, 5));
    ASSERT_EQ(77+3, copy(iCalcDomainSize-1, jCalcDomainSize-1, kCalcDomainSize-1));
}

TEST_F(DataFieldCUDAUnittest, InitWithExternalStorage)
{
    // define a data field without alignment
    typedef DataFieldCUDA<
        Real, 
        DataFieldStorageFormat<
            CUDAIJBoundary, 
            StorageOrder::KJI, 
            DataFieldAlignment<cDimI, 1> 
        > 
    > IJKRealFieldNoAlignment;

    // compute the size of the external using a dummy field
    IJKRealFieldNoAlignment sizeField;
    sizeField.Init("sizeField", calcDomainSize_, boundary_);
    
    IJKSize paddedSize = sizeField.storage().paddedSize();
    int size = paddedSize.iSize() * paddedSize.jSize() * paddedSize.kSize();

    // allocate the external stroage
    std::vector<Real> memory;
    memory.resize(size, 33.0);

    Real* pDeviceMemory;
    cudaMalloc((void**)&pDeviceMemory, size * sizeof(Real));
    cudaMemcpy(pDeviceMemory, &memory[0], size * sizeof(Real), cudaMemcpyHostToDevice);
    
    ExternalStorage<Real> externalStorage;
    externalStorage.Init(pDeviceMemory, paddedSize);

    // initialize the external memory field
    IJKRealFieldNoAlignment dataField;
    dataField.Init("dataField", externalStorage, calcDomainSize_, boundary_);

    // check the pointer
    ASSERT_EQ(pDeviceMemory, dataField.storage().pStorageBase()); 
    ASSERT_EQ(33.0, dataField(0, 0, 0)); 

    // allocate new memory for the field
    dataField.Init("dataField", calcDomainSize_, boundary_);
    ASSERT_NE(pDeviceMemory, dataField.storage().pStorageBase()); 

    cudaFree(pDeviceMemory);
}

TEST_F(DataFieldCUDAUnittest, SetExternalStorage)
{
    // define a data field without alignment
    typedef DataFieldCUDA<
        Real, 
        DataFieldStorageFormat<
            CUDAIJBoundary, 
            StorageOrder::KJI, 
            DataFieldAlignment<cDimI, 1> 
        > 
    > IJKRealFieldNoAlignment;

    // compute the size of the external using a dummy field
    IJKRealFieldNoAlignment sizeField;
    sizeField.Init("sizeField", calcDomainSize_, boundary_);
    
    IJKSize paddedSize = sizeField.storage().paddedSize();
    int size = paddedSize.iSize() * paddedSize.jSize() * paddedSize.kSize();

    // allocate the external stroage
    std::vector<Real> memory1, memory2;
    memory1.resize(size, 33.0);
    memory2.resize(size, 77.0);

    Real *pDeviceMemory1, *pDeviceMemory2;
    cudaMalloc((void**)&pDeviceMemory1, size * sizeof(Real));
    cudaMemcpy(pDeviceMemory1, &memory1[0], size * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMalloc((void**)&pDeviceMemory2, size * sizeof(Real));
    cudaMemcpy(pDeviceMemory2, &memory2[0], size * sizeof(Real), cudaMemcpyHostToDevice);
    
    // define external storage
    ExternalStorage<Real> externalStorage1, externalStorage2;
    externalStorage1.Init(pDeviceMemory1, paddedSize);
    externalStorage2.Init(pDeviceMemory2, paddedSize);

    // normal data field initialization
    IJKRealFieldNoAlignment dataField;
    dataField.Init("dataField", calcDomainSize_, boundary_);

    // test setting the first storage
    dataField.SetExternalStorage(externalStorage1);
    ASSERT_EQ(33.0, dataField(0,0,0));
    dataField(0,0,0) = 44.0; // set the originls

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

    // make sure destructor is not accissing external storage which is deleted below
    dataField.SynchronizeDeviceStorage();
    cudaFree(pDeviceMemory1);
    cudaFree(pDeviceMemory2);
}