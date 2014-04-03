#include <vector>
#include "gtest/gtest.h"
#include "DataFieldCUDAStorage.h"

typedef DataFieldAlignment<cDimK, 1> Alignment;
typedef DataFieldIJBoundary<-1, 2, -3, 4> Boundary;

// test DataFieldStorage and Iterators
class DataFieldCUDAStorageUnittest : public ::testing::Test
{
protected:
    // test data
    typedef double DataType3D;
    typedef float DataType2D;
    typedef int DataType1D;

    typedef DataFieldStorageFormat<Boundary, StorageOrder::IJK, Alignment> IJKStorageFormat;
    typedef DataFieldStorageFormat<Boundary, StorageOrder::IJ, Alignment> IJStorageFormat;
    typedef DataFieldStorageFormat<Boundary, StorageOrder::I, Alignment> IStorageFormat;

    DataFieldCUDAStorage<DataType3D, IJKStorageFormat> storage3D_;
    DataFieldCUDAStorage<DataType2D, IJStorageFormat> storage2D_;
    DataFieldCUDAStorage<DataType1D, IStorageFormat> storage1D_;
 
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
};

TEST_F(DataFieldCUDAStorageUnittest, Copy)
{
    // initialize data field with random numbers
    DataType3D value1, value2, value3;
    value1 = 3.0;
    cudaMemcpy(storage3D_.pStorageOrigin(), &value1, sizeof(DataType3D), cudaMemcpyHostToDevice);
    value1 = 2.0;
    cudaMemcpy(storage2D_.pStorageOrigin(), &value1, sizeof(DataType3D), cudaMemcpyHostToDevice);
    value1 = 1.0;
    cudaMemcpy(storage1D_.pStorageOrigin(), &value1, sizeof(DataType3D), cudaMemcpyHostToDevice);
    
    // copy fields
    DataFieldCUDAStorage<DataType3D, IJKStorageFormat> storage3DCopyConstruct(storage3D_);
    DataFieldCUDAStorage<DataType2D, IJStorageFormat> storage2DCopyConstruct(storage2D_);
    DataFieldCUDAStorage<DataType1D, IStorageFormat> storage1DCopyConstruct(storage1D_);

    DataFieldCUDAStorage<DataType3D, IJKStorageFormat> storage3DAssign;
    DataFieldCUDAStorage<DataType2D, IJStorageFormat> storage2DAssign;
    DataFieldCUDAStorage<DataType1D, IStorageFormat> storage1DAssign;

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

    cudaMemcpy(&value1, storage3D_.pStorageOrigin(), sizeof(DataType3D), cudaMemcpyDeviceToHost);
    cudaMemcpy(&value2, storage3DCopyConstruct.pStorageOrigin(), sizeof(DataType3D), cudaMemcpyDeviceToHost);
    cudaMemcpy(&value3, storage3DAssign.pStorageOrigin(), sizeof(DataType3D), cudaMemcpyDeviceToHost);
    ASSERT_EQ(value1, value2);
    ASSERT_EQ(value1, value3);

    cudaMemcpy(&value1, storage2D_.pStorageOrigin(), sizeof(DataType3D), cudaMemcpyDeviceToHost);
    cudaMemcpy(&value2, storage2DCopyConstruct.pStorageOrigin(), sizeof(DataType3D), cudaMemcpyDeviceToHost);
    cudaMemcpy(&value3, storage2DAssign.pStorageOrigin(), sizeof(DataType3D), cudaMemcpyDeviceToHost);
    ASSERT_EQ(value1, value2);
    ASSERT_EQ(value1, value3);

    cudaMemcpy(&value1, storage1D_.pStorageOrigin(), sizeof(DataType3D), cudaMemcpyDeviceToHost);
    cudaMemcpy(&value2, storage1DCopyConstruct.pStorageOrigin(), sizeof(DataType3D), cudaMemcpyDeviceToHost);
    cudaMemcpy(&value3, storage1DAssign.pStorageOrigin(), sizeof(DataType3D), cudaMemcpyDeviceToHost);
    ASSERT_EQ(value1, value2);
    ASSERT_EQ(value1, value3);
}

TEST_F(DataFieldCUDAStorageUnittest, ExternalStorage)
{
    // setup external memory
    IJKSize paddedSize = storage3D_.paddedSize();
    int size = paddedSize.iSize() * paddedSize.jSize() * paddedSize.kSize();
    DataType3D value1, value2;
    
    std::vector<DataType3D> memory1, memory2;
    memory1.resize(size, 11.0);
    memory2.resize(size, 22.0);
    
    DataType3D *pDeviceMemory1, *pDeviceMemory2;
    cudaMalloc((void**)&pDeviceMemory1, size * sizeof(DataType3D));
    cudaMemcpy(pDeviceMemory1, &memory1[0], sizeof(DataType3D) * size, cudaMemcpyHostToDevice);
    cudaMalloc((void**)&pDeviceMemory2, size * sizeof(DataType3D));
    cudaMemcpy(pDeviceMemory2, &memory2[0], sizeof(DataType3D) * size, cudaMemcpyHostToDevice);

    ExternalStorage<DataType3D> externalStorage1, externalStorage2;
    externalStorage1.Init(pDeviceMemory1, paddedSize);
    externalStorage2.Init(pDeviceMemory2, paddedSize);

    // storage with external memory
    DataFieldCUDAStorage<DataType3D, IJKStorageFormat> storage;
    storage.Init(externalStorage1, calculationDomain_, kBoundary_);
    
    // check the memory is ok
    ASSERT_EQ(pDeviceMemory1, storage.pStorageBase());
    cudaMemcpy(&value1, pDeviceMemory1, sizeof(DataType3D), cudaMemcpyDeviceToHost);
    ASSERT_EQ(11.0, value1);

    // set the second external storage
    storage.SetExternalStorage(externalStorage2);

    // check the memory is ok
    ASSERT_EQ(pDeviceMemory2, storage.pStorageBase());
    cudaMemcpy(&value2, pDeviceMemory2, sizeof(DataType3D), cudaMemcpyDeviceToHost);
    ASSERT_EQ(22.0, value2);

    cudaFree(pDeviceMemory1);
    cudaFree(pDeviceMemory2);
}

