#pragma once

#include <cstdlib>
#include <iostream>
#include <cuda_runtime_api.h>

/**
* Method asserting there is no CUDA error. 
* @descriptor descriptor used to generate an error message
*/
inline void assertNoCUDAError(const char* descriptor) 
{
    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess)
    {
        std::cerr << descriptor << " -> CUDA call failed with " << cudaGetErrorString(error) << std::endl;
        assert(false);
        exit(-1);
    }
    // useful for debugging purposes --> error shows up where it happens
    //error = cudaThreadSynchronize();
    //if(error != cudaSuccess)
    //{
    //    std::cerr << descriptor << " -> CUDA synchronization failed with " << cudaGetErrorString(error) << std::endl;
    //    assert(false);
    //    exit(-1);
    //}
}