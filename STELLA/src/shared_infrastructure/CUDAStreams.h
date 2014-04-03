#pragma once

#include "Definitions.h"
#include "CUDADefinitions.h"

/**
* @class CUDAStreams
* Static class providing a kernel stream and a copy stream reference
* (Note that reference counting is used for stream creation and destructions as static construction of CUDA objects is not safe)
*/
class CUDAStreams
{
    DISALLOW_COPY_AND_ASSIGN(CUDAStreams);
public:

    /**
    * Method incrementing the reference counter, if it is called for the first time the streams are initialized
    */
    __CPU__
    static void IncrementReferenceCounter();

    /**
    * Method decrementing the reference counter, if the counter is decremented to 0 the streams are deleted
    */
    __CPU__
    static void DecrementReferenceCounter();

    /**
    * @return kernel execution stream
    */
    __CPU__
    static cudaStream_t& kernelStream();

    /**
    * @return copy stream
    */
    __CPU__
    static cudaStream_t& copyStream();

private:
    static cudaStream_t kernelStream_;
    static cudaStream_t copyStream_;
    static int referenceCounter_;
};

  

