#include "CUDAStreams.h"

// initialize the reference counter to 0
int CUDAStreams::referenceCounter_ = 0;
// initialize stream structures with the default stream value 0
cudaStream_t CUDAStreams::kernelStream_ = 0;
cudaStream_t CUDAStreams::copyStream_ = 0;

void CUDAStreams::IncrementReferenceCounter()
{
    assert(referenceCounter_ >= 0);
    ++referenceCounter_;

    // if called for the first time create the streams
    if(referenceCounter_ == 1)
    {
#ifdef ENABLE_CUDA_STREAMS
        cudaStreamCreate(&kernelStream_);
        cudaStreamCreate(&copyStream_);
        assertNoCUDAError("CUDAStreams");
#endif
    }
}

void CUDAStreams::DecrementReferenceCounter()
{
    assert(referenceCounter_ > 0);
    --referenceCounter_;

    // if we are back at zero destroy the streams
    if(referenceCounter_ == 0)
    {
#ifdef ENABLE_CUDA_STREAMS
        cudaStreamDestroy(kernelStream_);
        cudaStreamDestroy(copyStream_);
        assertNoCUDAError("CUDAStreams");
#endif
    }
}

cudaStream_t& CUDAStreams::kernelStream()
{
    assert(referenceCounter_ > 0);
    return kernelStream_;
}

cudaStream_t& CUDAStreams::copyStream()
{
    assert(referenceCounter_ > 0);
    return copyStream_;
}
