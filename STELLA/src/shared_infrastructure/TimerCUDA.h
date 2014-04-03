#pragma once

#include "CUDADefinitions.h"
#include "CUDAStreams.h"
#include "Timer.h"

/**
* @class TimerCUDA
* CUDA implementation of the Timer interface
*/
class TimerCUDA : public Timer<TimerCUDA> // CRTP
{
public:
    __CPU__
    TimerCUDA() 
    { 
        // create the CUDA events
        cudaEventCreate(&start_); 
        cudaEventCreate(&stop_);
    }
    __CPU__
    ~TimerCUDA() 
    {
        // free the CUDA events 
        cudaEventDestroy(start_); 
        cudaEventDestroy(stop_); 
    }

    /**
    * Reset counters
    */
    __CPU__
    void ResetImpl() {}

    /**
    * Start the stop watch
    */
    __CPU__
    void StartImpl()
    {
        // insert a start event
        cudaEventRecord(start_, 0);
    }

    /**
    * Pause the stop watch
    */
    __CPU__
    double PauseImpl()
    {
        // insert stop event and wait for it
        cudaEventRecord(stop_, 0); 
        cudaEventSynchronize(stop_);
        
        // compute the timing
        float result;
        cudaEventElapsedTime(&result, start_, stop_);
        return static_cast<double>(result) * 0.001f; // convert ms to s
    }

private:
    cudaEvent_t start_;
    cudaEvent_t stop_;
};

  

