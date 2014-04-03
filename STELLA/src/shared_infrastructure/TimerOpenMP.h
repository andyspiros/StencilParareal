#pragma once

#include "Timer.h"

#ifdef __OPENMP_BACKEND__
#include <omp.h>
#else
#include <ctime>
#endif

/**
* @class TimerOpenMP
* OpenMP implementation of the Timer interface
*/
class TimerOpenMP : public Timer<TimerOpenMP> // CRTP
{
public:
    TimerOpenMP() { startTime_ = 0.0; }
    ~TimerOpenMP() {}

    /**
    * Reset counters
    */
    void ResetImpl()
    {
        startTime_ = 0.0;
    }

    /**
    * Start the stop watch
    */
    void StartImpl()
    {
#ifdef __OPENMP_BACKEND__
        startTime_ = omp_get_wtime();
#else
        startTime_ = (double)clock() / (double)CLOCKS_PER_SEC;
#endif
    }

    /**
    * Pause the stop watch
    */
    double PauseImpl()
    {
#ifdef __OPENMP_BACKEND__
        return omp_get_wtime() - startTime_;
#else
        return ((double)clock() / (double)CLOCKS_PER_SEC) - startTime_;
#endif
    }

private:
    double startTime_;
};

  
