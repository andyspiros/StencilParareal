#pragma once

#include "Timer.h"

/**
* @class TimerDummy
* Dummy timer implementation doing nothing in order to avoid runtime overhead
*/
class TimerDummy : public Timer<TimerDummy> // CRTP
{
public:
    TimerDummy() {}
    ~TimerDummy() {}

    /**
    * Reset counters
    */
    void ResetImpl() {}

    /**
    * Start the stop watch
    */
    void StartImpl() {}

    /**
    * Pause the stop watch
    */
    double PauseImpl() { return 0.0; }
};

  
