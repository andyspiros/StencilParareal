#pragma once

#include <sstream>
#include <string>
#include "Definitions.h"

/**
* @class Timer
* Measures total elapsed time between all start and stop calls
*/
template<typename TTimerImpl>
class Timer 
{
    DISALLOW_COPY_AND_ASSIGN(Timer);
protected:
    __CPU__
    Timer()
    {
        name_ = "";
        Reset();
    }
    __CPU__
    ~Timer() {}

public:
    /**
    * Initialize name and Reset()
    * @param name of meter for print out
    */
    __CPU__
    void Init(std::string name)
    {
        name_ = name;
        Reset();
    }

    /**
    * Reset counters
    */
    __CPU__
    void Reset()
    {
        totalTime_ = 0.0;
        static_cast<TTimerImpl*>(this)->ResetImpl();
    }

    /**
    * Start the stop watch
    */
    __CPU__
    void Start()
    {
        static_cast<TTimerImpl*>(this)->StartImpl();
    }

    /**
    * Pause the stop watch
    */
    __CPU__
    void Pause()
    {
        totalTime_ += static_cast<TTimerImpl*>(this)->PauseImpl();
    }
    
    /**
    * @return total elapsed time [s]
    */
    __CPU__
    double totalTime() const
    {
        return totalTime_;
    }

    /**
    * @return total elapsed time [s] as string
    */
    __CPU__
    std::string ToString() const
    {
        std::ostringstream out; 
        out << name_ << "\t[s]\t" << totalTime_; 
        return out.str();
    }

private:
    std::string name_;
    double totalTime_;
};

  
