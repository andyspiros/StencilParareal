#pragma once

#include <sstream>
#include <string>
#include "Definitions.h"

/**
* @class UsageMeter
* Measures the usage of a certain class, for instance number of apply calls or bytes transfered
*/
class UsageMeter 
{
    DISALLOW_COPY_AND_ASSIGN(UsageMeter);
public:
    __CPU__
    UsageMeter()
    {
        name_ = "";
        Reset();
    }
    __CPU__
    ~UsageMeter() {}

    /**
    * Initialize name and unit
    * @param name of meter for print out
    * @param unit unit of the usage counter
    */
    __CPU__
    void Init(std::string name, std::string unit)
    {
        name_ = name;
        unit_ = unit;
        Reset();
    }

    /**
    * Reset counters
    */
    __CPU__
    void Reset()
    {
        counter_ = 0.0;
    }

    /**
    * Increment the usage counter
    */
    __CPU__
    void Add(double increment)
    {
        counter_ += increment;
    }

    /**
    * @return usage counter
    */
    __CPU__
    double counter() const
    {
        return counter_;
    }

    /**
    * @return usage counter as string
    */
    __CPU__
    std::string ToString() const
    {
        std::ostringstream out; 
        out << name_ << "\t[" << unit_ << "]\t" << counter_; 
        return out.str();
    }

private:
    std::string name_;
    std::string unit_;
    double counter_;
};

  
