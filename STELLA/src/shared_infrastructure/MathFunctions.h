#pragma once

#include <cmath>
#include "Definitions.h"

#ifndef __CUDA_BACKEND__
#include <algorithm>
#endif

/**
* Namespace providing a set of math functions working for all backends
*/
namespace mf
{
    /**
    * fabs (see math.h)
    */
    template<typename T>
    __ACC_CPU__
    T fabs(const T x)
    {
        return ::fabs(x);
    }

    /**
    * floor (see math.h)
    */
    template<typename T>
    __ACC_CPU__
    T floor(const T x)
    {
        return ::floor(x);
    }

    /**
    * ceil (see math.h)
    */
    template<typename T>
    __ACC_CPU__
    T ceil(const T x)
    {
        return ::ceil(x);
    }

    /**
    * truncates x to an integer
    */
    template<typename T>
    __ACC_CPU__
    int trunc(const T x)
    {
         return static_cast<int>(x);
    }
    
    /**
    * fmod (see math.h)
    */
    template<typename T>
    __ACC_CPU__
    T fmod(const T x, const T y)
    {
        return ::fmod(x, y);
    }

    /**
    * pow (see math.h)
    *
    * x to the power of y
    */
    template<typename T>
    __ACC_CPU__
    T pow(const T x, const T y)
    {
        return ::pow(x, y);
    }
    
    /**
    * sqrt square root (see math.h)
    */
    template<typename T>
    __ACC_CPU__
    T sqrt(const T x)
    {
        return ::sqrt(x);
    }

    /**
    * Function computing the minimum given two inputs
    */
    template<typename T>
    __ACC_CPU__
    T min(const T x, const T y)
    {
#ifdef __CUDA_BACKEND__
        return x < y ? x : y;
#else
        return std::min(x, y);
#endif
    }

    /**
    * Function computing the maximum given two inputs
    */
    template<typename T>
    __ACC_CPU__
    T max(const T x, const T y)
    {
#ifdef __CUDA_BACKEND__
        return x > y ? x : y;
#else
        return std::max(x, y);
#endif
    }

    /**
    * Function computing the exponential
    */
    template<typename T>
    __ACC_CPU__
    T exp(const T x)
    {
        return ::exp(x);
    }

}
