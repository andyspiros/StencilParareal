#include "StencilFramework.h"
#include <cmath>

template<typename TEnv>
struct Advection2
{  
    STENCIL_FUNCTION(TEnv)

    FUNCTION_PARAMETER(0, data)
    FUNCTION_PARAMETER(1, _dx)
    FUNCTION_PARAMETER(2, _cx)
    FUNCTION_PARAMETER(3, _cy)
    FUNCTION_PARAMETER(4, _cz)

    __ACC__
    static T Do(Context ctx)
    {
#ifdef UPWINDADVECTION
        // Direction i
        T dq_di = 1.5 * ctx[data::Center()];
        if (ctx[_cx::Center()] > 0.)
            dq_di = + 1.5 * ctx[data::Center()]
                    -  2. * ctx[data::At(Offset<-1, 0, 0>())]
                    +  .5 * ctx[data::At(Offset<-2, 0, 0>())];
        else
            dq_di = - 1.5 * ctx[data::Center()]
                    +  2. * ctx[data::At(Offset< 1, 0, 0>())]
                    -  .5 * ctx[data::At(Offset< 2, 0, 0>())];
        dq_di *= ctx[_cx::Center()];

        // Direction j
        T dq_dj = 1.5 * ctx[data::Center()];
        if (ctx[_cy::Center()] > 0.)
            dq_dj = + 1.5 * ctx[data::Center()]
                    -  2. * ctx[data::At(Offset<0, -1, 0>())]
                    +  .5 * ctx[data::At(Offset<0, -2, 0>())];
        else                                           
            dq_dj = - 1.5 * ctx[data::Center()]        
                    +  2. * ctx[data::At(Offset<0,  1, 0>())]
                    -  .5 * ctx[data::At(Offset<0,  2, 0>())];
        dq_dj *= ctx[_cy::Center()];

        // Direction k
        T dq_dk = 1.5 * ctx[data::Center()];        
        if (ctx[_cz::Center()] > 0.)
            dq_dk =   1.5 * ctx[data::Center()]
                    -  2. * ctx[data::At(Offset<0, 0, -1>())]
                    +  .5 * ctx[data::At(Offset<0, 0, -2>())];
        else
            dq_dk = - 1.5 * ctx[data::Center()]
                    +  2. * ctx[data::At(Offset<0, 0,  1>())]
                    -  .5 * ctx[data::At(Offset<0, 0,  2>())];
        dq_dk *= ctx[_cz::Center()];

        T tensadv = -(dq_di + dq_dj + dq_dk) / ctx[_dx::Center()];
#else
        // Compute advection
        T dq_di = 
            +1. * ctx[data::At(Offset<-1, 0, 0>())]
            -1. * ctx[data::At(Offset< 1, 0, 0>())];
        T dq_dj = 
            +1. * ctx[data::At(Offset< 0, -1, 0>())]
            -1. * ctx[data::At(Offset< 0,  1, 0>())];
        T dq_dk = 
            +1. * ctx[data::At(Offset< 0, 0, -1>())]
            -1. * ctx[data::At(Offset< 0, 0,  1>())];
        T tensadv = 
            (
               ctx[_cx::Center()] * dq_di
             + ctx[_cy::Center()] * dq_dj 
             + ctx[_cz::Center()] * dq_dk
            ) / (2. * ctx[_dx::Center()]);
#endif

        return tensadv;
    }
};

template<typename TEnv>
struct Advection4
{  
    STENCIL_FUNCTION(TEnv)

    FUNCTION_PARAMETER(0, data)
    FUNCTION_PARAMETER(1, _dx)
    FUNCTION_PARAMETER(2, _cx)
    FUNCTION_PARAMETER(3, _cy)
    FUNCTION_PARAMETER(4, _cz)

    __ACC__
    static T Do(Context ctx)
    {
        // Compute advection
        T dq_di = 
            -1. * ctx[data::At(Offset<-2, 0, 0>())]
            +8. * ctx[data::At(Offset<-1, 0, 0>())]
            -8. * ctx[data::At(Offset< 1, 0, 0>())]
            +1. * ctx[data::At(Offset< 2, 0, 0>())];
        T dq_dj = 
            -1. * ctx[data::At(Offset< 0, -2, 0>())]
            +8. * ctx[data::At(Offset< 0, -1, 0>())]
            -8. * ctx[data::At(Offset< 0,  1, 0>())]
            +1. * ctx[data::At(Offset< 0,  2, 0>())];
        T dq_dk = 
            -1. * ctx[data::At(Offset< 0, 0, -2>())]
            +8. * ctx[data::At(Offset< 0, 0, -1>())]
            -8. * ctx[data::At(Offset< 0, 0,  1>())]
            +1. * ctx[data::At(Offset< 0, 0,  2>())];
        T tensadv = 
            (
               ctx[_cx::Center()] * dq_di
             + ctx[_cy::Center()] * dq_dj 
             + ctx[_cz::Center()] * dq_dk
            ) / (12. * ctx[_dx::Center()]);
        return tensadv;
    }
};

template<typename TEnv>
struct Laplace2
{  
    STENCIL_FUNCTION(TEnv)

    FUNCTION_PARAMETER(0, data)
    FUNCTION_PARAMETER(1, _dx2)

    __ACC__
    static T Do(Context ctx)
    {
        // Compute laplacian
        T d2q_di2 =
            +1. * ctx[data::At(Offset<-1, 0, 0>())]
            -2. * ctx[data::At(Offset< 0, 0, 0>())]
            +1. * ctx[data::At(Offset< 1, 0, 0>())];
        T d2q_dj2 =
            +1. * ctx[data::At(Offset< 0, -1, 0>())]
            -2. * ctx[data::At(Offset< 0,  0, 0>())]
            +1. * ctx[data::At(Offset< 0,  1, 0>())];
        T d2q_dk2 =
            +1. * ctx[data::At(Offset< 0, 0, -1>())]
            -2. * ctx[data::At(Offset< 0, 0,  0>())]
            +1. * ctx[data::At(Offset< 0, 0,  1>())];
        return (d2q_di2 + d2q_dj2 + d2q_dk2) / (ctx[_dx2::Center()]);
    }
};


template<typename TEnv>
struct Laplace4
{  
    STENCIL_FUNCTION(TEnv)

    FUNCTION_PARAMETER(0, data)
    FUNCTION_PARAMETER(1, _dx2)

    __ACC__
    static T Do(Context ctx)
    {
        // Compute laplacian
        T d2q_di2 =
            - 1. * ctx[data::At(Offset<-2, 0, 0>())]
            +16. * ctx[data::At(Offset<-1, 0, 0>())]
            -30. * ctx[data::At(Offset< 0, 0, 0>())]
            +16. * ctx[data::At(Offset< 1, 0, 0>())]
            - 1. * ctx[data::At(Offset< 2, 0, 0>())];
        T d2q_dj2 =
            - 1. * ctx[data::At(Offset< 0, -2, 0>())]
            +16. * ctx[data::At(Offset< 0, -1, 0>())]
            -30. * ctx[data::At(Offset< 0,  0, 0>())]
            +16. * ctx[data::At(Offset< 0,  1, 0>())]
            - 1. * ctx[data::At(Offset< 0,  2, 0>())];
        T d2q_dk2 =
            - 1. * ctx[data::At(Offset< 0, 0, -2>())]
            +16. * ctx[data::At(Offset< 0, 0, -1>())]
            -30. * ctx[data::At(Offset< 0, 0,  0>())]
            +16. * ctx[data::At(Offset< 0, 0,  1>())]
            - 1. * ctx[data::At(Offset< 0, 0,  2>())];
        return (d2q_di2 + d2q_dj2 + d2q_dk2) / (12. * ctx[_dx2::Center()]);
    }
};


