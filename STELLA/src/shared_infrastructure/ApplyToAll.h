#pragma once

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/next_prior.hpp>
#include <boost/mpl/deref.hpp>
#include "Definitions.h"
#include "ParameterType.h"

#define APPLY_TO_ALL_QUALIFIER __CPU__      
#include "ApplyToAllDefinition.h"
#undef APPLY_TO_ALL_QUALIFIER

// provide an accelerator specific implementation 
namespace acc 
{
    #define APPLY_TO_ALL_QUALIFIER __ACC__    
    #include "ApplyToAllDefinition.h"
    #undef APPLY_TO_ALL_QUALIFIER
}


  
