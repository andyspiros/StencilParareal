#pragma once

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/next_prior.hpp>
#include <boost/mpl/deref.hpp>
#include "Definitions.h"
#include "Array.h"
#include "ParameterType.h"

#define MODIFY_ARRAY_QUALIFIER __CPU__     
#include "ArrayAlgorithmsDefinition.h"
#undef MODIFY_ARRAY_QUALIFIER

// provide an accelerator specific implementation 
namespace acc 
{
    #define MODIFY_ARRAY_QUALIFIER __ACC__    
    #include "ArrayAlgorithmsDefinition.h"
    #undef MODIFY_ARRAY_QUALIFIER
}
  
