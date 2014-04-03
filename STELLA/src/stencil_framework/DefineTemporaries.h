#ifndef BOOST_PP_IS_ITERATING

    #ifndef DEFINE_TEMPORARIES_INCLUDED
    #define DEFINE_TEMPORARIES_INCLUDED

    #include <boost/preprocessor/repetition.hpp>
    #include <boost/preprocessor/iteration/iterate.hpp>
    #include <boost/preprocessor/cat.hpp>
    #include "Definitions.h"

    #include <boost/mpl/vector/vector20.hpp> // depends on MAX_TEMP_COUNT

    // generate specializations
    #define BOOST_PP_ITERATION_PARAMS_1 (3, (1, MAX_TEMP_COUNT, "DefineTemporaries.h"))
    #include BOOST_PP_ITERATE()

    #endif // DEFINE_TEMPORARIES_INCLUDED

#else // BOOST_PP_IS_ITERATING

    #define ITERATION_INDEX BOOST_PP_ITERATION()

    // specialization pattern
    template<BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, typename T)>
    __CPU__
    BOOST_PP_CAT(boost::mpl::vector, ITERATION_INDEX)<BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, T)>
    define_temporaries(BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, T))
    {
        BOOST_PP_CAT(boost::mpl::vector, ITERATION_INDEX)<BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, T)> result;
        return result;
    };

    #undef ITERATION_INDEX

#endif // BOOST_PP_IS_ITERATING
  
