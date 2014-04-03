#ifndef BOOST_PP_IS_ITERATING

    #ifndef DEFINE_SWITCH_INCLUDED
    #define DEFINE_SWITCH_INCLUDED

    #include <boost/preprocessor/repetition.hpp>
    #include <boost/preprocessor/iteration/iterate.hpp>
    #include <boost/preprocessor/repetition/enum.hpp>
    #include <boost/preprocessor/cat.hpp>
    #include <boost/mpl/pair.hpp>
    #include "Definitions.h"
    #include "StencilSweepGroupDescriptor.h"

    #include <boost/mpl/vector/vector30.hpp> // depends on MAX_GROUP_COUNT

    #define TEXT_STENCIL_SWEEP_GROUP(z, n, data) \
        StencilSweepGroupDescriptor< \
            BOOST_PP_CAT(TS, n), \
            BOOST_PP_CAT(TC, n), \
            data \
        >

    // generate specializations
    #define BOOST_PP_ITERATION_PARAMS_1 (3, (1, MAX_SWITCH_COUNT, "DefineSwitch.h"))
    #include BOOST_PP_ITERATE()

    #undef TEXT_STENCIL_SWEEP_GROUP

    #endif // DEFINE_SWITCH_INCLUDED

#else // BOOST_PP_IS_ITERATING

    #define ITERATION_INDEX BOOST_PP_ITERATION()

    // specialization pattern
    template<
        int VParameterIndex,
        BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, typename TS),
        BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, typename TC)>
    __CPU__
    BOOST_PP_CAT(boost::mpl::vector, ITERATION_INDEX)<
        BOOST_PP_ENUM(ITERATION_INDEX, TEXT_STENCIL_SWEEP_GROUP, VParameterIndex)
    >
    define_switch(BOOST_PP_ENUM(ITERATION_INDEX, TEXT_STENCIL_SWEEP_GROUP, 0))
    {
        // define a vector of stencil sweep groups from the given sweep group parameters
        // (note that the parameter index of the sweep groups is updated with the switch parameter index)
        BOOST_PP_CAT(boost::mpl::vector, ITERATION_INDEX)<
            BOOST_PP_ENUM(ITERATION_INDEX, TEXT_STENCIL_SWEEP_GROUP, VParameterIndex)
        > result;
        return result;
    };

    #undef ITERATION_INDEX

#endif // BOOST_PP_IS_ITERATING
  
