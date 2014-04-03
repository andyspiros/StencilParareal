#ifndef BOOST_PP_IS_ITERATING

    #ifndef DEFINE_CASE_INCLUDED
    #define DEFINE_CASE_INCLUDED

    #include <boost/preprocessor/repetition.hpp>
    #include <boost/preprocessor/iteration/iterate.hpp>
    #include <boost/preprocessor/cat.hpp>
    #include "Enums.h"
    #include "Definitions.h"
    #include "StencilSweepDescriptor.h"
    #include "StencilSweepGroupDescriptor.h"

    #include <boost/mpl/vector/vector30.hpp> // depends on MAX_GROUP_COUNT

    #define TEXT_STENCIL_SWEEP(z, n, data) \
        StencilSweepDescriptor< \
            BOOST_PP_CAT(TC, n), \
            BOOST_PP_CAT(TS, n), \
            BOOST_PP_CAT(VD, n) \
        >

    // generate specializations
    #define BOOST_PP_ITERATION_PARAMS_1 (3, (1, MAX_CASE_COUNT, "DefineCase.h"))
    #include BOOST_PP_ITERATE()

    #undef TEXT_STENCIL_SWEEP

    #endif // DEFINE_CASE_INCLUDED

#else // BOOST_PP_IS_ITERATING

    #define ITERATION_INDEX BOOST_PP_ITERATION()

    // specialization pattern
    template<
        typename TValue,
        TValue VValue,
        BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, typename TC),
        BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, typename TS),
        BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, KLoopDirection VD)>
    __CPU__
    StencilSweepGroupDescriptor<
        BOOST_PP_CAT(boost::mpl::vector, ITERATION_INDEX)<BOOST_PP_ENUM(ITERATION_INDEX, TEXT_STENCIL_SWEEP, void)>,
        ComparisonValue<TValue, VValue>,
        0 
    >
    define_case(BOOST_PP_ENUM(ITERATION_INDEX, TEXT_STENCIL_SWEEP, void))
    {
        StencilSweepGroupDescriptor<
            BOOST_PP_CAT(boost::mpl::vector, ITERATION_INDEX)<BOOST_PP_ENUM(ITERATION_INDEX, TEXT_STENCIL_SWEEP, void)>,
            ComparisonValue<TValue, VValue>,
            0
        > result;
        return result;
    };

    #undef ITERATION_INDEX

#endif // BOOST_PP_IS_ITERATING
  
