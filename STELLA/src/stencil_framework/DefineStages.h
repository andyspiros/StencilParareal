#ifndef BOOST_PP_IS_ITERATING

    #ifndef MERGE_SUB_STAGES_INCLUDED
    #define MERGE_SUB_STAGES_INCLUDED

    #include <boost/preprocessor/repetition.hpp>
    #include <boost/preprocessor/iteration/iterate.hpp>
    #include <boost/preprocessor/cat.hpp>
    #include "Definitions.h"
    #include "StencilStage.h"

    #include <boost/mpl/vector/vector10.hpp> // depends on MAX_MERGE_COUNT

    #define TEXT_STENCIL_STAGE(z, n, data) \
        StencilStage< \
            BOOST_PP_CAT(TS, n), \
            BOOST_PP_CAT(TIJR, n), \
            BOOST_PP_CAT(TKR, n) \
        >

    // generate specializations
    #define BOOST_PP_ITERATION_PARAMS_1 (3, (1, MAX_MERGE_COUNT, "DefineStages.h"))
    #include BOOST_PP_ITERATE()

    #undef TEXT_STENCIL_STAGE

    #endif // MERGE_SUB_STAGES_INCLUDED

#else // BOOST_PP_IS_ITERATING

    #define ITERATION_INDEX BOOST_PP_ITERATION()

    // specialization pattern
    template<
        BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, template<typename> class TS),
        BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, typename TIJR),
        BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, typename TKR)>
    __CPU__
    BOOST_PP_CAT(boost::mpl::vector, ITERATION_INDEX)<BOOST_PP_ENUM(ITERATION_INDEX, TEXT_STENCIL_STAGE, void)>
    define_stages(BOOST_PP_ENUM(ITERATION_INDEX, TEXT_STENCIL_STAGE, void))
    {
        BOOST_PP_CAT(boost::mpl::vector, ITERATION_INDEX)<BOOST_PP_ENUM(ITERATION_INDEX, TEXT_STENCIL_STAGE, void)> result;
        return result;
    };

    #undef ITERATION_INDEX

#endif // BOOST_PP_IS_ITERATING
  
