#ifndef BOOST_PP_IS_ITERATING

    #ifndef PACK_PARAMETERS_INCLUDED
    #define PACK_PARAMETERS_INCLUDED

    #include <boost/preprocessor/repetition.hpp>
    #include <boost/preprocessor/cat.hpp>
    #include <boost/preprocessor/iteration/iterate.hpp>
    #include "Enums.h"   
    #include "ParameterWrapper.h"
    #include "Tuple.h"

    #define TEXT_PARAMETER_TYPE(z, n, data) \
        ParameterWrapper< \
            BOOST_PP_CAT(TV, n), \
            BOOST_PP_CAT(VI, n), \
            BOOST_PP_CAT(VP, n) \
        >

    #define TEXT_PARAMETER_DEFINTION(z, n, data) \
        ParameterWrapper< \
            BOOST_PP_CAT(TV, n), \
            BOOST_PP_CAT(VI, n), \
            BOOST_PP_CAT(VP, n) \
        > BOOST_PP_CAT(val, n)

    // generate specializations
    #define BOOST_PP_ITERATION_PARAMS_1 (3, (1, MAX_TUPLE_SIZE, "PackParameters.h"))
    #include BOOST_PP_ITERATE()

    #undef TEXT_PARAMETER_TYPE
    #undef TEXT_PARAMETER_DEFINTION

    #endif // PACK_PARAMETERS_INCLUDED

#else // BOOST_PP_IS_ITERATING

    // define the iteration index
    #define ITERATION_INDEX BOOST_PP_ITERATION()

    template<
        BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, typename TV),
        BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, ParameterIntend VI),
        BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, int VP)>
    __CPU__
    Tuple< 
        TupleElements<
            BOOST_PP_CAT(BOOST_PP_CAT(boost::mpl::vector, ITERATION_INDEX), _c)<int, BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, VP)>,    
            BOOST_PP_CAT(boost::mpl::vector, ITERATION_INDEX)<BOOST_PP_ENUM(ITERATION_INDEX, TEXT_PARAMETER_TYPE, void)>            
        >
    > 
    pack_parameters(BOOST_PP_ENUM(ITERATION_INDEX, TEXT_PARAMETER_DEFINTION, void))
    {        
        // define the result tuple
        Tuple<
            TupleElements<
                BOOST_PP_CAT(BOOST_PP_CAT(boost::mpl::vector, ITERATION_INDEX), _c)<int, BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, VP)>,    
                BOOST_PP_CAT(boost::mpl::vector, ITERATION_INDEX)<BOOST_PP_ENUM(ITERATION_INDEX, TEXT_PARAMETER_TYPE, void)>                
            >
        > result;

        // initialize the result tuple
        result.Init(BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, val));

        return result;
    };

#endif // BOOST_PP_IS_ITERATING
