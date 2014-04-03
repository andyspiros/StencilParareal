#ifndef BOOST_PP_IS_ITERATING

    #ifndef WITH_WRAPPER_INCLUDED
    #define WITH_WRAPPER_INCLUDED

    #include <boost/mpl/assert.hpp>
    #include <boost/mpl/void.hpp>    
    #include <boost/mpl/integral_c.hpp>
    #include <boost/preprocessor/repetition.hpp>
    #include <boost/preprocessor/iteration/iterate.hpp>
    #include <boost/preprocessor/repetition/enum_params.hpp>
    #include <boost/preprocessor/cat.hpp>
    #include "FunctionParameter.h"
    #include "Definitions.h"
    #include "DummyContext.h"
    
    #include <boost/mpl/vector/vector20.hpp> // depends on MAX_PARAM_COUNT

    /**
    * @struct WithWrapper
    * Structure used to wrap stencil function in the definition of update functions.
    * It provides the With function which creates a function parameter using TFunc and the arguments handed over to the With function
    */
    template< 
        template<typename> class TStencilFunction,
        typename TDomain,
        typename TNumOfParameters>
    struct WithWrapper;

    // stencil function parameter definition wrong, couldn't count number of parameters 
    template< 
        template<typename> class TStencilFunction,
        typename TDomain>
    struct WithWrapper<TStencilFunction, TDomain, boost::mpl::void_> 
    { 
        BOOST_MPL_ASSERT_MSG(
            false, 
            CHECK_STENCIL_FUNCTION_PARAMETER_DEFINITION, 
            (TStencilFunction<StencilFunctionEnvironment<DummyContext, boost::mpl::vector0<> > >)
        );    
    };

    // generate specializations
    #define BOOST_PP_ITERATION_PARAMS_1 (3, (1, MAX_PARAM_COUNT, "WithWrapper.h"))
    #include BOOST_PP_ITERATE()

    #endif // WITH_WRAPPER_INCLUDED

#else // BOOST_PP_IS_ITERATING

    #define ITERATION_INDEX BOOST_PP_ITERATION()

    // specialization pattern
    template<
        template<typename> class TStencilFunction,
        typename TDomain> 
    struct WithWrapper<TStencilFunction, TDomain, boost::mpl::integral_c<int, ITERATION_INDEX> > 
    { 
        template<BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, typename T)> 
        __ACC__
        static FunctionParameter<
            TStencilFunction, 
            TDomain,
            Offset<0,0,0>, 
            BOOST_PP_CAT(boost::mpl::vector, ITERATION_INDEX)<BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, T)> > 
        With(BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, T)) 
        { 
            return FunctionParameter<
                TStencilFunction, 
                TDomain,
                Offset<0,0,0>, 
                BOOST_PP_CAT(boost::mpl::vector, ITERATION_INDEX)<BOOST_PP_ENUM_PARAMS(ITERATION_INDEX, T)> 
            >(); 
        } 
    };

    #undef ITERATION_INDEX

#endif // BOOST_PP_IS_ITERATING

  
