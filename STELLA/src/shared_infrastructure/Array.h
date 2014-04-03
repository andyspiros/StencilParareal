#ifndef BOOST_PP_IS_ITERATING

    #ifndef ARRAY_INCLUDED
    #define ARRAY_INCLUDED

    #include <boost/type_traits/remove_const.hpp>
    #include <boost/mpl/void.hpp>
    #include <boost/mpl/bool.hpp>
    #include <boost/preprocessor/repetition.hpp>
    #include <boost/preprocessor/enum.hpp>
    #include <boost/preprocessor/iteration/iterate.hpp>
    #include <boost/preprocessor/cat.hpp>
    #include "Definitions.h"

    #include <boost/mpl/vector/vector10_c.hpp> // depends on MAX_ARRAY_SIZE

    // define the Array implementation class
    #define DEFINE_ARRAY template<typename TValue, int VSize> class Array;

    // define an is array trait
    #define DEFINE_IS_ARRAY template<typename T> struct is_array : boost::mpl::false_ {};
    #define DEFINE_IS_ARRAY_SPECIALIZATION template<typename TValue, int VSize> struct is_array<Array<TValue, VSize> > : boost::mpl::true_ {};

    DEFINE_ARRAY
    DEFINE_IS_ARRAY
    DEFINE_IS_ARRAY_SPECIALIZATION

    // provide an accelerator specific array
    namespace acc 
    {
        DEFINE_ARRAY
        DEFINE_IS_ARRAY
        DEFINE_IS_ARRAY_SPECIALIZATION
    }

    #undef DEFINE_ARRAY
    #undef DEFINE_IS_ARRAY
    #undef DEFINE_IS_ARRAY_SPECIALIZATION

    // generate specializations
    #define BOOST_PP_ITERATION_PARAMS_1 (3, (0, MAX_ARRAY_SIZE, "Array.h"))
    #include BOOST_PP_ITERATE()

    #endif // ARRAY_INCLUDED
#else // BOOST_PP_IS_ITERATING
    #define ITERATION_INDEX BOOST_PP_ITERATION()

    #define TEXT_NUMBER(z, n, data) n

    #define TEXT_ASSIGN(z, n, data) \
        BOOST_PP_CAT(BOOST_PP_CAT(val, n), _) = BOOST_PP_CAT(BOOST_PP_CAT(other.val, n), _);

    #define TEXT_AT(z, n, data) \
        ARRAY_QUALIFIER \
        ValueType & At(boost::mpl::integral_c<int, n >*) { return BOOST_PP_CAT(BOOST_PP_CAT(val, n), _); }

    #define TEXT_CONST_AT(z, n, data) \
        ARRAY_QUALIFIER \
        typename boost::add_const<ValueType>::type & At(boost::mpl::integral_c<int, n >*) const { return BOOST_PP_CAT(BOOST_PP_CAT(val, n), _); }

    #define TEXT_MEMBER(z, n, data) \
        ValueType BOOST_PP_CAT(BOOST_PP_CAT(val,n),_);

    // define the standard array implementation
    #define ARRAY_QUALIFIER __CPU__     
    #include "ArrayDefinition.h"
    #undef ARRAY_QUALIFIER

    namespace acc
    {
        #define ARRAY_QUALIFIER __ACC__
        #include "ArrayDefinition.h"
        #undef ARRAY_QUALIFIER
    }

    #undef TEXT_NUMBER
    #undef TEXT_ASSIGN
    #undef TEXT_AT
    #undef TEXT_CONST_AT
    #undef TEXT_MEMBER

    #undef ITERATION_INDEX
#endif // BOOST_PP_IS_ITERATING

