#ifndef BOOST_PP_IS_ITERATING

    #ifndef TUPLE_INCLUDED
    #define TUPLE_INCLUDED

    #include <boost/type_traits/add_const.hpp>
    #include <boost/mpl/bool.hpp>
    #include <boost/mpl/integral_c.hpp>
    #include <boost/mpl/size.hpp>
    #include <boost/mpl/void.hpp>
    #include <boost/preprocessor/repetition.hpp>
    #include <boost/preprocessor/arithmetic/sub.hpp>
    #include <boost/preprocessor/iteration/iterate.hpp>
    #include <boost/preprocessor/cat.hpp>
    #include <boost/preprocessor/if.hpp>
    #include "Definitions.h"
    #include "TupleElements.h"
    
    // define the TupleImpl implementation class
    #define DEFINE_TUPLE_IMPL template<typename TTupleElements, int VSize> class TupleImpl;

    // define the tuple class used to store elements of arbitrary type
    // The elements can be accessed via a bracket operator
    #define DEFINE_TUPLE template<typename TTupleElements> class Tuple : public TupleImpl<TTupleElements, tuple_size<TTupleElements>::value> {};
    
    // define an is tuple trait
    #define DEFINE_IS_TUPLE template<typename T> struct is_tuple : boost::mpl::false_ {};
    #define DEFINE_IS_TUPLE_SPECIALIZATION template<typename TTupleElements> struct is_tuple<Tuple<TTupleElements> > : is_tuple_elements<TTupleElements> {};

    DEFINE_TUPLE_IMPL
    DEFINE_TUPLE
    DEFINE_IS_TUPLE
    DEFINE_IS_TUPLE_SPECIALIZATION
 
    // provide an accelerator specific tuple 
    namespace acc 
    {
        DEFINE_TUPLE_IMPL
        DEFINE_TUPLE
        DEFINE_IS_TUPLE
        DEFINE_IS_TUPLE_SPECIALIZATION
    }

    #undef DEFINE_TUPLE_IMPL
    #undef DEFINE_TUPLE

    // generate specializations
    #define BOOST_PP_ITERATION_PARAMS_1 (3, (0, MAX_TUPLE_SIZE, "Tuple.h"))
    #include BOOST_PP_ITERATE()
 
    #endif // TUPLE_INCLUDED
#else // BOOST_PP_IS_ITERATING
    #define ITERATION_INDEX BOOST_PP_ITERATION()

    // define repeat print functions
    #define TEXT_DEFINE_INDEX_AND_ELEMENT_TYPES(z, n, data) \
        typedef typename element_index<TTupleElements, n >::type BOOST_PP_CAT(I,n); \
        typedef typename element_type<TTupleElements, n >::type BOOST_PP_CAT(T,n);

    #define TEXT_INIT(z, n, data) \
        BOOST_PP_CAT(BOOST_PP_CAT(val, n), _) = BOOST_PP_CAT(val, n);

    #define TEXT_ASSIGN(z, n, data) \
        BOOST_PP_CAT(BOOST_PP_CAT(val, n), _) = BOOST_PP_CAT(BOOST_PP_CAT(other.val, n), _);

    #define TEXT_OPERATOR(z, n, data) \
        TUPLE_QUALIFIER \
        BOOST_PP_CAT(T,n)& operator() (BOOST_PP_CAT(I,n)*) \
        { \
        return BOOST_PP_CAT(BOOST_PP_CAT(val,n), _); \
        }; 

    #define TEXT_CONST_OPERATOR(z, n, data) \
        TUPLE_QUALIFIER \
        typename boost::add_const<BOOST_PP_CAT(T,n)>::type & operator() (BOOST_PP_CAT(I,n)*) const  \
        { \
        return BOOST_PP_CAT(BOOST_PP_CAT(val,n), _); \
        };

    #define TEXT_MEMBER(z, n, data) \
        BOOST_PP_CAT(T,n) BOOST_PP_CAT(BOOST_PP_CAT(val,n),_);

    // define the standard tuple implementation
    #define TUPLE_QUALIFIER __CPU__     
    #include "TupleDefinition.h"
    #undef TUPLE_QUALIFIER

    namespace acc
    {
        #define TUPLE_QUALIFIER __ACC__
        #include "TupleDefinition.h"
        #undef TUPLE_QUALIFIER
    }
  
    #undef TEXT_DEFINE_INDEX_AND_ELEMENT_TYPES
    #undef TEXT_INIT
    #undef TEXT_ASSIGN
    #undef TEXT_OPERATOR
    #undef TEXT_CONST_OPERATOR
    #undef TEXT_MEMBER

    #undef ITERATION_INDEX
#endif // BOOST_PP_IS_ITERATING

