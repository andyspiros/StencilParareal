#pragma once

#include <boost/type_traits/is_pointer.hpp>
#include <boost/type_traits/is_reference.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/type_traits/add_reference.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/or.hpp>

/**
* @struct parameter_type
* Meta function returning the parameter type of a type T.
* Pointers, references and arithmetic values are passed by value (everything else is passed by reference).
*/
template<typename T>
struct parameter_type : 
    boost::mpl::eval_if<
        boost::mpl::or_<boost::is_pointer<T>, boost::is_reference<T>, boost::is_arithmetic<T> >, 
        boost::mpl::identity<T>, 
        boost::add_reference<T> 
    >
{};



