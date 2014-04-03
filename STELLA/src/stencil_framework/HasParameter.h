#pragma once

#include <boost/config.hpp> 
#include <boost/mpl/bool.hpp>
#include <boost/mpl/void.hpp>
#include "YesNo.h"

/**
* @struct has_parameter_member
* Meta function testing if a functor has a Parameter member (function, variable or typedef)
*/
template<typename TFunctor>
struct has_parameter_member
{
    // define a MixIn class providing a Parameter member
    struct MixIn
    {
        void Parameter() {}
    };
    // derive from MixIn and TFunctor
    struct Derived : public TFunctor, public MixIn {};
    
    // define an SFINAE structure
    template <typename TParameterFunc, TParameterFunc VFunc>
    struct SFINAE{};

    // if TStencilStageDefinition has no Parameter method the first test method matches
    // otherwise the ellipsis test method matches due to failing template parameter substitution. 
    // The template parameter substitution fails due to ambiguous Parameter methods
    template<typename TDerived>
    static no test(SFINAE<void (MixIn::*) (), &TDerived::Parameter>*);
    template<typename TDerived>
    static yes test(...);
    
    // use the sizeof trick in order to check which overload matches
    BOOST_STATIC_CONSTANT(bool, value = (sizeof(test<Derived>(0)) == sizeof(yes)) );
    typedef boost::mpl::bool_<bool(value)> type;
};

// setup the HasParameterDetails namespace which provides a check function searching for a specific Parameter method signature
// (SFINAE is used with a comma operator matching the return value of the Parameter function in case
// void is return we fall back to the default comma operator returning the element after the comma)
namespace HasParameterDetails
{
    template<typename T>
    T& operator, (T&, boost::mpl::void_);
    template<typename T>
    const T& operator, (const T&, boost::mpl::void_);

    // check if Parameter method returning yes was found
    struct check_result
    {
        static yes test(yes); // matching Parameter method found which returns yes
        static no test(no); // no matching Parameter method found except no Parameter(...)
        static no test(boost::mpl::void_); // should not happen as Parameter methods never return void 
        static no test(...);           
    };

    // if there is no Parameter member 
    template<
        bool VCondition, 
        typename TDerived, 
        typename TParameterIndex> 
    struct check_if
    {
        BOOST_STATIC_CONSTANT(bool, value = false);
        typedef boost::mpl::bool_<bool(value)> type;
    };

    // if there is a Parameter member check for yes Parameter(TParameterIndex*)
    template<
        typename TFunctor, 
        typename TParameterIndex> 
    struct check_if<true, TFunctor, TParameterIndex>
    {   
        // extend TFunctor with a fall back Parameter method and import all TFunctor Parameter methods
        // (Note that it is crucial to import the TFunctor Parameter methods as the overload
        // resolution only considers Parameter methods imported in the scope of Derived)
        struct Derived : public TFunctor
        { 
            using TFunctor::Parameter; 
            static no Parameter(...); 
        }; 

        // SFINAE using the comma operator!
        // if the template parameter substitution for the comma operator fails (as Parameter returns void)
        // the default comma operator is used which returns boost::mpl::void_
        BOOST_STATIC_CONSTANT(bool, value = (
            sizeof(            
                check_result::test(
                    (
                        Derived::Parameter((TParameterIndex*)0)
                        , 
                        boost::mpl::void_()
                    )
                )
            ) == sizeof(yes))
        );
        typedef boost::mpl::bool_<bool(value)> type;
    };
}

/**
* @struct has_parameter
* Meta function testing if a functor has a specific Parameter method signature
*/
template<
    typename TFunctor,
    typename TParameterIndex> 
struct has_parameter
{ 
    BOOST_STATIC_CONSTANT(bool, value = (
        HasParameterDetails::check_if<
            has_parameter_member<TFunctor>::value,
            TFunctor,
            TParameterIndex  
        >::type::value)
    );
    typedef boost::mpl::bool_<bool(value)> type;
}; 
