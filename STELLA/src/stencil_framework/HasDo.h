#pragma once

#include <boost/config.hpp> 
#include <boost/mpl/void.hpp>
#include "YesNo.h"

/**
* @struct has_do_member
* Meta function testing if a functor has a Do member (function, variable or typedef)
*/
template<typename TFunctor>
struct has_do_member
{
    // define a MixIn class providing a Do member
    struct MixIn
    {
        void Do() {}
    };
    // derive from MixIn and TStencilStageDefinition
    struct Derived : public TFunctor, public MixIn {};

    // define an SFINAE structure
    template <typename TDoFunc, TDoFunc VFunc>
    struct SFINAE{};

    // if TFunctor has no Do method the first test method matches
    // otherwise the ellipsis test method matches due to failing template parameter substitution. 
    // The template parameter substitution fails due to ambiguous Do methods
    template<typename TDerived>
    static no test(SFINAE<void (MixIn::*) (), &TDerived::Do>*);
    template<typename TDerived>
    static yes test(...); 

    // use the sizeof trick in order to check which overload matches
    BOOST_STATIC_CONSTANT(bool, value = (sizeof(test<Derived>(0)) == sizeof(yes)) );
    typedef boost::mpl::bool_<bool(value)> type;
};

// setup the HasDoDetails namespace which provides a check function searching for a specific Do method signature
// (SFINAE is used with a comma operator matching the return value of the Do function in case
// void is return we fall back to the default comma operator returning the element after the comma)
namespace HasDoDetails
{
    // define specific comma operators used for SFINAE
    template<typename T>
    T& operator, (T&, boost::mpl::void_);
    template<typename T>
    const T& operator, (const T&, boost::mpl::void_);

    // if the result is not void make sure the specific comma operator matched and returned TResult
    template<typename TResult> 
    struct check_result
    {
        static yes test(TResult); // matching Do method found which returns TResult
        static no test(boost::mpl::void_); // matching Do method found which returns void
        static no test(no); // no matching Do method found except no Do(...)
        static no test(...); 
    };

    // if the result has to be void make sure the default comma operator matched and returned boost::mpl::void__
    template<> 
    struct check_result<void>
    {        
        static yes test(boost::mpl::void_); // matching Do method found which returns void
        static no test(no); // no matching Do method found except no Do(...)
        static no test(...); 
    };

    // if there is no Do member 
    template<
        bool VCondition, 
        typename TFunctor, 
        typename TResult,
        typename TContext, 
        typename TDomain> 
    struct check_if
    {
        BOOST_STATIC_CONSTANT(bool, value = false);
        typedef boost::mpl::bool_<bool(value)> type;
    };

    // if there is a Do member but no domain parameter check for TResult Do(TContext)
    template<
        typename TFunctor,
        typename TResult,
        typename TContext>
    struct check_if<true, TFunctor, TResult, TContext, boost::mpl::void_>
    {
        // extend TFunctor with a fall back Do method and import all TFunctor Do methods
        // (Note that it is crucial to import the TFunctor Do methods as the overload
        // resolution only considers Do methods imported in the scope of Derived)
        struct Derived : public TFunctor
        { 
            using TFunctor::Do; 
            static no Do(...); 
        }; 

        // SFINAE using the comma operator!
        // if the template parameter substitution for the comma operator fails (as Do returns void)
        // the default comma operator is used which returns boost::mpl::void_
        BOOST_STATIC_CONSTANT(bool, value = (
            sizeof(
                check_result<TResult>::test(
                    (
                        Derived::Do(*(TContext*)0)
                        ,
                        boost::mpl::void_()
                    )
                )
            ) == sizeof(yes))
        );
        typedef boost::mpl::bool_<bool(value)> type;
    };

    // if there is a Do member with a domain parameter check for TResult Do(TContext, TDomain)
    template<
        typename TFunctor, 
        typename TResult,
        typename TContext, 
        typename TDomain> 
    struct check_if<true, TFunctor, TResult, TContext, TDomain>
    {   
        // extend TFunctor with a fall back Do method and import all TFunctor Do methods
        // (Note that it is crucial to import the TFunctor Do methods as the overload
        // resolution only considers Do methods imported in the scope of Derived)
        struct Derived : public TFunctor
        { 
            using TFunctor::Do; 
            static no Do(...); 
        }; 

        // SFINAE using the comma operator!
        // if the template parameter substitution for the comma operator fails (as Do returns void)
        // the default comma operator is used which returns boost::mpl::void_
        BOOST_STATIC_CONSTANT(bool, value = (
            sizeof(            
                check_result<TResult>::test(
                    (
                        Derived::Do(*(TContext*)0, *(TDomain*)0)
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
* @struct has_do
* Meta function testing if a functor has a specific Do methods signature
*/
template<
    typename TFunctor, 
    typename TResult,
    typename TContext, 
    typename TDomain> 
struct has_do
{ 
    BOOST_STATIC_CONSTANT(bool, value = (
        HasDoDetails::check_if<
            has_do_member<TFunctor>::value, 
            TFunctor,
            TResult,
            TContext,
            TDomain            
        >::type::value)
    );
    typedef boost::mpl::bool_<bool(value)> type;
}; 
