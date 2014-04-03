#pragma once

#include <cassert>
#include <boost/type_traits/add_const.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/type_traits/is_enum.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/and.hpp>
#include "Enums.h"

/**
* @class ParameterWrapper
* Class storing a pointer to a stencil parameter. The parameter is usually a data field or a scalar parameter.
* The access policy parameter defines if the parameter is a read only or a read write parameter.
*/
template<
    typename TValue,
    ParameterIntend VParameterIntend,
    int VParameterIndex>
class ParameterWrapper 
{
public:
    // define the parameter type info
    typedef boost::mpl::integral_c<int, VParameterIndex> ParameterIndex;
    
    // define a value type depending on access policy
    typedef typename boost::mpl::if_c<
        VParameterIntend == cInOut, 
        TValue, 
        typename boost::add_const<TValue>::type
    >::type ValueType;
    typedef ValueType& ReferenceType;
    typedef ValueType* PointerType;
    
    __CPU__
    ParameterWrapper() { pValue_ = NULL; }
    __CPU__
    ~ParameterWrapper() {}

    __CPU__
    explicit ParameterWrapper(ReferenceType value) { pValue_ = &value; }

    __CPU__
    ParameterWrapper(const ParameterWrapper& other) { *this = other; }
    __CPU__
    ParameterWrapper& operator= (const ParameterWrapper& other) 
    { 
        pValue_ = other.pValue_; 
        return *this; 
    }

    /**
    * @return dereferenced value
    */
    __CPU__
    ReferenceType Unwrap()
    {
        assert(pValue_);
        return *pValue_;         
    }

private:
    PointerType pValue_;
};

/**
* @struct parameter_wrapper_type
* Meta function returning the type of the wrapped element if T is a ParameterWrapper
*/
template<typename T>
struct parameter_wrapper_type;

template<
    typename TValue,
    ParameterIntend VParameterIntend,
    int VParameterIndex>
struct parameter_wrapper_type<ParameterWrapper<TValue, VParameterIntend, VParameterIndex> >
{
    typedef TValue type;
};

/**
* @struct parameter_wrapper_index
* Meta function returning the index of the wrapped parameter
*/
template<typename T>
struct parameter_wrapper_index;

template<
    typename TValue,
    ParameterIntend VParameterIntend,
    int VParameterIndex>
struct parameter_wrapper_index<ParameterWrapper<TValue, VParameterIntend, VParameterIndex> > :
    boost::mpl::integral_c<int, VParameterIndex>
{};

/**
* @struct parameter_wrapper_intend
* Meta function returning the intend of the wrapped parameter
*/
template<typename T>
struct parameter_wrapper_intend;

template<
    typename TValue,
    ParameterIntend VParameterIntend,
    int VParameterIndex>
struct parameter_wrapper_intend<ParameterWrapper<TValue, VParameterIntend, VParameterIndex> > :
    boost::mpl::integral_c<ParameterIntend, VParameterIntend>
{};

/**
* Method creating a ParameterWrapper
* @param value wrapped
*/
template<
    int VParameterIndex, 
    ParameterIntend VParameterIntend, 
    typename TValue>
__CPU__
ParameterWrapper<typename boost::remove_const<TValue>::type, VParameterIntend, VParameterIndex> Param(TValue& value)
{
    // check the parameter type
    BOOST_MPL_ASSERT_MSG(
        (
            boost::mpl::or_<
                boost::mpl::and_<
                    boost::mpl::or_<
                        is_data_field<TValue>, 
                        is_joker_data_field<TValue> 
                    >,
                    boost::mpl::bool_<VParameterIntend != cScalar>
                >,
                boost::mpl::and_<
                    boost::mpl::or_<
                        boost::is_arithmetic<TValue>, 
                        boost::is_enum<TValue> 
                    >,
                    boost::mpl::bool_<VParameterIntend == cScalar>
                >
            >::value
        ),
        PARAM_PARAMETER_WRAPPER_TYPE_OR_INTEND_NOT_VALID,
        (TValue) 
    );

    return ParameterWrapper<typename boost::remove_const<TValue>::type, VParameterIntend, VParameterIndex>(value);
}

