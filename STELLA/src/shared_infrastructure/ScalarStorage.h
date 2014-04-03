#pragma once

#include <boost/type_traits/is_pod.hpp>
#include <boost/type_traits/add_const.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/bool.hpp>
#include "Enums.h"
#include "ParameterTraits.h"

/**
* @class ScalarStorage
* Class storing a scalar value 
*/
template<
    typename TValue, 
    AccessPolicy VAccessPolicy>
class ScalarStorage
{
public:
    typedef typename boost::mpl::if_c<
        VAccessPolicy == cReadOnly, 
        typename boost::add_const<TValue>::type,
        typename boost::remove_const<TValue>::type
    >::type ValueType;
    typedef ValueType& ReferenceType;
    typedef ValueType* PointerType;
   
    // define the return type
    typedef typename boost::mpl::if_c<
        boost::is_pod<ValueType>::value && VAccessPolicy == cReadOnly, 
        ValueType,
        ReferenceType
    >::type ReturnType;

    __ACC_CPU__
    ScalarStorage() {}
    __ACC_CPU__
    ~ScalarStorage() {}
    
    // copy constructor and assignment
    __ACC_CPU__
    ScalarStorage(const ScalarStorage& other) { *this = other; }
    __ACC_CPU__
    ScalarStorage& operator= (const ScalarStorage& other)
    {
        value_ = other.value_;
        // by convention, always return *this
        return *this;
    }

    /**
    * Set the value
    */
    __ACC_CPU__
    void set_value(ReferenceType value) { value_ = value; }

    /**
    * @return the value
    */
    __ACC_CPU__
    ReturnType value() { return value_; }

    /**
    * @return the value
    */
    __ACC_CPU__
    ReturnType value() const { return value_; } 

private:
    typename boost::remove_const<ValueType>::type value_;
};

// return type specialization
template<
    typename TValue, 
    AccessPolicy VAccessPolicy>
struct return_type<ScalarStorage<TValue, VAccessPolicy> >
{
    typedef typename ScalarStorage<TValue, VAccessPolicy>::ReturnType type;
};

// value type specialization
template<
    typename TValue, 
    AccessPolicy VAccessPolicy>
struct value_type<ScalarStorage<TValue, VAccessPolicy> >
{
    typedef typename ScalarStorage<TValue, VAccessPolicy>::ValueType type;
};

// is iterable specialization
template<
    typename TValue, 
    AccessPolicy VAccessPolicy>
struct is_iterable<ScalarStorage<TValue, VAccessPolicy> > : boost::mpl::false_ {};

