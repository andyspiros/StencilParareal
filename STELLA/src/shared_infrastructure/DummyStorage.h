#pragma once

#include <boost/mpl/bool.hpp>
#include "ParameterTraits.h"

/**
* @class DummyStorage
* Dummy storage class used to implement temporary fields which are not accessed due to caching
*/
template<typename TValue>
class DummyStorage
{
public:
    typedef TValue ValueType;
    typedef ValueType& ReferenceType;
    typedef ValueType* PointerType;
   
    // define the return type
    typedef ReferenceType ReturnType;

    __ACC_CPU__
    DummyStorage() {}
    __ACC_CPU__
    ~DummyStorage() {}
    
    // copy constructor and assignment
    __ACC_CPU__
    DummyStorage(const DummyStorage& other) { *this = other; }
    __ACC_CPU__
    DummyStorage& operator= (const DummyStorage& other)
    {
        // by convention, always return *this
        return *this;
    }
};

/**
* @struct is_dummy_storage
* Meta function returning true the parameter is a dummy storage
*/
template<typename T>
struct is_dummy_storage : boost::mpl::false_ {};

template<typename TValue>
struct is_dummy_storage<DummyStorage<TValue> > : boost::mpl::true_ {};

// return type specialization
template<typename TValue>
struct return_type<DummyStorage<TValue> >
{
    typedef typename DummyStorage<TValue>::ReturnType type;
};

// value type specialization
template<typename TValue>
struct value_type<DummyStorage<TValue> >
{
    typedef typename DummyStorage<TValue>::ValueType type;
};

// is iterable specialization
template<typename TValue>
struct is_iterable<DummyStorage<TValue> > : boost::mpl::false_ {};

