#pragma once

#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/if.hpp>
#include "KRange.h"
//#include "KWindow.h"

/**
* @struct StageVariable
* Structure holding the type information needed to define a stage variable.
* @tparam VFieldIndex   unique integer index that identifies this temporary field
* @tparam TValue        storage type of variable (double, int, etc.)
*/
template<
    int VFieldIndex,
    typename TValue>
struct StageVariable
{
    typedef TValue ValueType;
    typedef boost::mpl::integral_c<int, VFieldIndex> FieldIndex;
};

/**
* @struct is_stage_variable
* meta function to determine if a type is a stage variable
*/
template<typename T>
struct is_stage_variable : boost::mpl::false_ {};

// specialization for stage variables
template<
    int VFieldIndex,
    typename TValue>
struct is_stage_variable<StageVariable<VFieldIndex, TValue> > : boost::mpl::true_ {};

/**
* @struct StencilBuffer
* Structure holding the type information needed to define a sweep buffer.
* @tparam VFieldIndex   unique integer index that identifies this temporary field
* @tparam TValue        storage type of variable (double, int, etc.)
* @tparam TKRange      KRange type that specifies the KRange of the buffer
*/
template<
    int VFieldIndex,
    typename TValue,
    typename TKRange>
struct StencilBuffer
{
    BOOST_STATIC_ASSERT(is_k_range<TKRange>::value);
    typedef TValue ValueType;
    typedef TKRange KRange;
    typedef boost::mpl::integral_c<int, VFieldIndex> FieldIndex; 
};

/**
* @struct is_temporary_field
* Meta function returning true if T is a temporary field
*/
template<typename T>
struct is_temporary_field : boost::mpl::false_ {};

// specialization for stencil buffers
template<
    int VFieldIndex,
    typename TValue,
    typename TKRange>
struct is_temporary_field<StencilBuffer<VFieldIndex, TValue, TKRange> > : boost::mpl::true_ {};

// specialization for stage variables
template<
    int VFieldIndex,
    typename TValue>
struct is_temporary_field<StageVariable<VFieldIndex, TValue> > : boost::mpl::true_ {};

/**
* @struct temporary_field_index
* Meta function returning the temporary field index
*/
template<typename T>
struct temporary_field_index;

// specialization for stencil buffers
template<
    int VFieldIndex,
    typename TValue,
    typename TKRange>
struct temporary_field_index<StencilBuffer<VFieldIndex, TValue, TKRange> >
    : boost::mpl::integral_c<int, VFieldIndex>
{};

// specialization for stage variables
template<
    int VFieldIndex,
    typename TValue>
struct temporary_field_index<StageVariable<VFieldIndex, TValue> >
    : boost::mpl::integral_c<int, VFieldIndex>
{};


/**
* @struct SweepBuffer
* Structure holding the type information needed to define a sweep buffer.
* @tparam VFieldIndex   unique integer index that identifies this temporary field
* @tparam TValue        storage type of variable (double, int, etc.)
* @tparam TKRange       KRange of the buffer
* @tparam TKWindow      KWindow of the buffer
*/
/*
template <
    int VFieldIndex,
    typename TValue,
    typename TKRange,
    typename TKWindow>
struct SweepBuffer
{
    BOOST_STATIC_ASSERT((is_sweep_buffer_storage_format<TKRange, TKWindow>::value));
    typedef TValue ValueType;
    typedef TKWindow WindowType;
    typedef boost::mpl::integral_c<int, VFieldIndex> FieldIndex; 
};
*/
/*
// specialization for sweep buffers
template<
    int VFieldIndex,
    typename TValue,
    typename TKRange,
    typename TKWindow>
struct temporary_field_index<SweepBuffer<VFieldIndex, TValue, TKRange, TKWindow> >
    : boost::mpl::integral_c<int, VFieldIndex>
{};
*/
/*
// specialization for sweep buffers
template<
    int VFieldIndex,
    typename TValue,
    typename TKRange,
    typename TKWindow>
struct is_temporary_field<SweepBuffer<VFieldIndex, TValue, TKRange, TKWindow> >
    : boost::mpl::true_
{};
*/

