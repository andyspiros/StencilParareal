#pragma once

#include <boost/config.hpp>
#include <boost/static_assert.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/bool.hpp>
#include "Enums.h"
#include "Definitions.h"

/**
* @struct KPosition
* Structure defining a k position
*/
template<
    KReferencePositions VKReference,
    int VKOffset>
struct KPosition
{
    // define reference position and offset
    typedef boost::mpl::integral_c<KReferencePositions, VKReference> KReference;
    typedef boost::mpl::integral_c<int, VKOffset> KOffset;
};

/**
* @struct is_k_position
* Meta function returning true if the parameter is a k position
*/
template<typename TKPosition>
struct is_k_position : boost::mpl::false_ {};

template<
    KReferencePositions VKReference,
    int VKOffset>
struct is_k_position<KPosition<VKReference, VKOffset> > : boost::mpl::true_ {};

/**
* @struct shift_k_position_by
* Meta function shifting k position
*/
template<
    typename TKPosition,
    int VKShift>
struct shift_k_position_by;

template<
    KReferencePositions VKReference,
    int VKOffset,
    int VKShift>
struct shift_k_position_by<KPosition<VKReference, VKOffset>, VKShift>
{
    typedef KPosition<VKReference, VKOffset + VKShift> type;
};

/**
* @struct k_position_compile_time_origin_offset
* Meta function returning the compile time k origin offset of a k position
* (assumption the run time values might be different but the order remains valid, 
* e.g. the k size might be different from cDefault size but not smaller than the flat limit!)
*/
template<typename TKPosition>
struct k_position_compile_time_origin_offset;

template<
    KReferencePositions VKReference,
    int VKOffset>
struct k_position_compile_time_origin_offset<KPosition<VKReference, VKOffset> >
{
    BOOST_STATIC_CONSTANT(int, value = 
        (
            VKOffset +
            (VKReference == cKMaximumTerrain ? cDefaultKSize - 1 : 0) +
            (VKReference == cKMinimumTerrain ? cFlatLimit : 0) +
            (VKReference == cKMaximumFlat ? cFlatLimit - 1 : 0) 
        )
    );
    typedef boost::mpl::integral_c<int, int(value)> type;
};

/**
* @struct k_positions_direction
* Meta function returning the direction of two positions
* (note that zero is returned if the positions are equivalent)
*/
template<
    typename TKPosition1,
    typename TKPosition2>
struct k_positions_direction
{
    BOOST_STATIC_CONSTANT(int, value = 
        (
            k_position_compile_time_origin_offset<TKPosition2>::value > 
            k_position_compile_time_origin_offset<TKPosition1>::value ? 
            1 : -1
        )
    );
    typedef boost::mpl::integral_c<int, int(value)> type;
};

// specialization returning 0 in case the k positions are equivalent
template<
    KReferencePositions VKReference,
    int VKOffset>
struct k_positions_direction<KPosition<VKReference, VKOffset>, KPosition<VKReference, VKOffset> > : boost::mpl::integral_c<int, 0> {};

/**
* @struct k_positions_compile_time_length_component
* Meta function returning the compile time length of the range defined by two k positions
*/
template<
    typename TKPosition1,
    typename TKPosition2>
struct k_positions_compile_time_length_component
{
    BOOST_STATIC_CONSTANT(int, value = 
        (
            k_positions_direction<TKPosition1, TKPosition2>::value * 
            (
                k_position_compile_time_origin_offset<TKPosition2>::value -
                k_position_compile_time_origin_offset<TKPosition1>::value
            )
        )
    );
    typedef boost::mpl::integral_c<int, int(value)> type;
};

/**
* @struct have_k_positions_run_time_length_component
* Meta function returning true if the range defined by two k positions has a runtime component
*/
template<
    typename TKPosition1,
    typename TKPosition2>
struct have_k_positions_run_time_length_component;

template<
    KReferencePositions VKReference1,
    int VKOffset1,
    KReferencePositions VKReference2,
    int VKOffset2>
struct have_k_positions_run_time_length_component<KPosition<VKReference1, VKOffset1>, KPosition<VKReference2, VKOffset2> >
{
    BOOST_STATIC_CONSTANT(bool, value = 
        (
            VKReference1 != VKReference2 && 
            (VKReference1 == cKMaximumTerrain || VKReference2 == cKMaximumTerrain)
        )
    );
    typedef boost::mpl::integral_c<bool, bool(value)> type;
};


