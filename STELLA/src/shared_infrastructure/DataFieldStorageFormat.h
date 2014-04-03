#pragma once

#include <boost/config.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/mpl/empty.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/if.hpp>
#include "Definitions.h"
#include "DataFieldIJBoundary.h"
#include "DataFieldStorageOrder.h"
#include "DataFieldAlignment.h"

/**
* @struct DataFieldStorageFormat
* Structure defining the format of a data field.
* The IJBoundary parameter specifies the data field boundary in i and j dimension.
* The StorageOrder is a list of length rank which defines the storage order of the different dimensions.
* The TAlignment parameter defines the alignment of the origin in a certain alignment dimension (as a multiple of sizeof(TValue)).
* Note that there is only an alignment if alignment dimension and constant stride dimension match. 
*/
template<
    typename TIJBoundary,
    typename TStorageOrder,
    typename TAlignment>
struct DataFieldStorageFormat
{
    BOOST_STATIC_ASSERT(is_data_field_ij_boundary<TIJBoundary>::value);
    BOOST_STATIC_ASSERT(is_data_field_storage_order<TStorageOrder>::value);
    BOOST_STATIC_ASSERT(is_data_field_alignment<TAlignment>::value);
 
    typedef TIJBoundary IJBoundary;
    typedef TStorageOrder StorageOrder;
    typedef TAlignment Alignment;
};

// define data field type traits

/**
* @struct is_data_field_storage_format
* Meta function returning true the parameter is a data field storage format type
*/
template<typename TStorageFormat>
struct is_data_field_storage_format : boost::mpl::false_ {};

template<
    typename TIJBoundary,
    typename TStorageOrder,
    typename TAlignment>
struct is_data_field_storage_format<DataFieldStorageFormat<TIJBoundary, TStorageOrder, TAlignment> > : boost::mpl::true_ {};

/**
* @struct storage_rank
* Meta function computing the rank of the data field storage
*/
template<typename TStorageFormat>
struct storage_rank : boost::mpl::size<typename TStorageFormat::StorageOrder> {};

/**
* @struct has_storage_in
* Meta function returning true if a data field type has storage in a certain dimension
*/
template<
    typename TStorageFormat,
    Dimension VDimension>
struct has_storage_in : 
    boost::mpl::contains<
        typename TStorageFormat::StorageOrder, 
        boost::mpl::integral_c<Dimension, VDimension> 
    > 
{};

/**
* @struct is_storage_continuous_in
* Meta function returning true if the storage is continuous in a certain dimension
*/
template<
    typename TStorageFormat,
    Dimension VDimension>
struct is_storage_continuous_in 
{
    BOOST_STATIC_CONSTANT(bool, value = 
        (
            boost::mpl::eval_if<
                    boost::mpl::empty<typename TStorageFormat::StorageOrder>,
                    boost::mpl::false_,     
                    boost::is_same<
                        typename boost::mpl::front<typename TStorageFormat::StorageOrder>::type, 
                        boost::mpl::integral_c<Dimension, VDimension>
                    >
            >::type::value
        )
    );
    typedef boost::mpl::integral_c<bool, bool(value)> type;
};

/**
* @struct storage_alignment
* Meta function calculating the number of alignment elements used by a certain data field type
* (note that the minimal alignment is 1)
*/
template<typename TStorageFormat>
struct storage_alignment
{
    BOOST_STATIC_CONSTANT(int, value = 
        (
            boost::mpl::if_<
                is_storage_continuous_in<
                    TStorageFormat,
                    TStorageFormat::Alignment::Dimension::value
                >,
                typename TStorageFormat::Alignment::Elements,
                boost::mpl::integral_c<int, 1>
            >::type::value
        )
    );
    typedef boost::mpl::integral_c<int, int(value)> type;
};

// data field type helper methods

/**
* Method computing the allocation size for a certain dimension
* @param size size input
* @return allocation size
*/
template<
    typename TStorageFormat, 
    Dimension VDimension>
__CPU__
int allocation_size(const int size)
{
    return has_storage_in<TStorageFormat, VDimension>::value ? size : 1;
}

/**
* Method computing the padded size
* @param size size input
* @return padded size
*/
template<
    typename TStorageFormat, 
    Dimension VDimension>
__ACC_CPU__
int padded_size(const int size)
{
    int result = size; 
    const int alignment = storage_alignment<TStorageFormat>::value;
    if(is_storage_continuous_in<TStorageFormat, VDimension>::value)
    {
        result = ((size + alignment - 1) / alignment) * alignment;
    }
    return result;
}

/**
* Method computing an aligned pointer
* @param pointer pointer to align
* @return aligned pointer
*/
template<
    typename TStorageFormat,
    typename TValue>
__CPU__
static TValue* aligned_pointer(TValue* pointer)
{
    const int alignment = storage_alignment<TStorageFormat>::value * sizeof(TValue);
    return reinterpret_cast<TValue*>(((reinterpret_cast<size_t>(pointer) + alignment - 1) / alignment) * alignment); 
}
