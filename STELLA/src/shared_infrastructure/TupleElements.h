#pragma once

#include <boost/static_assert.hpp>

#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/at.hpp> 
#include <boost/mpl/apply.hpp>

#include <boost/mpl/vector/vector40_c.hpp> // depends on MAX_TUPLE_SIZE
#include <boost/mpl/vector/vector40.hpp> // depends on MAX_TUPLE_SIZE

/**
* @struct TupleElements
* Structure combining a vector of all element indexes with a map defining the associated element types
*/
template<
    typename TElementIndexes,
    typename TElementTypes>
struct TupleElements
{   
    BOOST_STATIC_ASSERT(boost::mpl::is_sequence<TElementIndexes>::value);
    BOOST_STATIC_ASSERT(boost::mpl::is_sequence<TElementTypes>::value);
    BOOST_STATIC_ASSERT(boost::mpl::size<TElementIndexes>::value == boost::mpl::size<TElementTypes>::value);

    // define an index and an element type list
    typedef TElementIndexes ElementIndexes;
    typedef TElementTypes ElementTypes;
};

/**
* @struct is_tuple_elements
* Meta function returning true if the parameter is a TupleElements
*/
template<typename T>
struct is_tuple_elements : boost::mpl::false_ {};

template<
    typename TElementIndexes,
    typename TElementTypes>
struct is_tuple_elements<TupleElements<TElementIndexes, TElementTypes> > : boost::mpl::true_ {};

/**
* @struct element_type
* Meta function returning the element type at a certain index 
*/
template<
    typename TTupleElements,
    int VElementIndex>
struct element_type;

template<
    typename TElementIndexes,
    typename TElementTypes,
    int VElementIndex>
struct element_type<TupleElements<TElementIndexes, TElementTypes>, VElementIndex>
{
    typedef typename boost::mpl::at_c<TElementTypes, VElementIndex>::type type; 
};

/**
* @struct element_index
* Meta function returning the element index at a certain index 
*/
template<
    typename TTupleElements,
    int VElementIndex>
struct element_index;

template<
    typename TElementIndexes,
    typename TElementTypes,
    int VElementIndex>
struct element_index<TupleElements<TElementIndexes, TElementTypes>, VElementIndex>
{
    typedef typename boost::mpl::at_c<TElementIndexes, VElementIndex>::type type; 
};

/**
* @struct tuple_size
* Meta function returning tuple size
*/
template<typename TTupleElements>
struct tuple_size;

template<
    typename TElementIndexes,
    typename TElementTypes>
struct tuple_size<TupleElements<TElementIndexes, TElementTypes> > : boost::mpl::size<TElementIndexes> {};

/**
* @struct merge_tuple_elements
* Meta function merging two tuple elements,
* transform operations allow modifying the element types at the same time
*/
template<
    typename TTupleElements1,
    typename TTupleElements2,
    typename TTransformOp1,
    typename TTransformOp2>
struct merge_tuple_elements
{
    // compute the new indexes
    typedef typename boost::mpl::fold<
        typename TTupleElements2::ElementIndexes,
        typename TTupleElements1::ElementIndexes,
        boost::mpl::push_back<
            boost::mpl::_1, 
            boost::mpl::_2
        >
    >::type ElementIndexes;

    // compute the new element types
    typedef typename boost::mpl::fold<
        typename TTupleElements2::ElementTypes,
        typename boost::mpl::transform<
            typename TTupleElements1::ElementTypes,
            TTransformOp1
        >::type,
        boost::mpl::push_back<
            boost::mpl::_1, 
            typename boost::mpl::apply<TTransformOp2, boost::mpl::_2>::type // apply the parameter 2 placeholder to the transform operation
        >
    >::type ElementTypes;

    // assemble the resulting tuple elements structure
    typedef TupleElements<ElementIndexes, ElementTypes> type;
};



