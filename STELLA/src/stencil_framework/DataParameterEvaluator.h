#pragma once

#include <boost/mpl/assert.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/has_key.hpp>
#include <boost/mpl/at.hpp>

/**
* @struct DataParameterEvaluator
* Structure which provides the type info used by the context to evaluate a data parameter
*/
template<
    typename TIndexToElementMap,
    int VParameterIndex>
struct DataParameterEvaluator
{
    // check that the parameter index is not -1
    // this happens if the parameter is invalid in the domain
    // or if not enough parameters where handed over to a stencil function
    BOOST_MPL_ASSERT_MSG(
        VParameterIndex != -1,
        STAGE_PARAMETER_DOMAIN_OR_FUNCTION_STAGE_PARAMETERLIST_INCORRECT,
        (TIndexToElementMap)
    );
    
    // check that the parameter index is part of the element index
    BOOST_MPL_ASSERT_MSG(
        (
            boost::mpl::has_key<
                TIndexToElementMap, 
                boost::mpl::integral_c<int, VParameterIndex> 
            >::value
        ),
        PARAMETER_INDEX_FOUND_MORE_THAN_ONCE_OR_NOT_AT_ALL,
        (TIndexToElementMap)
    );

    // compute the return type
    typedef typename return_type<
        typename boost::mpl::at<
            TIndexToElementMap, 
            boost::mpl::integral_c<int, VParameterIndex> 
        >::type
    >::type ReturnType;
};

