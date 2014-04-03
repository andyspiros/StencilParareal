#pragma once

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/push_front.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/if.hpp>
#include "Enums.h"
#include "Domains.h"
#include "ParameterType.h"

#include <boost/mpl/vector/vector10.hpp> // depends on MAX_BOUNDARY_LEVELS

/**
* @struct compute_begin_boundary_domains
* Meta function computing a list of begin boundary domains
*/
template<
    typename TBaseDomain,
    KLoopDirection VKLoopDirection,
    int VBeginKOffset,
    int VBeginKBoundarySize>
struct compute_begin_boundary_domains
{
    typedef typename boost::mpl::eval_if_c<
        VKLoopDirection == cKIncrement,
        boost::mpl::fold<
            boost::mpl::range_c<int, VBeginKOffset, VBeginKOffset + VBeginKBoundarySize>,
            boost::mpl::vector0<>,
            boost::mpl::push_back<
                boost::mpl::_1, 
                compute_k_minimum_boundary_domain<
                    TBaseDomain,
                    boost::mpl::_2
                >
            >
        >,
        boost::mpl::fold<
            boost::mpl::range_c<int, VBeginKOffset - VBeginKBoundarySize + 1, VBeginKOffset + 1>, // shift range by one
            boost::mpl::vector0<>,
            boost::mpl::push_front<
                boost::mpl::_1, 
                compute_k_maximum_boundary_domain<
                    TBaseDomain,
                    boost::mpl::_2
                >
            >
        >
    >::type type;
};

/**
* @struct compute_end_boundary_domains
* Meta function computing a list of end boundary domains
*/
template<
    typename TBaseDomain,
    KLoopDirection VKLoopDirection,
    int VEndKOffset,
    int VEndKBoundarySize>
struct compute_end_boundary_domains
{
    typedef typename boost::mpl::eval_if_c<
        VKLoopDirection == cKIncrement,
        boost::mpl::fold<
            boost::mpl::range_c<int, VEndKOffset - VEndKBoundarySize, VEndKOffset>,
            boost::mpl::vector0<>,
            boost::mpl::push_back<
                boost::mpl::_1, 
                compute_k_maximum_boundary_domain<
                    TBaseDomain,
                    boost::mpl::_2
                >
            >
        >,
        boost::mpl::fold<
            boost::mpl::range_c<int, VEndKOffset + 1, VEndKOffset + VEndKBoundarySize + 1>, // shift range by one
            boost::mpl::vector0<>,
            boost::mpl::push_front<
                boost::mpl::_1, 
                compute_k_minimum_boundary_domain<
                    TBaseDomain,
                    boost::mpl::_2
                >
            >
        >
    >::type type;
};

/**
* @struct KLoopBoundaryLevelFunctor
* Functor used to apply a boundary level
*/
template<
    typename TIterationMask,
    typename TStencilSweepFunctor,
    KLoopDirection VKLoopDirection>
struct KLoopBoundaryLevelFunctor
{
    template<
        typename TContext, 
        typename TDomain>
    __ACC__
    static void Do(typename parameter_type<TContext>::type context)
    {
        TStencilSweepFunctor::template Do<TDomain, false, false>(context);
        context.template Advance<TIterationMask, 0, 0, VKLoopDirection>();
    }
};





