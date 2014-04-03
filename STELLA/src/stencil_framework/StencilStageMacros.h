#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/or.hpp>
#include "Definitions.h"
#include "YesNo.h"
#include "DataParameter.h"
#include "Domains.h"

// macro used to define a stencil stage
#define STENCIL_STAGE(environment) \
    BOOST_STATIC_ASSERT(is_stencil_stage_environment<environment>::value); \
    typedef environment StencilStageEnvironment; \
    typedef typename environment::ValueType T; \
    typedef typename environment::Context Context; 

// macro used for the stage parameter definition
// it provides a type of name param which is used to access the parameter inside the context tuple
#define STAGE_PARAMETER(domain, param) \
    BOOST_MPL_ASSERT_MSG( \
        is_2d_domain<domain>::value || is_3d_domain<domain>::value, \
        CHECK_IF_THE_DOMAIN_##domain##_OF_##param##_IS_VALID, \
        (domain) \
    ); \
    typedef typename boost::mpl::if_< \
        is_sub_domain<typename StencilStageEnvironment::Domain, domain>, \
        DataParameter<(int)::param, Offset<0,0,0> >, \
        DataParameter<-1, Offset<0,0,0> > \
    >::type param; \
    typedef typename boost::mpl::if_< \
        boost::mpl::or_< \
            is_sub_domain<typename StencilStageEnvironment::Domain, domain>, \
            is_sub_domain<domain, typename StencilStageEnvironment::Domain> \
        >, \
        yes, \
        no \
    >::type has##param; \
    static has##param Parameter(boost::mpl::integral_c<int, (int)::param>*); 

// -> parameters are only accessible if the stage domain is sub domain of the parameter domain
// -> the Parameter() method returns yes if the stage and the parameter domain overlap
// (note that this is necessary as the Parameter() method is use to check if a parameter is used somewhere inside the stage domain)



  
