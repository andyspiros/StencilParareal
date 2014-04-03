#pragma once

#include <boost/mpl/assert.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/transform.hpp>
#include "StencilFunctionEnvironment.h"
#include "HasDo.h"
#include "DataParameter.h"
#include "FunctionParameter.h"
#include "Offset.h"

/**
* @struct shift_parameter
* Meta function shifting a parameter by an offset
*/
template<typename TParameter, typename TParameterOffset>
struct shift_parameter
{
    // by default do nothing e.g. finite difference offset parameters
    typedef TParameter type;
};

template<
    int VParameterIndex,
    typename TOffset,
    typename TParameterOffset>
struct shift_parameter<DataParameter<VParameterIndex, TOffset>, TParameterOffset>
{
    typedef DataParameter<
        VParameterIndex, 
        typename add_offsets<
            TOffset, 
            TParameterOffset
        >::type
    > type;
};

template<
    template<typename> class TStencilFunction,
    typename TDomain,
    typename TOffset,
    typename TParameterList,
    typename TParameterOffset>
struct shift_parameter<FunctionParameter<TStencilFunction, TDomain, TOffset, TParameterList>, TParameterOffset>
{
    typedef FunctionParameter<
        TStencilFunction, 
        TDomain,
        typename add_offsets<
            TOffset, 
            TParameterOffset
        >::type, 
        TParameterList
    > type;
};

/**
* @struct FunctionParameterEvaluator
* Structure which provides the type info used by the context to evaluate a function parameter
*/
template<
    typename TContext,
    template<typename> class TStencilFunction,
    typename TDomain,
    typename TOffset,
    typename TParameterList>
struct FunctionParameterEvaluator
{
    // define the function
    typedef TStencilFunction<
        StencilFunctionEnvironment<
            TContext,
            typename boost::mpl::transform<
                TParameterList,
                shift_parameter<boost::mpl::_, TOffset>
            >::type
        >
    > StencilFunction;

    // check if the method returns a value type or void
    typedef typename has_do<StencilFunction, typename TContext::ValueType, TContext, TDomain>::type ReturnValueType;
    typedef typename has_do<StencilFunction, void, TContext, TDomain>::type ReturnVoid;

    // check that a do method with a correct signature is provided
    // static T Do(Context) -> if TDomain == boost::mpl::void_
    // static T Do(Context, Domain)
    // static void Do(Context) -> if TDomain == boost::mpl::void_
    // static void Do(Context, Domain)
    BOOST_MPL_ASSERT_MSG(
        ReturnValueType::value || ReturnVoid::value,
        STENCIL_FUNCTION_DOES_NOT_PROVIDE_CORRECT_DO_METHOD_SIGNATURE,
        (StencilFunction)
    );

    // set the type to value type or void
    typedef typename boost::mpl::if_<
        ReturnValueType,
        typename TContext::ValueType,
        void
    >::type ReturnType;

    // call the stencil function do method
    // (note in case the TDomain is boost::mpl::void_ use the lean Do signature taking only the context parameter
    __ACC__
    static ReturnType Evaluate(TContext& ctx)
    {
        return evaluate(ctx, static_cast<typename boost::mpl::is_void_<TDomain>::type*>(0));
    }

private:
    // call do implementations with and without domain parameter
    __ACC__
    static ReturnType evaluate(TContext& ctx, boost::mpl::true_*) { return StencilFunction::Do(ctx); }
    __ACC__
    static ReturnType evaluate(TContext& ctx, boost::mpl::false_*) { return StencilFunction::Do(ctx, TDomain()); }
};

  
