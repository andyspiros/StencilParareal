#pragma once

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/empty.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/if.hpp>
#include "Enums.h"
#include "Definitions.h"
#include "KPosition.h"
#include "KRange.h"

/**
* @struct KLoopAdvance
* Structure used by the k loop to advance the parameter positions form a certain position to another position
*/
template<
    typename TFromKPosition,
    typename TToKPosition,
    typename TIterationMask,
    int VIOffset, 
    int VJOffset>
struct KLoopAdvance
{
    BOOST_STATIC_ASSERT(is_k_position<TFromKPosition>::value);
    BOOST_STATIC_ASSERT(is_k_position<TToKPosition>::value);

    template<typename TContext>
    __ACC__
    static void Do(TContext& context)
    {
        // select the implementation
        // if there is no work (empty iteration mask or advance by zero) select the empty implementation
        // otherwise select the runtime respectively the compile-time implementation
        typedef typename boost::mpl::eval_if_c<
            (
                boost::mpl::empty<TIterationMask>::value ||
                (boost::is_same<TFromKPosition, TToKPosition>::value && VIOffset == 0 && VJOffset == 0)
            ),
            boost::mpl::integral_c<ImplementationPolicy, cEmpty>, 
            boost::mpl::if_<
                have_k_positions_run_time_length_component<TFromKPosition, TToKPosition>,
                boost::mpl::integral_c<ImplementationPolicy, cRuntime>, 
                boost::mpl::integral_c<ImplementationPolicy, cCompileTime>
            >
        >::type Implementation;

        // call the selected implementation
        doImpl(context, static_cast<Implementation*>(0));
    }

private:
    // empty implementation called if no work is needed -> e.g. advance by 0
    template<typename TContext>
    __ACC__
    static void doImpl(TContext& context, boost::mpl::integral_c<ImplementationPolicy, cEmpty>*) {}

    // implementation working with a compile-time offset
    template<typename TContext>
    __ACC__
    static void doImpl(TContext& context, boost::mpl::integral_c<ImplementationPolicy, cCompileTime>*) 
    {
        // compute the offset
        const int kOffset =
            k_positions_direction<TFromKPosition, TToKPosition>::value * 
            k_positions_compile_time_length_component<TFromKPosition, TToKPosition>::value; 

        context.template Advance<
            TIterationMask,
            VIOffset,
            VJOffset,
            kOffset
        >();   
    }

    // implementation computing a runtime offset
    template<typename TContext>
    __ACC__
    static void doImpl(TContext& context, boost::mpl::integral_c<ImplementationPolicy, cRuntime>*) 
    {
        // compute the run time offset 
        // --> add the deviation between default k size and runtime k size
        const int kOffset = 
            k_positions_direction<TFromKPosition, TToKPosition>::value * 
            (
                k_positions_compile_time_length_component<TFromKPosition, TToKPosition>::value + 
                (context.kSize() - cDefaultKSize)
            );

        context.template Advance<
            TIterationMask,
            VIOffset,
            VJOffset
        >(kOffset);
    }
};

