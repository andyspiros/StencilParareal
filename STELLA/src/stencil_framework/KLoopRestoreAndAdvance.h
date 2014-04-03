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
* @struct KLoopRestoreAndAdvance
* Structure used by the k loop to restore the parameter positions and move them to the start position of the loop
*/
template<
    typename TToKPosition,
    typename TIterationMask>
struct KLoopRestoreAndAdvance
{
    BOOST_STATIC_ASSERT(is_k_position<TToKPosition>::value);
    
    // define the origin k position
    typedef KPosition<cKMinimumFlat, 0> OriginKPosition;

    template<typename TContext>
    __ACC__
    static void Do(TContext& context, const int iOffset, const int jOffset)
    {
        // select the implementation
        // if there is no work (empty iteration mask) select the empty implementation
        // otherwise select the runtime respectively the compile-time implementation
        typedef typename boost::mpl::eval_if<
            boost::mpl::empty<TIterationMask>,
            boost::mpl::integral_c<ImplementationPolicy, cEmpty>, 
            boost::mpl::if_<
                have_k_positions_run_time_length_component<OriginKPosition, TToKPosition>,
                boost::mpl::integral_c<ImplementationPolicy, cRuntime>, 
                boost::mpl::integral_c<ImplementationPolicy, cCompileTime>
            >
        >::type Implementation;
        
        // call the selected implementation
        doImpl(context, iOffset, jOffset, static_cast<Implementation*>(0));
    }

private:
    // empty implementation called if no work is needed -> e.g. advance by 0
    template<typename TContext>
    __ACC__
    static void doImpl(TContext& context, const int iOffset, const int jOffset, boost::mpl::integral_c<ImplementationPolicy, cEmpty>*) {}

    // implementation working with a compile-time offset
    template<typename TContext>
    __ACC__
    static void doImpl(TContext& context, const int iOffset, const int jOffset, boost::mpl::integral_c<ImplementationPolicy, cCompileTime>*) 
    {
        // compute the compile time offset
        const int kOffset =
            k_positions_direction<OriginKPosition, TToKPosition>::value * 
            k_positions_compile_time_length_component<OriginKPosition, TToKPosition>::value; 

        context.template RestoreAndAdvance<
            TIterationMask, 
            kOffset
        >(iOffset, jOffset);
    }

    // implementation computing a runtime offset
    template<typename TContext>
    __ACC__
    static void doImpl(TContext& context, const int iOffset, const int jOffset, boost::mpl::integral_c<ImplementationPolicy, cRuntime>*) 
    {
        // compute the run time offset 
        // --> add the deviation between default k size and runtime k size
        const int kOffset = 
            k_positions_direction<OriginKPosition, TToKPosition>::value * 
            (
                k_positions_compile_time_length_component<OriginKPosition, TToKPosition>::value + 
                (context.kSize() - cDefaultKSize)
            );

        context.template RestoreAndAdvance<
            TIterationMask
        >(iOffset, jOffset, kOffset);
    }
};
