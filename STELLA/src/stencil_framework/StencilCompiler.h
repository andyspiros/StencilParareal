#pragma once

#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/vector/vector0.hpp>
#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/insert.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/back_inserter.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/plus.hpp>
#include <boost/mpl/next.hpp>
#include <boost/mpl/min_max.hpp>
#include "SharedInfrastructure.h"
#include "StencilConfiguration.h"
#include "Temporaries.h"
#include "StencilSweepGroupDescriptor.h"
#include "Stencil.h"
#include "StencilCompilerFunctors.h"

// include the wrapper matching the back end
#ifdef __CUDA_BACKEND__
#include "StencilCUDA.h"
#else
#include "StencilOpenMP.h"
#endif

#include <boost/mpl/set/set40.hpp> // depends on MAX_TUPLE_SIZE
#include <boost/mpl/map/map40.hpp> // depends on MAX_TUPLE_SIZE

/**
* @struct StencilCompiler
* Structure used to create / initialize a stencil object
*/
struct StencilCompiler 
{
    /**
    * Method initializing the stencil (signature without temporary field support)
    * @param stencil stencil to initialize
    * @param name name of the stencil
    * @param calculationDomain size of the calculation domain
    * @param configuration type holding general stencil information as calculation type and block size
    * @param parameterTuple tuple holding the stencil parameters
    * @param stencilSweepOrSweepGroupDescriptors type list containing all the sweeps or sweep groups which are executed consecutively
    */
    template<
        typename TStencilConfiguration,
        typename TParameterTupleElements,
        typename TStencilSweepOrSweepGroupDescriptors>
    __CPU__
    static void Build(
        Stencil& stencil,
        std::string name, 
        const IJKSize& calculationDomain, 
        TStencilConfiguration configuration, 
        Tuple<TParameterTupleElements> parameterTuple,
        TStencilSweepOrSweepGroupDescriptors stencilSweepOrSweepGroupDescriptors)
    {
        // call init with an empty temporary field list
        StencilCompiler::Build(stencil, name, calculationDomain, configuration, parameterTuple, boost::mpl::vector0<>(), stencilSweepOrSweepGroupDescriptors);
    }

    /**
    * Method initializing the stencil (signature with column buffer support
    * @param stencil stencil to initialize
    * @param name name of the stencil
    * @param calculationDomain size of the calculation domain
    * @param configuration type holding general stencil information as calculation type and block size
    * @param parameterTuple tuple holding the stencil parameters
    * @param temporaryFields type list containing the definition of all temporary variables and buffers needed by the stencil
    * @param stencilSweepOrSweepGroupDescriptors type list containing all the sweeps or sweep groups which are executed consecutively
    */
    template<
        typename TStencilConfiguration,
        typename TParameterTupleElements,
        typename TTemporaryFields,
        typename TStencilSweepOrSweepGroupDescriptors>
    __CPU__
    static void Build(
        Stencil& stencil,
        std::string name, 
        const IJKSize& calculationDomain, 
        TStencilConfiguration configuration, 
        Tuple<TParameterTupleElements> parameterTuple,
        TTemporaryFields temporaryFields,
        TStencilSweepOrSweepGroupDescriptors stencilSweepOrSweepGroupDescriptors) 
    {
        // check that a stencil sweep or sweep group vector is handed over
        BOOST_MPL_ASSERT_MSG( 
            boost::mpl::is_sequence<TStencilSweepOrSweepGroupDescriptors>::value,
            STENCIL_INIT_CALLED_WITH_INVALID_LOOPS_NO_SEQUENCE,
            (TStencilSweepOrSweepGroupDescriptors) 
        );

        // check the stencil configuration
        BOOST_MPL_ASSERT_MSG(
            is_stencil_configuration<TStencilConfiguration>::value,
            STENCIL_INIT_CALLED_WITH_INVALID_STENCIL_CONFIGURATION,
            (TStencilConfiguration) 
        );

        // check that a tuple is handed over
        BOOST_MPL_ASSERT_MSG(
            is_tuple_elements<TParameterTupleElements>::value,
            STENCIL_INIT_CALLED_WITH_INVALID_PARAMETER_TUPLE,
            (TParameterTupleElements) 
        );

        // check the temporary fields
        BOOST_MPL_ASSERT_MSG( 
            (
                boost::mpl::is_sequence<TTemporaryFields>::value &&
                boost::mpl::fold<
                    TTemporaryFields,
                    boost::mpl::true_,
                    boost::mpl::if_<
                        is_temporary_field<boost::mpl::_2>, 
                        boost::mpl::_1,
                        boost::mpl::false_
                    >
                >::type::value
            ),
            STENCIL_INIT_CALLED_WITH_INVALID_TEMPORARY_FIELD_LIST,
            (TTemporaryFields) 
        );

        // check the number of parameters and temporary fields is smaller than max tuple size
        BOOST_MPL_ASSERT_MSG( 
            (
                boost::mpl::size<TTemporaryFields>::value + 
                boost::mpl::size<typename TParameterTupleElements::ElementTypes>::value <= 
                MAX_TUPLE_SIZE
            ),
            STENCIL_INIT_CALLED_WITH_TOO_MANY_PARAMETERS,
            (TTemporaryFields, TParameterTupleElements) 
        );

        // check that the parameter indices used for the parameters and the temporary fields are unique
        // 1. create a set and insert all parameter and temporary indices
        // 2. check the size of the set is equivalent to number of parameters plus number of temporaries
        BOOST_MPL_ASSERT_MSG( 
            (
                boost::mpl::size<
                    typename boost::mpl::fold<
                        TTemporaryFields,
                        typename boost::mpl::fold<
                            typename TParameterTupleElements::ElementIndexes,
                            boost::mpl::set0<>,
                            boost::mpl::insert<boost::mpl::_1, boost::mpl::_2>
                        >::type,
                        boost::mpl::insert<
                            boost::mpl::_1, 
                            temporary_field_index<boost::mpl::_2> 
                        >
                    >::type
                >::value == 
                boost::mpl::size<TTemporaryFields>::value +
                boost::mpl::size<typename TParameterTupleElements::ElementTypes>::value
            ),
            STENCIL_INIT_PARAMETER_INDEXES_ARE_NOT_UNIQUE,
            (TTemporaryFields, TParameterTupleElements) 
        );

        // count number of sweep groups which is equivalent to the number of
        // sweeps outside a switch case statement plus the number of case statements
        // (note later on succeeding sweep groups with identical conditions are merged which is not taken into account here)
        BOOST_MPL_ASSERT_MSG( 
            (
                boost::mpl::fold<
                    TStencilSweepOrSweepGroupDescriptors,
                    boost::mpl::integral_c<int, 0>,
                    boost::mpl::if_<
                        boost::mpl::is_sequence<boost::mpl::_2>,
                        boost::mpl::plus<
                            boost::mpl::_1,
                            boost::mpl::size<boost::mpl::_2>
                        >,
                        boost::mpl::next<boost::mpl::_1>
                    >
                >::type::value <= MAX_GROUP_COUNT
            ),
            STENCIL_INIT_CALLED_WITH_INVALID_LOOPS_TOO_MANY_SWEEPS,
            (TStencilSweepOrSweepGroupDescriptors) 
        );
                 
        // expand all switch case statements stored as vectors 
        // embedded in the stencil sweep or sweep group vector
        typedef typename boost::mpl::fold<
            TStencilSweepOrSweepGroupDescriptors,
            boost::mpl::vector0<>,
            boost::mpl::if_<
                boost::mpl::is_sequence<boost::mpl::_2>,
                boost::mpl::copy<
                    boost::mpl::_2, 
                    boost::mpl::back_inserter<boost::mpl::_1> 
                >,
                boost::mpl::push_back<
                    boost::mpl::_1, 
                    boost::mpl::_2
                >
            >
        >::type ExpandedStencilSweepOrSweepGroupDescriptors;

        // meta creating a list of sweep group descriptors
        // first all sweep descriptors are converted into sweep descriptors
        // the resulting sweep group descriptors are inserted in a vector
        // (note that sweep group descriptors with identical condition are merged)
        typedef typename boost::mpl::fold<
            ExpandedStencilSweepOrSweepGroupDescriptors,
            boost::mpl::vector0<>,
            stencil_sweep_group_descriptors_add_stencil_sweep_group_descriptor<
                boost::mpl::_1, 
                boost::mpl::if_<
                    is_stencil_sweep_descriptor<boost::mpl::_2>,
                    create_stencil_sweep_group_descriptor_from_sweep<boost::mpl::_2>, 
                    boost::mpl::_2
                >
            >
        >::type StencilSweepGroupDescriptors;
        
        // check that the sweep or sweep group vector is valid
        // make sure define loops contains only sweeps or switch case statements
        BOOST_MPL_ASSERT_MSG( 
            (
                boost::mpl::fold<
                    StencilSweepGroupDescriptors,
                    boost::mpl::true_,
                    boost::mpl::if_<
                        is_stencil_sweep_group_descriptor<boost::mpl::_2>,
                        boost::mpl::_1,
                        boost::mpl::false_
                    >
                >::type::value
            ),
            STENCIL_INIT_CALLED_WITH_INVALID_LOOPS_EXPANSION_FAILED,
            (StencilSweepGroupDescriptors) 
        );

        // create a map which associates a scalar index with a parameter type
        // (note that these map holds all parameters which potentially can be used by a switch case)
        typedef typename boost::mpl::fold<
            typename TParameterTupleElements::ElementTypes,
            boost::mpl::map0<>,
            boost::mpl::if_<
                boost::is_same<
                    parameter_wrapper_intend<boost::mpl::_2>,
                    boost::mpl::integral_c<ParameterIntend, cScalar>
                >,
                boost::mpl::insert<
                    boost::mpl::_1,
                    boost::mpl::pair<
                        parameter_wrapper_index<boost::mpl::_2>,
                        parameter_wrapper_type<boost::mpl::_2>
                    >
                >,
                boost::mpl::_1
            >
        >::type IndexToScalarMap;

        // check that all sweep group parameter indices point to scalar parameters
        // make sure the scalar type matches the comparison value type
        // (note if the comparison value is void no comparison is performed and therefore no check is necessary)
        BOOST_MPL_ASSERT_MSG( 
            (
                boost::mpl::fold<
                    StencilSweepGroupDescriptors,
                    boost::mpl::true_,
                    boost::mpl::if_<
                        boost::mpl::is_void_<
                            stencil_sweep_group_descriptor_comparison_value<boost::mpl::_2> 
                        >,
                        boost::mpl::_1, 
                        boost::mpl::if_<
                            comparison_value_is_comparison_possible<
                                stencil_sweep_group_descriptor_comparison_value<boost::mpl::_2>,
                                boost::mpl::at<
                                    IndexToScalarMap, 
                                    stencil_sweep_group_descriptor_parameter_index<boost::mpl::_2> 
                                >
                            >,
                            boost::mpl::_1,
                            boost::mpl::false_                             
                        >
                    >
                >::type::value
            ),
            STENCIL_INIT_CALLED_WITH_INVALID_SWITCH_CASE_PARAMETERS,
            (StencilSweepGroupDescriptors) 
        );

        // check that no StageVariable is used in more than one stage of each sweep
        typedef typename boost::mpl::fold< // first, define list of stage variables
            TTemporaryFields,
            boost::mpl::vector0<>,
            boost::mpl::if_<
                is_stage_variable<boost::mpl::_2>,
                boost::mpl::push_back<boost::mpl::_1, temporary_field_index<boost::mpl::_2> >,
                boost::mpl::_1
            >
        >::type StageVariableIndexes;

        // check that each stage variable is used in at most one stage of each stencil sweep
        BOOST_MPL_ASSERT_MSG(
            (
                boost::mpl::fold<
                    StageVariableIndexes,
                    boost::mpl::true_,
                    boost::mpl::if_<
                        stencil_sweep_group_descriptors_is_parameter_usage_unique<
                            StencilSweepGroupDescriptors, 
                            boost::mpl::_2
                        >,
                        boost::mpl::_1,
                        boost::mpl::false_
                    >
                >::type::value
            ),
            STAGE_VARIABLE_USED_IN_MORE_THAN_ONE_STAGE_OF_A_SWEEP,
            (StageVariableIndexes)
        );
        
        // define the data field parameter index type list
        typedef typename data_field_parameters<
            typename TParameterTupleElements::ElementTypes
        >::type DataFieldParameterIndexes;

        // check the calculation domain sizes of the parameters are ok
        modify_tuple<
            CalculationDomainCheckFunctor, 
            DataFieldParameterIndexes,
            const IJKSize
        >(parameterTuple, calculationDomain);

        // compute a pair containing the ij parameter boundary size
        // and a boolean set to true if the parameter boundary is unique for all parameters
        typedef typename boost::mpl::fold<
            typename TParameterTupleElements::ElementTypes, 
            boost::mpl::pair<boost::mpl::void_, boost::mpl::true_>,
            boost::mpl::if_<
                is_data_field_parameter<boost::mpl::_2>,
                boost::mpl::pair<
                    data_field_parameter_ij_boundary<boost::mpl::_2>, 
                    boost::mpl::or_<
                        boost::mpl::is_void_<
                            boost::mpl::first<boost::mpl::_1> 
                        >,
                        boost::is_same<
                            boost::mpl::first<boost::mpl::_1>, 
                            data_field_parameter_ij_boundary<boost::mpl::_2> 
                        >
                    >
                >,
                boost::mpl::_1
            >
        >::type ParameterIJBoundaryInformation;

        // define parameter ij boundary
        typedef typename boost::mpl::first<
            ParameterIJBoundaryInformation
        >::type ParameterIJBoundary;

        // make sure there is only one parameter ij boundary size
        BOOST_MPL_ASSERT_MSG(
            boost::mpl::second<ParameterIJBoundaryInformation>::type::value, 
            STENCIL_INIT_CANNOT_DETERMINE_UNIQUE_IJ_BOUNDARY,
            (ParameterIJBoundaryInformation) 
        );

        // define the stencil wrapper type depending on the back end
        typedef typename create_stencil<
            TStencilConfiguration, 
            ParameterIJBoundary,
            TParameterTupleElements, 
            TTemporaryFields, 
            StencilSweepGroupDescriptors
        >::type StencilWrapperType;
        
        // create the stencil wrapper
        StencilWrapperType* pStencilWrapper = new StencilWrapperType();

        // depending on the data field types 
        pStencilWrapper->Init(calculationDomain, parameterTuple);

        // setup the stencil base class
        stencil.Init(name, pStencilWrapper);
    }
};
  
