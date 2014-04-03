#pragma once

#include <cassert>
#include <iostream>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/if.hpp>
#include "SharedInfrastructure.h"
#include "StencilStageEnvironment.h"
#include "StencilStage.h"

/**
* @struct is_data_field_parameter
* Meta function returning true if the parameter is a data field parameter
*/
template<typename TWrappedParameter>
struct is_data_field_parameter : 
    boost::mpl::or_<
        is_data_field<typename parameter_wrapper_type<TWrappedParameter>::type>,
        is_joker_data_field<typename parameter_wrapper_type<TWrappedParameter>::type>
    >
{};

/**
* @struct data_field_parameters
* Meta function returning an index type list containing all data field and joker data field parameters in TElementTypes
*/
template<typename TParameterElementTypes>
struct data_field_parameters
{
    typedef typename boost::mpl::fold<
        TParameterElementTypes,
        boost::mpl::vector0<>,
        boost::mpl::if_<
            is_data_field_parameter<boost::mpl::_2>,
            boost::mpl::push_back<boost::mpl::_1, parameter_wrapper_index<boost::mpl::_2> >, 
            boost::mpl::_1 
        >
    >::type type;
};

/**
* @struct data_field_parameter_ij_boundary
* Meta function returning ij boundary of a data field parameter
*/
template<typename T>
struct data_field_parameter_ij_boundary;

// return data field ij boundary
template<
    typename TDataField,
    ParameterIntend VParameterIntend,
    int VParameterIndex>
struct data_field_parameter_ij_boundary<ParameterWrapper<TDataField, VParameterIntend, VParameterIndex> >
{
    typedef typename TDataField::StorageFormat::IJBoundary type;
};

// return joker data field ij boundary
template<
    typename TDataField,
    ParameterIntend VParameterIntend,
    int VParameterIndex>
struct data_field_parameter_ij_boundary<ParameterWrapper<JokerDataField<TDataField>, VParameterIntend, VParameterIndex> >
{
    typedef typename TDataField::StorageFormat::IJBoundary type;
};

/**
* @struct CalculationDomainCheckFunctor
* Functor checking the calculation domain size of the data fields handed over to the stencil.
*/
struct CalculationDomainCheckFunctor
{
    /**
    * Check the calculation domain size of data fields.
    */
    template<
        typename TDataField,
        ParameterIntend VParameterIntend,
        int VParameterIndex>
    static void Do(
        ParameterWrapper<TDataField, VParameterIntend, VParameterIndex>& wrappedField, 
        parameter_type<const IJKSize>::type calculationDomain)
    {
        sizeCheck(wrappedField.Unwrap().name(), wrappedField.Unwrap().calculationDomain(), calculationDomain);
    }

    /**
    * Check the calculation domain size of joker data fields.
    */
    template<
        typename TDataField,
        ParameterIntend VParameterIntend,
        int VParameterIndex>
    static void Do(
        ParameterWrapper<JokerDataField<TDataField>, VParameterIntend, VParameterIndex>& wrappedField, 
        parameter_type<const IJKSize>::type calculationDomain)
    {
        sizeCheck(wrappedField.Unwrap().name(), wrappedField.Unwrap().dataField().calculationDomain(), calculationDomain);
    }

private:
    static void sizeCheck(std::string name, const IJKSize& size, parameter_type<const IJKSize>::type calculationDomain)
    {
        if(size.iSize() != calculationDomain.iSize())
        {
            std::cerr << "CalculationDomainCheckFunctor -> i size doesn't match for field " << name << std::endl;
            assert(false); 
            exit(-1);
        }
        if(size.jSize() != calculationDomain.jSize())
        {       
            std::cerr << "CalculationDomainCheckFunctor -> j size doesn't match for field " << name << std::endl;
            assert(false);
            exit(-1);
        }
        if(size.kSize() != calculationDomain.kSize())
        {
            std::cerr  << "CalculationDomainCheckFunctor -> k sizes doesn't match for field " << name << std::endl;
            assert(false);
            exit(-1);
        }
    }
};
  
