#pragma once

#include <cassert>
#include <iostream>
#include <string>
#include <boost/type_traits/add_const.hpp>
#include <boost/mpl/bool.hpp>
#include "Definitions.h"
#include "IJKSize.h"
#include "IJKIndex.h"
#include "IJKBoundary.h"
#include "KBoundary.h"
#include "ExternalStorage.h"
#include "DataFieldStorageFormat.h"

/**
* @class DataField
* Base class of different data field implementations
*/
template<
    typename TStorage,
    typename TDataFieldImpl>
class DataField
{
    DISALLOW_COPY_AND_ASSIGN(DataField);
protected: // hide the constructors
    __CPU__
    DataField() { isInitialized_ = false; }
    __CPU__
    ~DataField() {}

    // method used by the actual data field implementation to implement copy constructor and assignment operator
    __CPU__
    void Assign(const DataField& other) 
    {
        // copy values
        isInitialized_ = other.isInitialized_;
        name_ = other.name_;
        calculationDomain_ = other.calculationDomain_;
        boundary_ = other.boundary_; 
    }

public: 
    // define the data field interface
    typedef TStorage StorageType; 
    typedef typename StorageType::StorageFormat StorageFormat;
    typedef typename StorageType::ValueType ValueType;
    typedef ValueType* PointerType;
    typedef ValueType& ReferenceType;
    typedef typename boost::add_const<ValueType>::type& ConstReferenceType;
    typedef ExternalStorage<ValueType> ExternalStorageType;    

    /**
    * Init the data field to a certain dimension
    * @param name name of the data field
    * @param calculationDomain size of the field without any halo or boundary
    * @param kBoundary boundary around the calculation domain in k direction
    */
    __CPU__
    void Init(std::string name, const IJKSize& calculationDomain, const KBoundary& kBoundary)
    {
        // initialize data field
        initDataField(name, calculationDomain, kBoundary);

        // initialize the data field implementation
        static_cast<TDataFieldImpl*>(this)->InitImpl(kBoundary);
    }

    /**
    * Init the data field to a certain dimension
    * @param name name of the data field
    * @param externalStorage external storage
    * @param calculationDomain size of the field without any halo or boundary
    * @param kBoundary boundary around the calculation domain in k direction
    */
    __CPU__
    void Init(std::string name, const ExternalStorageType& externalStorage, const IJKSize& calculationDomain, const KBoundary& kBoundary)
    {
        // initialize data field
        initDataField(name, calculationDomain, kBoundary);
      
        // initialize the data field implementation 
        static_cast<TDataFieldImpl*>(this)->InitImpl(externalStorage, kBoundary);
    }
    
    /**
    * @return true if the field is initialized
    */
    __CPU__
    bool isInitialized() const { return isInitialized_; }

    /**
    * Set data field storage to external storage
    * @param externalStorage external storage
    */
    __CPU__
    void SetExternalStorage(const ExternalStorageType& externalStorage)
    {
        // call the set external storage implementation
        static_cast<TDataFieldImpl*>(this)->SetExternalStorageImpl(externalStorage);
    }

    /**
    * Swap storage with other data field (used by SwapDataField)
    * @param other data field
    */
    __CPU__
    void SwapWith(TDataFieldImpl& other)
    {
        assert(calculationDomain_ == other.calculationDomain_);
        assert(isInitialized_ == other.isInitialized_);
        assert(boundary_ == other.boundary_); 
        static_cast<TDataFieldImpl*>(this)->SwapWithImpl(other);
    }

    /**
    * Access operator (calculation domain origin is i=0, j=0, k=0)
    * @param i position in i dimension
    * @param j position in j dimension
    * @param k position in k dimension
    * @return a reference to the value at the position
    */
    __CPU__
    ReferenceType operator() (const int i, const int j, const int k)
    {
        return static_cast<TDataFieldImpl*>(this)->AccessImpl(i,j,k);
    }

    /**
    * Access operator (calculation domain origin is i=0, j=0, k=0)
    * @param i position in i dimension
    * @param j position in j dimension
    * @param k position in k dimension
    * @return a reference to the value at the position
    */
    __CPU__
    ConstReferenceType operator() (const int i, const int j, const int k) const
    {
        return static_cast<const TDataFieldImpl*>(this)->AccessImpl(i,j,k);
    }
    
    /**
    * @return the name of the data field
    */
    __CPU__
    std::string name() const { return name_; }

    /**
    * @return size of calculation domain
    */
    __CPU__
    const IJKSize& calculationDomain() const { return calculationDomain_; }

    /**
    * @return maximal allowed boundary
    */
    __CPU__
    const IJKBoundary& boundary() const { return boundary_; }
    
    /**
    * @return data field storage
    */
    __CPU__
    const StorageType& storage() const { return static_cast<const TDataFieldImpl*>(this)->storageImpl(); }

private:
    // init data field members
    void initDataField(std::string name, const IJKSize& calculationDomain, const KBoundary& kBoundary)
    {
        isInitialized_ = true;
        name_ = name;
        calculationDomain_ = calculationDomain;

        // define the boundary type
        typedef typename StorageFormat::IJBoundary IJBoundary;

        // initialize the maximal boundary
        boundary_.Init(
            IJBoundary::IMinusOffset::value, IJBoundary::IPlusOffset::value,
            IJBoundary::JMinusOffset::value, IJBoundary::JPlusOffset::value,
            kBoundary.kMinusOffset(), kBoundary.kPlusOffset()
        );
    }

    bool isInitialized_;
    std::string name_;
    IJKSize calculationDomain_;
    IJKBoundary boundary_; 
};

/**
* @struct is_data_field
* Meta function returning true the parameter is a data field
*/
template<typename T>
struct is_data_field : boost::mpl::false_ {};
