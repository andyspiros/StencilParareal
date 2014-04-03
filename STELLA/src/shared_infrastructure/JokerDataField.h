#pragma once

#include <cassert>
#include <boost/static_assert.hpp>
#include <boost/mpl/bool.hpp>
#include "Definitions.h"
#include "DataField.h"

/**
* @class JokerDataField
* Data container for one data field storage
*/
template<typename TDataField>
class JokerDataField
{
    DISALLOW_COPY_AND_ASSIGN(JokerDataField);
    BOOST_STATIC_ASSERT(is_data_field<TDataField>::value);
public:
    typedef TDataField DataField;
    __CPU__
    JokerDataField() { pDataField_ = NULL; }
    __CPU__
    ~JokerDataField() {}
    
    /**
    * Init the data field to a certain dimension
    * @param name name of the data field
    * @param dataField initial data field (might be a dummy field)
    */
    __CPU__
    void Init(std::string name, TDataField& dataField)
    {
        name_ = name;
        pDataField_ = &dataField;
    }
   
    /**
    * Sets the data field
    * @param dataField data field the joker field points to
    */
    __CPU__
    void set_dataField(TDataField& dataField)
    {
        assert(pDataField_);
        assert(pDataField_->calculationDomain() == dataField.calculationDomain());
        assert(pDataField_->boundary() == dataField.boundary());
        pDataField_ = &dataField;
    }

    /**
    * @return actual data field
    */
    __CPU__
    TDataField& dataField()
    {
        assert(pDataField_);
        return *pDataField_;
    }
    
    /**
    * @return actual data field
    */
    __CPU__
    const TDataField& dataField() const
    {
        assert(pDataField_);
        return *pDataField_;
    }

    /**
    * @return the name of the data field
    */
    __CPU__
    std::string name() const { return name_; }

private:
    std::string name_;
    TDataField* pDataField_;
};

/**
* @struct is_joker_data_field
* Meta function returning true the parameter is a joker field
*/
template<typename T>
struct is_joker_data_field : boost::mpl::false_ {};

template<typename TDataField>
struct is_joker_data_field<JokerDataField<TDataField> > : boost::mpl::true_ {};

