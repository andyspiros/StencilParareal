#ifndef PERIODICITY_H_
#define PERIODICITY_H_

#include "SharedInfrastructure.h"

template<typename TDataField>
class Periodicity
{
    typedef JokerDataField<TDataField> JokerType;

public:
    Periodicity(TDataField& field)
    {
        joker_.Init(field.name(), field);
        pJoker_ = &joker_;
        Init(field);
    }

    Periodicity(JokerType& joker)
    {
        pJoker_ = &joker;
        const TDataField& field = joker.dataField();
        Init(field);
    }

    Periodicity(const Periodicity& other)
    {
        if (other.pJoker_ == &other.joker_)
        {
            // Other owns the joker field: take its field
            joker_.Init(other.joker_.dataField());
            pJoker_ = &joker_;
        }
        else
        {
            // Other does not own the joker data field: take the joker
            pJoker_ = other.pJoker_;
        }
        Init(pJoker_->dataField());
    }

    Periodicity& operator=(const Periodicity& other)
    {
        if (other.pJoker_ == &other.joker_)
        {
            // Other owns the joker field: take its field
            joker_.Init(other.joker_.dataField());
            pJoker_ = &joker_;
        }
        else
        {
            // Other does not own the joker data field: take the joker
            pJoker_ = other.pJoker_;
        }
        Init(pJoker_->dataField());

        return *this;
    }

    void Apply()
    {
        //std::cout << "Applying periodicity on "
        //    << pJoker_->dataField().storage().pStorageBase()
        //    << std::endl;
        ApplyImpl(pJoker_->dataField().storage().pStorageBase(), dsize_, boundary_, strides_);
    }


private:
    void Init(const TDataField& field)
    {
        dsize_ = field.calculationDomain();
        boundary_ = field.boundary();
        IJKSize psize_ = field.storage().paddedSize();

        // Compute strides
        typedef typename TDataField::StorageFormat::StorageOrder StorageOrder;

        int psizes_[3] = { psize_.iSize(), psize_. jSize(), psize_.kSize() };
        int tmpsize = 1;
        int c;

        c = boost::mpl::at<StorageOrder, boost::mpl::int_<0> >::type::value;
        strides_[c] = tmpsize;
        tmpsize *= psizes_[c];

        c = boost::mpl::at<StorageOrder, boost::mpl::int_<1> >::type::value;
        strides_[c] = tmpsize;
        tmpsize *= psizes_[c];

        c = boost::mpl::at<StorageOrder, boost::mpl::int_<2> >::type::value;
        strides_[c] = tmpsize;
    }

    JokerType joker_;
    JokerType *pJoker_;

    IJKSize dsize_;
    IJKBoundary boundary_;
    int strides_[3];
};

void ApplyImpl(double *pStorageBase, const IJKSize& dsize, const IJKBoundary& boundary, int strides[3]);

#endif // PERIODICITY_H

