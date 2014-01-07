#ifndef KFIELD_H_
#define KFIELD_H_

#include "SharedInfrastructure.h"
#include <iostream>

template<typename DataFieldType>
class KField
{
public:
    KField(const DataFieldType& data)
    {
        IJKSize domain = data.calculationDomain();
        IJKBoundary boundary = data.boundary();
        KBoundary kboundary;
        kboundary.Init(boundary.kMinusOffset(), boundary.kPlusOffset());

        k_ .Init("k" , domain, kboundary);
        k1_.Init("k1", domain, kboundary);
        k2_.Init("k2", domain, kboundary);
        k3_.Init("k3", domain, kboundary);
        k4_.Init("k4", domain, kboundary);
        current_ = 0;
    }

    DataFieldType& k() { return k_; }
    DataFieldType& k1() { return k1_; }
    DataFieldType& k2() { return k2_; }
    DataFieldType& k3() { return k3_; }
    DataFieldType& k4() { return k4_; }

    void set(int id)
    {
        restore();
        switch (id)
        {
        case 1:
            k_.SwapWith(k1_);
            break;
        case 2:
            k_.SwapWith(k2_);
            break;
        case 3:
            k_.SwapWith(k3_);
            break;
        case 4:
            k_.SwapWith(k4_);
            break;
        default:
            std::cerr << "Invalid id for KField: " << id;
        }
        current_ = id;
    }

    void restore()
    {
        switch (current_)
        {
        case 1:
            k_.SwapWith(k1_);
            break;
        case 2:
            k_.SwapWith(k2_);
            break;
        case 3:
            k_.SwapWith(k3_);
            break;
        case 4:
            k_.SwapWith(k4_);
            break;
        default:
            break;
        }
    }

private:
    DataFieldType k_, k1_, k2_, k3_, k4_;

    int current_;
};


#endif  // KFIELD_H_

