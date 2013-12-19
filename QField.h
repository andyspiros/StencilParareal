#ifndef QFIELD_H_
#define QFIELD_H_

#include "SharedInfrastructure.h"
#include <iostream>

template<typename DataFieldType>
class QField
{
public:

    enum CurrentLaplace
    {
        main, temp, free
    };

    // Constructor
    QField(DataFieldType& q)
        : qmain_(q)
    {
        IJKSize domain = q.calculationDomain();
        IJKBoundary boundary = q.boundary();
        KBoundary kboundary;
        kboundary.Init(boundary.kMinusOffset(), boundary.kPlusOffset());

        qtemp_.Init("qtemp", domain, kboundary);
        qlaplace_.Init("qlaplace", domain, kboundary);
        currentLaplace_ = free;
    }
    
    void setLaplaceMain()
    {
        restore();
        qlaplace_.SwapWith(qmain_);
        currentLaplace_ = main;
    }

    void setLaplaceTemp()
    {
        restore();
        qlaplace_.SwapWith(qtemp_);
        currentLaplace_ = temp;
    }

    void restore()
    {
        switch (currentLaplace_)
        {
        case main:
            qlaplace_.SwapWith(qmain_);
            break;
        case temp:
            qlaplace_.SwapWith(qtemp_);
            break;
        default:
            break;
        }
        currentLaplace_ = free;
    }

    DataFieldType& qlaplace() { return qlaplace_; }
    DataFieldType& qmain() { return qmain_; }
    DataFieldType& qtemp() { return qtemp_; }


private:
    // qmain is the data field provided by the user
    // qtemp is the data field used in the inner part of the Runge-Kutta integration
    // qlapl is the input data field of the laplacian
    DataFieldType& qmain_;
    DataFieldType qtemp_;
    DataFieldType qlaplace_;

    CurrentLaplace currentLaplace_;
};

#endif // QFIELD_H_
