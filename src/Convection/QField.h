#ifndef QFIELD_H_
#define QFIELD_H_

#include "SharedInfrastructure.h"
#include <iostream>

template<typename DataFieldType>
class QField
{
public:

    enum CurrentRHS
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
        qrhs_.Init("qrhs", domain, kboundary);
        currentRHS_ = free;
    }
    
    void setRHSMain()
    {
        restore();
        qrhs_.SwapWith(qmain_);
        currentRHS_ = main;
    }

    void setRHSTemp()
    {
        restore();
        qrhs_.SwapWith(qtemp_);
        currentRHS_ = temp;
    }

    void restore()
    {
        switch (currentRHS_)
        {
        case main:
            qrhs_.SwapWith(qmain_);
            break;
        case temp:
            qrhs_.SwapWith(qtemp_);
            break;
        default:
            break;
        }
        currentRHS_ = free;
    }

    DataFieldType& qrhs() { return qrhs_; }
    DataFieldType& qmain() { return qmain_; }
    DataFieldType& qtemp() { return qtemp_; }


private:
    // qmain is the data field provided by the user
    // qtemp is the data field used in the inner part of the Runge-Kutta integration
    // qrhs is the input data field of the right hand side
    DataFieldType& qmain_;
    DataFieldType qtemp_;
    DataFieldType qrhs_;

    CurrentRHS currentRHS_;
};

#endif // QFIELD_H_
