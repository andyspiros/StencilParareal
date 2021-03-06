#ifndef CONVECTIONSOLUTION_H_
#define CONVECTIONSOLUTION_H_

#include <cmath>
#include "SharedInfrastructure.h"


inline double exactQ(double nu0, double nufreq, double cx, double cy, double cz,
                     double x, double y, double z, double t)
{
    const double pi = 3.14159265358979;
    const double a = nufreq > 0 ? nu0 * t - (nu0 / (2.*nufreq)) * (std::cos(nufreq*t)-1.) : nu0*t;
    const double res = sin(2.*pi*(x-cx*t))*sin(2.*pi*(y-cy*t))*sin(2.*pi*(z-cz*t)) * exp(-12.*pi*pi*a);

#ifndef NDEBUG
    std::cout << "Exact at t = " << t << ", (" << x << ", " << y << ", " << z << ") = " << res << std::endl;
#endif

    return res;
}


template<typename TDataField>
void fillQ(TDataField& q,
           double nu0, double nufreq, double cx, double cy, double cz,
           double t,
           double xstart, double xend,
           double ystart, double yend,
           double zstart, double zend
          )
{
    IJKSize domain = q.calculationDomain();
    const int iSize = domain.iSize();
    const int jSize = domain.jSize();
    const int kSize = domain.kSize();

    const double dxhalf = (xend-xstart)/iSize / 2.;
    const double dyhalf = (yend-ystart)/jSize / 2.;
    const double dzhalf = (zend-zstart)/kSize / 2.;

    double x, y, z;

    for (int i = 0; i < iSize; ++i)
        for (int j = 0; j < jSize; ++j)
            for (int k = 0; k < kSize; ++k)
            {
                x = xstart + (2*i+1)*dxhalf;
                y = ystart + (2*j+1)*dyhalf;
                z = zstart + (2*k+1)*dzhalf;
                q(i, j, k) = exactQ(nu0, nufreq, cx, cy, cz, x, y, z, t);
            }
}


template<typename TDataField>
double computeError(const TDataField& q,
                double nu0, double nufreq, double cx, double cy, double cz,
                double t,
                double xstart, double xend,
                double ystart, double yend,
                double zstart, double zend,
                double& exact_inf, double &error_inf,
                TDataField* errfield=0
            )
{
    IJKSize domain = q.calculationDomain();
    const int iSize = domain.iSize();
    const int jSize = domain.jSize();
    const int kSize = domain.kSize();

    const double dxhalf = (xend-xstart)/iSize / 2.;
    const double dyhalf = (yend-ystart)/jSize / 2.;
    const double dzhalf = (zend-zstart)/kSize / 2.;

    double x, y, z;

    error_inf = 0.;
    exact_inf = 0.;
    double exact, e;

    for (int i = 0; i < iSize; ++i)
        for (int j = 0; j < jSize; ++j)
            for (int k = 0; k < kSize; ++k)
            {
                // Coordinates
                x = xstart + (2*i+1)*dxhalf;
                y = ystart + (2*j+1)*dyhalf;
                z = zstart + (2*k+1)*dzhalf;

                // Exact solution
                exact = exactQ(nu0, nufreq, cx, cy, cz, x, y, z, t);
                exact_inf = std::max(std::abs(exact), exact_inf);

                // Error
                e = q(i,j,k) - exact;
                error_inf = std::max(std::abs(e), error_inf);

                // Error field
                if (errfield)
                    (*errfield)(i, j, k) = e;
            }
    return error_inf / exact_inf;
}

#endif // CONVECTIONSOLUTION_H_

