#include "Periodicity.h"

inline double& access(double *pField,
        const int& i, const int& istride, const int& j, const int& jstride, const int& k, const int& kstride
    )
{
    return *(pField + i*istride + j*jstride + k*kstride);
}

void ApplyPeriodicI(double *pField,
        int haloim, int haloip, int dsizei, int stridei,
        int halojm, int halojp, int dsizej, int stridej,
        int halokm, int halokp, int dsizek, int stridek
    )
{
    const int jStart = halojm,
              jEnd   = halojm + dsizej;
    const int kStart = halokm,
              kEnd   = halokm + dsizek;

    for (int k = kStart; k < kEnd; ++k)
        for (int j = jStart; j < jEnd; ++j)
        {
            for (int i = 0; i < haloim; ++i)
            {
                access(pField, i         , stridei, j, stridej, k, stridek)
                    =
                access(pField, i + dsizei, stridei, j, stridej, k, stridek);
            }
            const int iEnd = haloim + haloip;
            for (int i = haloim; i < iEnd; ++i)
            {
                access(pField, i + dsizei, stridei, j, stridej, k, stridek)
                    =
                access(pField, i         , stridei, j, stridej, k, stridek);
            }
        }
}

void ApplyPeriodicJ(double *pField,
        int haloim, int haloip, int dsizei, int stridei,
        int halojm, int halojp, int dsizej, int stridej,
        int halokm, int halokp, int dsizek, int stridek
    )
{
    const int iStart = 0,
              iEnd   = haloim + dsizei + haloip;
    const int kStart = halokm,
              kEnd   = halokm + dsizek;

    for (int i = iStart; i < iEnd; ++i)
        for (int k = kStart; k < kEnd; ++k)
        {
            for (int j = 0; j < halojm; ++j)
            {
                access(pField, i, stridei, j         , stridej, k, stridek)
                    =
                access(pField, i, stridei, j + dsizej, stridej, k, stridek);
            }
            const int jEnd = halojm + halojp;
            for (int j = halojm; j < jEnd; ++j)
            {
                access(pField, i, stridei, j + dsizej, stridej, k, stridek)
                    =
                access(pField, i, stridei, j         , stridej, k, stridek);
            }
        }
}

void ApplyPeriodicK(double *pField,
        int haloim, int haloip, int dsizei, int stridei,
        int halojm, int halojp, int dsizej, int stridej,
        int halokm, int halokp, int dsizek, int stridek
    )
{
    const int iStart = 0,
              iEnd   = haloim + dsizei + haloip;
    const int jStart = 0,
              jEnd   = halojm + dsizej + halojp;

    for (int i = iStart; i < iEnd; ++i)
        for (int j = jStart; j < jEnd; ++j)
        {
            for (int k = 0; k < halokm; ++k)
            {
                access(pField, i, stridei, j, stridej, k         , stridek)
                    =                                            
                access(pField, i, stridei, j, stridej, k + dsizek, stridek);
            }
            const int kEnd = halokm + halokp;
            for (int k = halokm; k < kEnd; ++k)
            {
                access(pField, i, stridei, j, stridej, k + dsizek, stridek)
                    =                                            
                access(pField, i, stridei, j, stridej, k         , stridek);
            }
        }
}

void ApplyImpl(double *pStorageBase, const IJKSize& dsize, const IJKBoundary& boundary, int strides[3])
{
    const int dsizei = dsize.iSize();
    const int dsizej = dsize.jSize();
    const int dsizek = dsize.kSize();

    const int haloim = -boundary.iMinusOffset();
    const int halojm = -boundary.jMinusOffset();
    const int halokm = -boundary.kMinusOffset();
    const int haloip = boundary.iPlusOffset();
    const int halojp = boundary.jPlusOffset();
    const int halokp = boundary.kPlusOffset();

    const int stridei = strides[0];
    const int stridej = strides[1];
    const int stridek = strides[2];

    ApplyPeriodicI(pStorageBase,
        haloim, haloip, dsizei, stridei,
        halojm, halojp, dsizej, stridej,
        halokm, halokp, dsizek, stridek
    );

    ApplyPeriodicJ(pStorageBase,
        haloim, haloip, dsizei, stridei,
        halojm, halojp, dsizej, stridej,
        halokm, halokp, dsizek, stridek
    );

    ApplyPeriodicK(pStorageBase,
        haloim, haloip, dsizei, stridei,
        halojm, halojp, dsizej, stridej,
        halokm, halokp, dsizek, stridek
    );
 
}

