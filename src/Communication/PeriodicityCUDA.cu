#include "Periodicity.h"
#include <cstdio>

__device__ double& access(double *pField,
        int i, int istride, int j, int jstride, int k, int kstride
    )
{
    return *(pField + i*istride + j*jstride + k*kstride);
}

__global__ void ApplyPeriodicI(double *pField,
        int haloim, int haloip, int dsizei, int stridei,
        int halojm, int halojp, int dsizej, int stridej,
        int halokm, int halokp, int dsizek, int stridek
    )
{
    int tj = blockIdx.x*blockDim.x + threadIdx.x;
    int tk = blockIdx.y*blockDim.y + threadIdx.y;

    if (tj >= dsizej || tk >= dsizek)
        return;

    //if (tj == 0 && tk == 0 && blockIdx.z == 0)
    //{
    //    std::printf("Value at (0,0,0): %6d  |  %6d\n", static_cast<int>(access(pField, 0, stridei, 0, stridej, 0, stridek)), static_cast<int>(*(pField+  0)));
    //    std::printf("Value at (1,0,0): %6d  |  %6d\n", static_cast<int>(access(pField, 1, stridei, 0, stridej, 0, stridek)), static_cast<int>(*(pField+  1)));
    //    std::printf("Value at (0,1,0): %6d  |  %6d\n", static_cast<int>(access(pField, 0, stridei, 1, stridej, 0, stridek)), static_cast<int>(*(pField+ 32)));
    //    std::printf("Value at (0,0,1): %6d  |  %6d\n", static_cast<int>(access(pField, 0, stridei, 0, stridej, 1, stridek)), static_cast<int>(*(pField+320)));
    //    std::printf("Strides: %d  %d  %d\n", stridei, stridej, stridek);

    //    std::printf("First 3200 values:\n");
    //    for (int i = 0; i < 3200; ++i)
    //        std::printf("%3d     %6d\n", i, static_cast<int>(*(pField+i)));
    //}


    tj += halojm;
    tk += halokm;

    if (blockIdx.z == 0)
    {
        for (int i = 0; i < haloim; ++i)
        {
            //double *dest = &(access(pField, i,        stridei, tj, stridej, tk, stridek));
            //double *orig = &(access(pField, i+dsizei, stridei, tj, stridej, tk, stridek));
            //std::printf("Setting %p (%2d, %2d, %2d) to value of %p (%2d, %2d, %2d)\n",
            //        dest, i, tj, tk, orig, i+dsizei, tj, tk);
            access(pField, i,        stridei, tj, stridej, tk, stridek)
                =
            access(pField, i+dsizei, stridei, tj, stridej, tk, stridek);
        }
    }
    else if (blockIdx.z == 1)
    {
        const int istop = haloim + haloip;
        for (int i = haloim; i < istop; ++i)
        {
            access(pField, i+dsizei, stridei, tj, stridej, tk, stridek)
                =
            access(pField, i,        stridei, tj, stridej, tk, stridek);
        }
    }
}

__global__ void ApplyPeriodicJ(double *pField,
        int haloim, int haloip, int dsizei, int stridei,
        int halojm, int halojp, int dsizej, int stridej,
        int halokm, int halokp, int dsizek, int stridek
    )
{
    int ti = blockIdx.x*blockDim.x + threadIdx.x;
    int tk = blockIdx.y*blockDim.y + threadIdx.y;

    const int tsizei = haloim + haloip + dsizei;

    if (ti >= tsizei || tk >= dsizek)
        return;

    tk += halokm;

    if (blockIdx.z == 0)
    {
        for (int j = 0; j < halojm; ++j)
        {
            access(pField, ti, stridei, j       , stridej, tk, stridek)
                =                                
            access(pField, ti, stridei, j+dsizej, stridej, tk, stridek);
        }
    }
    else if (blockIdx.z == 1)
    {
        const int jstop = halojm + halojp;
        for (int j = halojm; j < jstop; ++j)
        {
            access(pField, ti, stridei, j+dsizej, stridej, tk, stridek)
                =                                
            access(pField, ti, stridei, j       , stridej, tk, stridek);
        }
    }
}

__global__ void ApplyPeriodicK(double *pField,
        int haloim, int haloip, int dsizei, int stridei,
        int halojm, int halojp, int dsizej, int stridej,
        int halokm, int halokp, int dsizek, int stridek
    )
{
    int ti = blockIdx.x*blockDim.x + threadIdx.x;
    int tj = blockIdx.y*blockDim.y + threadIdx.y;

    const int tsizei = haloim + haloip + dsizei;
    const int tsizej = halojm + halojp + dsizej;

    if (ti >= tsizei || tj >= tsizej)
        return;

    if (blockIdx.z == 0)
    {
        for (int k = 0; k < halokm; ++k)
        {
            access(pField, ti, stridei, tj, stridej, k       , stridek)
                =                                             
            access(pField, ti, stridei, tj, stridej, k+dsizek, stridek);
        }
    }
    else if (blockIdx.z == 1)
    {
        const int kstop = halokm + halokp;
        for (int k = halokm; k < kstop; ++k)
        {
            access(pField, ti, stridei, tj, stridej, k+dsizek, stridek)
                =                                             
            access(pField, ti, stridei, tj, stridej, k       , stridek);
        }
    }

    //if (ti == 0 && tj == 0 && blockIdx.z == 0)
    //{
    //    std::printf("First 3200 values after update:\n");
    //    for (int i = 0; i < 3200; ++i)
    //        std::printf("%3d     %6d\n", i, static_cast<int>(*(pField+i)));
    //}
}

__host__ static inline void ApplyPeriodic(double *pField,
        const IJKSize& dsize, const IJKBoundary& boundary, int strides[3]
    )
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

    const int tsizei = dsizei + haloim + haloip;
    const int tsizej = dsizej + halojm + halojp;
    //const int tsizek = dsizek + halokm + halokp;

    const int stridei = strides[0];
    const int stridej = strides[1];
    const int stridek = strides[2];

    // 32 x 4 blocksize
    const int blockSizex = 32;
    const int blockSizey = 4;
    int iblocks, jblocks, kblocks;

    // Periodicity in I
    jblocks = dsizej / blockSizex + !!(dsizej % blockSizex);
    kblocks = dsizek / blockSizey + !!(dsizek % blockSizey);
    //std::cout << "Dsize: " << dsizei << " " << dsizej << " " << dsizek << std::endl;
    //std::cout << "PeriodicI with " << jblocks << "x" << kblocks << " blocks\n";
    ApplyPeriodicI<<< dim3(jblocks, kblocks, 2), dim3(blockSizex, blockSizey) >>>(
            pField,
            haloim, haloip, dsizei, stridei,
            halojm, halojp, dsizej, stridej,
            halokm, halokp, dsizek, stridek
        );

    // Periodicity in J
    iblocks = tsizei / blockSizex + !!(tsizei % blockSizex);
    kblocks = dsizek / blockSizey + !!(dsizek % blockSizey);
    //std::cout << "PeriodicJ with " << iblocks << "x" << kblocks << " blocks\n";
    ApplyPeriodicJ<<< dim3(iblocks, kblocks, 2), dim3(blockSizex, blockSizey) >>>(
            pField,
            haloim, haloip, dsizei, stridei,
            halojm, halojp, dsizej, stridej,
            halokm, halokp, dsizek, stridek
        );

    // Periodicity in K
    iblocks = tsizei / blockSizex + !!(tsizei % blockSizex);
    jblocks = tsizej / blockSizey + !!(tsizej % blockSizey);
    //std::cout << "PeriodicK with " << iblocks << "x" << jblocks << " blocks\n";
    ApplyPeriodicK<<< dim3(iblocks, jblocks, 2), dim3(blockSizex, blockSizey) >>>(
            pField,
            haloim, haloip, dsizei, stridei,
            halojm, halojp, dsizej, stridej,
            halokm, halokp, dsizek, stridek
        );
}


void ApplyImpl(double *pStorageBase, const IJKSize& dsize, const IJKBoundary& boundary, int strides[3])
{
    ApplyPeriodic(pStorageBase, dsize, boundary, strides);
}

