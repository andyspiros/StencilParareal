#include "Periodicity.h"
#include "mpi.h"

#include "SharedInfrastructure.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

template <typename TDataField>
void fillField(TDataField& field)
{
    IJKSize d = field.calculationDomain();
    IJKBoundary b = field.boundary();
    const int iStart = b.iMinusOffset();
    const int jStart = b.jMinusOffset();
    const int kStart = b.kMinusOffset();
    const int iEnd = d.iSize() + b.iPlusOffset();
    const int jEnd = d.jSize() + b.jPlusOffset();
    const int kEnd = d.kSize() + b.kPlusOffset();

    for (int i = iStart; i < iEnd; ++i)
        for (int j = jStart; j < jEnd; ++j)
            for (int k = kStart; k < kEnd; ++k)
            {
                int i_ = i - iStart;
                int j_ = j - jStart;
                int k_ = k - kStart;
                int value = k_ + 100*j_ + 10000*i_;
                field(i, j, k) = value;
            }
}

template <typename TDataField>
void print(const TDataField& field)
{
    IJKSize d = field.calculationDomain();
    IJKBoundary b = field.boundary();
    const int iStart = b.iMinusOffset();
    const int jStart = b.jMinusOffset();
    const int kStart = b.kMinusOffset();
    const int iEnd = d.iSize() + b.iPlusOffset();
    const int jEnd = d.jSize() + b.jPlusOffset();
    const int kEnd = d.kSize() + b.kPlusOffset();

    for (int i = iStart; i < iEnd; ++i)
    {
        for (int j = jStart; j < jEnd; ++j)
        {
            for (int k = kStart; k < kEnd; ++k)
            {
                std::cout << std::setw(6) << std::setfill('0') << field(i, j, k) << "  ";
            }
            std::cout << "\n";
        }
        std::cout << "\n\n";
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

#ifdef __CUDA_BACKEND__
    std::cout << "CUDA version\n";
#else
    std::cout << "CPU version\n";
#endif

    // Initialize field
    IJKSize domain;
    domain.Init(4, 4, 4);
    KBoundary kboundary;
    kboundary.Init(-3, 3);
    IJKRealField field;
    field.Init("field", domain, kboundary);

    std::cerr << "Field at " << field.storage().pStorageBase() << std::endl;

    fillField(field);
    field.SynchronizeDeviceStorage();
    print(field);

    Periodicity<IJKRealField> p(field);
    p.Apply();
    field.SynchronizeHostStorage();

    std::cout << "\n\nAfter periodicity:\n\n";
    print(field);

    MPI_Finalize();
}

