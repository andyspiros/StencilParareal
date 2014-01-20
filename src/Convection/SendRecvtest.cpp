#include "GCL.h"
#include "utils/layout_map.h"
#include "utils/boollist.h"
#include "L2/include/halo_exchange.h"
#include "L2/include/descriptors.h"
#include "L3/include/proc_grids_3D.h"
#include "L3/include/Halo_Exchange_3D.h"

#include "SharedInfrastructure.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

template <typename TDataField>
void fillField(TDataField& field, double value)
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
                field(i, j, k) = value;
}

template <typename TDataField>
void print(std::ofstream& fs, const TDataField& field)
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
                fs << std::setw(4) << field(i, j, k) << "  ";
            }
            fs << "\n";
        }
        fs << "\n\n";
    }
}

template<typename TDataField>
void getStrides(const TDataField& field, int strides[3])
{
    typedef typename TDataField::StorageFormat::StorageOrder StorageOrder;
    const IJKSize& psize = field.storage().paddedSize();

    int psizes[3] = { psize.iSize(), psize. jSize(), psize.kSize() };

    int tmpsize = 1;

    int c;

    c = boost::mpl::at<StorageOrder, boost::mpl::int_<0> >::type::value;
    strides[c] = tmpsize;
    tmpsize *= psizes[c];

    c = boost::mpl::at<StorageOrder, boost::mpl::int_<1> >::type::value;
    strides[c] = tmpsize;
    tmpsize *= psizes[c];

    c = boost::mpl::at<StorageOrder, boost::mpl::int_<2> >::type::value;
    strides[c] = tmpsize;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int pid, commsize;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    const bool isRoot = !pid;
    const bool isLast = pid == commsize - 1;

    if(isRoot)
    {
#ifdef __CUDA_BACKEND__
        std::cout << "CUDA version\n";
#else
        std::cout << "CPU version\n";
#endif
    }


    // Open my file
    std::ostringstream fname;
    fname << "Output_" << pid << ".log";
    std::ofstream fs(fname.str().c_str());

    // Create my data field
    typedef IJKRealField DataFieldType;
    IJKSize domain; domain.Init(4, 4, 4);
    KBoundary kboundary; kboundary.Init(-3, 3);
    DataFieldType myfield;
    myfield.Init("myfield", domain, kboundary);

    // Padded size
    IJKSize paddedSize = myfield.storage().paddedSize();
    const int pSize = paddedSize.iSize() * paddedSize.jSize() * paddedSize.kSize();

    if (isRoot)
    {
        std::cout << "Padded size: "
            << paddedSize.iSize() << "x"
            << paddedSize.jSize() << "x"
            << paddedSize.kSize() << "  =  " << pSize;
    }


    // Fill field
    fillField(myfield, pid);
    for (int i = 0; i < domain.iSize(); ++i)
        for (int j = 0; j < domain.jSize(); ++j)
            for (int k = 0; k < domain.kSize(); ++k)
            {
                myfield(i, j, k) = 1000*pid + 100*i + 10*j + k;
            }
#ifdef __CUDA_BACKEND__
    myfield.SynchronizeDeviceStorage();
#endif
    fs << "My field:\n";
    print(fs, myfield);
    fs << "\nEnd of my field...\n\n" << std::endl;

    // Prepare receive data field
    IJKRealField recvField;
    recvField.Init("recv", domain, kboundary);

    // Do asynchrnous receive
    MPI_Request req[2];
    if (!isRoot)
        MPI_Irecv(recvField.storage().pStorageBase(), pSize, MPI_DOUBLE, pid-1, 0, MPI_COMM_WORLD, req);

    // Send
    if (!isLast)
        MPI_Isend(myfield.storage().pStorageBase(), pSize,  MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD, req+1);

    // Wait
    MPI_Request* rreq = req;
    MPI_Status status[2];
    int rreqcount = 2;
    if (isRoot) { rreq=req+1; rreqcount = 1; }
    if (isLast) { rreq = req; rreqcount = 1; }
    MPI_Waitall(rreqcount, rreq, status);

    // Print recvfield
    print(fs, recvField);

    MPI_Finalize();
}
