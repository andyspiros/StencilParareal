#ifndef MATFILE_H_
#define MATFILE_H_

#include <iostream>
#include <fstream>
#include <string>
#include "stdint.h"
#include "SharedInfrastructure.h"

class MatFile
{
public:
    MatFile(const std::string& fname);

    void startCell(const std::string&, int n);
    void endCell();

    template<typename TDataField>
    void addField(const TDataField& field, int inc = -1);

    template<typename TDataField>
    void addField(const std::string& fieldName, const TDataField& field, int inc = -1);

    void close();

private:
    std::ofstream fs_;
    int lengthPos_;

    void writeHeader();

};


/////////////
// Helpers //
/////////////

template<typename T>
struct MATTraits;

template<>
struct MATTraits<float>
{
    static const int size = 4;
    static const int mitype = 7;
};

template<>
struct MATTraits<double>
{
    static const int size = 8;
    static const int mitype = 9;
};



////////////////////
// Implementation //
////////////////////

template<typename TDataField>
void MatFile::addField(const std::string& fieldName, const TDataField& field, int inc)
{
#ifdef __CUDA_BACKEND__
    const bool synchronize = field.isDeviceUpToDate();
    if (synchronize)
        field.SynchronizeHostStorage();
#endif

    std::ostringstream nameStream;
    nameStream << fieldName;
    if (inc >= 0)
        nameStream << "_" << inc;
    const std::string name = nameStream.str();

    IJKSize size = field.calculationDomain();
    IJKBoundary boundary = field.boundary();

    int iSize, jSize, kSize,
        iStart, jStart, kStart,
        iEnd, jEnd, kEnd;
    bool withHalo = false;

    if (withHalo)
    {
        iSize = size.iSize() - boundary.iMinusOffset() + boundary.iPlusOffset();
        jSize = size.jSize() - boundary.jMinusOffset() + boundary.jPlusOffset();
        kSize = size.kSize() - boundary.kMinusOffset() + boundary.kPlusOffset();

        iStart = boundary.iMinusOffset();
        jStart = boundary.iMinusOffset();
        kStart = boundary.kMinusOffset();
    }
    else
    {
        iSize = size.iSize();
        jSize = size.jSize();
        kSize = size.kSize();

        iStart = 0;
        jStart = 0;
        kStart = 0;
    }

    iEnd = iStart + iSize;
    jEnd = jStart + jSize;
    kEnd = kStart + kSize;

    int totsize = iSize * jSize * kSize;

    // Prepare header
    uint32_t header[2] = {14, 0};
    uint32_t flags[4] = {6, 8, 6, 0};

    // Sizes
    int32_t sizes[6] = {
        5,
        12,
        iSize,
        jSize,
        kSize,
        0
    };

    // Name
    int32_t nameheader[2] = {1, name.size()};
    int npd = (name.size()+7)/8*8 - name.size();
    const char padding[] = "        ";

    // Data length
    int32_t datalength[2] = {
        MATTraits<Real>::mitype,
        MATTraits<Real>::size * totsize
    };

    // Total size
    header[1] =
        /* flags */ 16 +
        /* sizes */ 24 +
        /* name  */ 8 + (name.size()+7)/8*8 +
        /* data  */ 8 + datalength[1];

    // Start writing
    fs_.write(reinterpret_cast<const char*>(header), 8);
    fs_.write(reinterpret_cast<const char*>(flags), 16);
    fs_.write(reinterpret_cast<const char*>(sizes), 24);
    fs_.write(reinterpret_cast<const char*>(nameheader), 8);
    fs_.write(name.c_str(), name.size());
    fs_.write(padding, npd);
    fs_.write(reinterpret_cast<const char*>(datalength), 8);

    // Write the data
    for (int k = kStart; k < kEnd; ++k)
        for (int j = jStart; j < jEnd; ++j)
            for (int i = iStart; i < iEnd; ++i)
            {
                fs_.write(
                    reinterpret_cast<const char*>(&field(i,j,k)),
                    MATTraits<Real>::size
                    );
            }

#ifdef __CUDA_BACKEND__
    if (synchronize)
        field.SynchronizeDeviceStorage();
#endif
}

template<typename TDataField>
void MatFile::addField(const TDataField& field, int inc)
{
    addField(field.name(), field, inc);
}

#endif // MATFILE_H_

