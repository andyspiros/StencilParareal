#ifndef MATFILE_H_
#define MATFILE_H_

#include <iostream>
#include <fstream>
#include <string>
#include <stdint.h>
#include "SharedInfrastructure.h"

class MatFile
{
public:
    MatFile(const std::string& fname);

    void startCell(const std::string&, int n);
    void endCell();

    template<typename TDataField>
    void addField(const TDataField& field, int inc = -1);

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

// Implementation of template methods

template<typename TDataField>
void MatFile::addField(const TDataField& field, int inc)
{
    std::ostringstream fieldNameStream;
    if (inc >= 0)
        fieldNameStream << field.name() << "_" << inc;
    else
        fieldNameStream << field.name();
    std::string fieldName = fieldNameStream.str();

    IJKSize size = field.calculationDomain();
    IJKBoundary boundary = field.boundary();
    int iSize = size.iSize() - boundary.iMinusOffset() + boundary.iPlusOffset();
    int jSize = size.jSize() - boundary.jMinusOffset() + boundary.jPlusOffset();
    int kSize = size.kSize() - boundary.kMinusOffset() + boundary.kPlusOffset();
    int totsize = iSize * jSize * kSize;

    int iStart = boundary.iMinusOffset();
    int jStart = boundary.iMinusOffset();
    int kStart = boundary.kMinusOffset();
    int iEnd = iStart + iSize;
    int jEnd = jStart + jSize;
    int kEnd = kStart + kSize;


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
    int32_t nameheader[2] = {1, fieldName.size()};
    int npd = (fieldName.size()+7)/8*8 - fieldName.size();
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
        /* name  */ 8 + (fieldName.size()+7)/8*8 +
        /* data  */ 8 + datalength[1];

    // Start writing
    fs_.write(reinterpret_cast<const char*>(header), 8);
    fs_.write(reinterpret_cast<const char*>(flags), 16);
    fs_.write(reinterpret_cast<const char*>(sizes), 24);
    fs_.write(reinterpret_cast<const char*>(nameheader), 8);
    fs_.write(fieldName.c_str(), fieldName.size());
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
}

#endif // MATFILE_H_
