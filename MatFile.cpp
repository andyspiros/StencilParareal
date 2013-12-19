#include "MatFile.h"
#include <stdint.h>
#include "sstream"

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

MatFile::MatFile(const std::string& fname)
{
    fs_.open(
        fname.c_str(), 
        std::ofstream::out | std::ofstream::binary | std::ofstream::trunc
    );
    writeHeader();
}

void MatFile::startCell(const std::string& fieldName, int n)
{
    int32_t header[] = {
        // Header
        14, 0,
        
        // Flags
        6, 8, 1, 0,

        // Sizes
        5, 4, n, 0,

        // Name
        1, fieldName.size()
    };

    // Name padding
    const char padding[] = "        ";
    const int npd = (fieldName.size()+7)/8*8 - fieldName.size();

    // Write to file
    fs_.write(reinterpret_cast<const char*>(header), 4);
    lengthPos_ = fs_.tellp();
    fs_.write(reinterpret_cast<const char*>(header+1), 11*4);
    fs_.write(fieldName.c_str(), fieldName.size());
    fs_.write(padding, npd);

    // Ready for writing fields...
}

void MatFile::endCell()
{
    int endPos = fs_.tellp();
    int32_t length = endPos - lengthPos_ - 4;
    fs_.seekp(lengthPos_);
    fs_.write(reinterpret_cast<const char*>(&length), 4);
    fs_.seekp(endPos);
}


void MatFile::addField(const IJKRealField& field, int inc)
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

void MatFile::writeHeader()
{
    static const char description[] =
        "Stencil-parareal project. Output data file. Good luck           "
        "                                                                ";
    const uint16_t version = 0x0100;
    const uint16_t endian = 0x4D49;

    fs_.seekp(0);
    fs_.write(description, 124);
    fs_.write(reinterpret_cast<const char*>(&version), 2);
    fs_.write(reinterpret_cast<const char*>(&endian), 2);
}
